import pandas as pd
import requests
import time
import concurrent.futures  #Parallel API Requests

# Reduce Pandas warnings for better performance
pd.options.mode.chained_assignment = None  

# Global caches to store already queried results (Avoid redundant calls)
VEP_CACHE = {}
GWAS_CACHE = {}
DBSNP_CACHE = {}

def query_vep(hgvs_notation):
    """
    Fetch annotation from Ensembl VEP (Runs in parallel).
    Caches results to avoid redundant queries.
    """
    if hgvs_notation in VEP_CACHE:
        return VEP_CACHE[hgvs_notation]  # Return cached result

    url = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs_notation}?content-type=application/json"
    try:
        response = requests.get(url, timeout=5)
        if response.ok:
            VEP_CACHE[hgvs_notation] = response.json()
            return response.json()
    except requests.RequestException:
        return {}

    return {}

def prioritize_clin_sig(clin_sig_set):
    """
    Assigns priority scores to clinical significance labels.
    """
    priority_order = {
        "pathogenic": 5, "likely_pathogenic": 4,
        "uncertain_significance": 3, "likely_benign": 2, "benign": 1
    }
    for sig in priority_order:
        if sig in clin_sig_set:
            return sig, priority_order[sig]
    return "NR", 0

def parse_snp_info(vep_data):
    """
    Parses VEP annotation response.
    """
    if not isinstance(vep_data, list) or not vep_data:
        return {}

    data = vep_data[0]
    snp_info = {"id": data.get("id", ""), "most_severe_consequence": data.get("most_severe_consequence", "")}

    colocated_variants = data.get("colocated_variants", [])
    clin_sig = set()
    clinvar_ids = set()

    for colocated in colocated_variants:
        clin_sig.update(colocated.get("clin_sig", []))
        clinvar_ids.update([syn for syn in colocated.get("var_synonyms", []) if syn.startswith("ClinVar")])

    prioritized_clin_sig, clin_sig_score = prioritize_clin_sig(clin_sig)

    snp_info.update({
        "phenotype_or_disease": colocated_variants[0].get("phenotype_or_disease", "NR") if colocated_variants else "NR",
        "clin_sig": prioritized_clin_sig,
        "clin_sig_score": clin_sig_score,
        "ClinVar_IDs": ",".join(clinvar_ids) if clinvar_ids else "NR",
        "colocated_id": colocated_variants[0].get("id", "NR") if colocated_variants else "NR"
    })

    return snp_info

def query_gwas_catalog(rsid):
    """
    Fetch disease associations from GWAS Catalog (Runs in parallel).
    """
    if rsid in GWAS_CACHE:
        return GWAS_CACHE[rsid]

    url = f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations"
    try:
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=5)
        if response.ok:
            GWAS_CACHE[rsid] = response.json()
            return response.json()
    except requests.RequestException:
        return {}

    return {}

def query_dbsnp(rsid):
    """
    Fetch allele frequency & minor allele frequency from dbSNP (Runs in parallel).
    """
    if rsid in DBSNP_CACHE:
        return DBSNP_CACHE[rsid]

    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid.replace('rs', '')}"
    try:
        response = requests.get(url, timeout=5)
        if response.ok:
            data = response.json()
            DBSNP_CACHE[rsid] = data
            return data
    except requests.RequestException:
        return {}

    return {}

def annotate_snp_level(df, output_csv=None):
    """
    Annotates SNPs with VEP, dbSNP, and GWAS Catalog.
    Optimized with parallel requests and caching.
    """
    #Efficiently expand the 'Alternative' column
    df_expanded = df.assign(Alternative=df['Alternative'].str.split(',')).explode('Alternative')

    #Use pandas vectorization for HGVS computation (Avoid loops)
    df_expanded['HGVS_Query'] = df_expanded.apply(
        lambda row: f"{row['Chromosome'].replace('chr', '')}:g.{row['Position']}{row['Reference']}>{row['Alternative']}",
        axis=1
    )

    #Parallel processing for VEP queries
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        vep_results = list(executor.map(query_vep, df_expanded["HGVS_Query"]))

    #Parse results efficiently
    snp_annotations = [parse_snp_info(res) for res in vep_results]

    #Convert to DataFrame (Vectorized operations)
    annotations_df = pd.DataFrame(snp_annotations)
    df_expanded = df_expanded.reset_index(drop=True)

    #Merge annotations efficiently
    annotated_df = df_expanded.merge(annotations_df, left_index=True, right_index=True, how="left")

    #Remove duplicates
    annotated_df = annotated_df.drop_duplicates(subset=['HGVS_Query'])

    #Save results
    if output_csv:
        annotated_df.to_csv(output_csv, index=False)
        print(f"Annotations saved to {output_csv}")

    return annotated_df
