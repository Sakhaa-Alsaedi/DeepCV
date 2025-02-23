import pandas as pd
import requests
import concurrent.futures  #Faster parallel API requests

#Global cache for API responses (to avoid duplicate calls)
GENE_LENGTH_CACHE = {}
OPEN_TARGETS_CACHE = {}

def preprocess_data(df):
    """Preprocess input DataFrame: convert scores to numeric and compute initial risk scores."""
    df['clin_sig_score'] = pd.to_numeric(df['clin_sig_score'], errors='coerce').fillna(0)
    df['phenotype_or_disease'] = pd.to_numeric(df['phenotype_or_disease'], errors='coerce').fillna(0)

    impact_scores = {'LOW': 1, 'MODIFIER': 2, 'MODERATE': 3, 'HIGH': 4}
    df['Impact_Score'] = df['Impact'].map(impact_scores).fillna(0)

    df['Risk_Score'] = df['clin_sig_score'] + df['phenotype_or_disease']

    return df.drop_duplicates(subset=['Chromosome', 'Position', 'Reference', 'Alternative'])

def calculate_gene_scores(df_unique):
    """Calculate weighted gene scores based on SNP risk scores."""
    gene_group = df_unique.groupby('Gene').agg(
        Num_SNPs=('ID', 'size'),
        Total_Risk_Score=('Risk_Score', 'sum')
    ).reset_index()

    gene_group['Weighted_Gene_Score'] = gene_group['Num_SNPs'] * gene_group['Total_Risk_Score']

    if 'Gene_ID' in df_unique.columns:
        gene_ids = df_unique[['Gene', 'Gene_ID']].drop_duplicates()
        gene_group = gene_group.merge(gene_ids, on='Gene', how='left')

    return gene_group.sort_values(by='Weighted_Gene_Score', ascending=False).reset_index(drop=True)

def query_gene_length(gene_id):
    """Fetch gene length from Ensembl API (optimized with caching)."""
    if gene_id in GENE_LENGTH_CACHE:
        return GENE_LENGTH_CACHE[gene_id]

    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    try:
        response = requests.get(url, headers=headers, timeout=5)
        if response.ok:
            gene_data = response.json()
            gene_length = sum(
                exon["end"] - exon["start"] + 1
                for transcript in gene_data.get("Transcript", [])
                for exon in transcript.get("Exon", [])
            )
            GENE_LENGTH_CACHE[gene_id] = gene_length if gene_length > 0 else None
            return GENE_LENGTH_CACHE[gene_id]
    except requests.RequestException:
        return None
    return None

def query_open_targets(gene_id):
    """Fetch disease associations from Open Targets API (optimized with caching)."""
    if gene_id in OPEN_TARGETS_CACHE:
        return OPEN_TARGETS_CACHE[gene_id]

    url = "https://api.platform.opentargets.org/api/v4/graphql"
    query = """
    query getAssociations($geneId: String!) {
      target(ensemblId: $geneId) {
        associatedDiseases {
          rows {
            disease { id name }
            score
          }
        }
      }
    }
    """
    try:
        response = requests.post(url, json={'query': query, 'variables': {"geneId": gene_id}}, timeout=5)
        if response.ok:
            data = response.json()
            rows = data.get('data', {}).get('target', {}).get('associatedDiseases', {}).get('rows', [])
            if rows:
                df = pd.DataFrame(rows)
                df["Gene_ID"] = gene_id
                df["disease_id"] = df["disease"].apply(lambda x: x['id'] if isinstance(x, dict) else None)
                df["disease_name"] = df["disease"].apply(lambda x: x['name'] if isinstance(x, dict) else None)
                df.drop(columns=["disease"], inplace=True)
                OPEN_TARGETS_CACHE[gene_id] = df
                return df
    except requests.RequestException:
        return None
    return None

def integrate_gene_information(df_unique, gene_group_sorted):
    """Fetch gene lengths in parallel and merge with gene scores."""
    gene_info = df_unique[['Gene', 'Gene_ID', 'Feature_Type']].drop_duplicates()

    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        gene_lengths = list(executor.map(query_gene_length, gene_info["Gene_ID"]))

    gene_info["Gene_Length_bp"] = gene_lengths

    return gene_group_sorted.merge(gene_info, on=['Gene', 'Gene_ID'], how='left')

def calculate_final_risk_scores(gene_group_extended):
    """Compute risk scores using Open Targets API data."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(query_open_targets, gene_group_extended['Gene_ID']))

    filtered_results = [df for df in results if df is not None]
    if not filtered_results:
        return None

    all_disease_associations = pd.concat(filtered_results, ignore_index=True)

    merged_df = all_disease_associations.merge(gene_group_extended, on="Gene_ID", how="left")

    merged_df["score"] = pd.to_numeric(merged_df["score"], errors='coerce').fillna(0)
    merged_df["Calculated_Risk_Score"] = (merged_df["Weighted_Gene_Score"] / merged_df["Gene_Length_bp"]) + merged_df["score"] * 10

    return merged_df

def filter_by_disease(merged_df, disease_name=None):
    """Filter genes associated with multiple diseases, remove duplicates, and save results."""
    if disease_name:
        if isinstance(disease_name, list):
            pattern = '|'.join(disease_name)
        else:
            pattern = disease_name
        filtered_df = merged_df[merged_df["disease_name"].astype(str).str.contains(pattern, case=False, na=False)]
    else:
        filtered_df = merged_df  

    expected_order = ["Gene", "Gene_ID", "disease_name", "disease_id", "Calculated_Risk_Score", "score"]
    filtered_df = filtered_df[expected_order].drop_duplicates()

    filtered_df.to_csv("results/filtered_by_disease.csv", index=False)
    print(f"Results saved as 'results/filtered_by_disease.csv'")

    return filtered_df

# ---- MAIN FUNCTION ---- #
def main(df, disease_name=None):
    """Main function to process input DataFrame and return filtered genes based on disease."""
    print("Starting Optimized DeepCV Risk Assessment...")

    df_unique = preprocess_data(df)
    print(f"Preprocessed {len(df_unique)} unique variants.")

    gene_group_sorted = calculate_gene_scores(df_unique)
    print(f"Calculated gene scores for {len(gene_group_sorted)} genes.")

    gene_group_extended = integrate_gene_information(df_unique, gene_group_sorted)
    print("Integrated gene length information.")

    merged_df = calculate_final_risk_scores(gene_group_extended)

    if merged_df is not None:
        print(f"Computed risk scores for {len(merged_df)} disease associations.")
        return filter_by_disease(merged_df, disease_name)
    else:
        print("No disease associations found.")
        return None

# Example usage:
# df = pd.read_csv("your_filtered_vcf_data.csv")
# final_results = main(df, disease_name=["type 2 diabetes mellitus", "hypertension"])
# print(final_results.head())
