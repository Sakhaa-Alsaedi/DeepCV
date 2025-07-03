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

def load_cache(cache_file):
    """Load cache from file if it exists"""
    import pickle
    import os
    try:
        if os.path.exists(cache_file):
            with open(cache_file, 'rb') as f:
                return pickle.load(f)
    except:
        pass
    return {}

def save_cache(cache_dict, cache_file):
    """Save cache to file"""
    import pickle
    import os
    try:
        os.makedirs(os.path.dirname(cache_file), exist_ok=True)
        with open(cache_file, 'wb') as f:
            pickle.dump(cache_dict, f)
    except:
        pass

# Load caches at startup
VEP_CACHE = load_cache('cache/vep_cache.pkl')
GWAS_CACHE = load_cache('cache/gwas_cache.pkl')
DBSNP_CACHE = load_cache('cache/dbsnp_cache.pkl')

print(f"Loaded caches: VEP={len(VEP_CACHE)}, GWAS={len(GWAS_CACHE)}, dbSNP={len(DBSNP_CACHE)}")

def create_variant_key(chrom, pos, ref, alt):
    """Create a unique key for variant caching"""
    return f"{chrom}:{pos}:{ref}:{alt}"

def query_vep_batch(hgvs_list, batch_size=200):
    """
    Query VEP API in batches for better performance
    """
    results = {}
    
    # Split into batches
    for i in range(0, len(hgvs_list), batch_size):
        batch = hgvs_list[i:i+batch_size]
        
        # Prepare batch request
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        data = {"hgvs_notations": batch}
        
        try:
            response = requests.post(
                "https://rest.ensembl.org/vep/human/hgvs",
                headers=headers,
                json=data,
                timeout=30
            )
            
            if response.ok:
                batch_results = response.json()
                
                # Map results back to HGVS notations
                for j, result in enumerate(batch_results):
                    if j < len(batch):
                        hgvs = batch[j]
                        results[hgvs] = result
                        # Cache the result
                        VEP_CACHE[hgvs] = result
            else:
                print(f"VEP batch request failed: {response.status_code}")
                
        except requests.RequestException as e:
            print(f"VEP batch request error: {e}")
            
        # Rate limiting
        time.sleep(0.1)
    
    return results

def query_vep_individual(hgvs_notation):
    """
    Fetch annotation from Ensembl VEP (individual query as fallback).
    Caches results to avoid redundant queries.
    """
    if hgvs_notation in VEP_CACHE:
        return VEP_CACHE[hgvs_notation]  # Return cached result

    url = f"https://rest.ensembl.org/vep/human/hgvs/{hgvs_notation}?content-type=application/json"
    try:
        response = requests.get(url, timeout=5)
        if response.ok:
            result = response.json()
            VEP_CACHE[hgvs_notation] = result
            return result
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

    # Initialize default values
    parsed_info = {
        "consequence_terms": "NR",
        "impact": "NR", 
        "gene_symbol": "NR",
        "gene_id": "NR",
        "transcript_id": "NR",
        "protein_id": "NR",
        "sift_prediction": "NR",
        "polyphen_prediction": "NR",
        "clinical_significance": "NR"
    }

    # Extract information from the first transcript
    transcript = vep_data[0]
    
    # Basic consequence information
    if 'consequence_terms' in transcript:
        parsed_info["consequence_terms"] = ','.join(transcript['consequence_terms'])
    
    if 'impact' in transcript:
        parsed_info["impact"] = transcript['impact']
    
    # Gene information
    if 'gene_symbol' in transcript:
        parsed_info["gene_symbol"] = transcript['gene_symbol']
    
    if 'gene_id' in transcript:
        parsed_info["gene_id"] = transcript['gene_id']
    
    if 'transcript_id' in transcript:
        parsed_info["transcript_id"] = transcript['transcript_id']
    
    if 'protein_id' in transcript:
        parsed_info["protein_id"] = transcript['protein_id']
    
    # Prediction scores
    if 'sift_prediction' in transcript:
        parsed_info["sift_prediction"] = transcript['sift_prediction']
    
    if 'polyphen_prediction' in transcript:
        parsed_info["polyphen_prediction"] = transcript['polyphen_prediction']
    
    # Clinical significance from colocated variants
    if 'colocated_variants' in transcript:
        clin_sigs = set()
        for variant in transcript['colocated_variants']:
            if 'clin_sig' in variant:
                clin_sigs.update(variant['clin_sig'])
        
        if clin_sigs:
            top_sig, _ = prioritize_clin_sig(clin_sigs)
            parsed_info["clinical_significance"] = top_sig

    return parsed_info

def annotate_snp_level_optimized(df, output_csv=None, batch_size=200, max_workers=50, save_interval=10000):
    """
    Optimized SNP-level annotation with batch processing and parallel execution.
    """
    print("Starting OPTIMIZED SNP annotation...")
    
    # Expand alternative alleles
    print("Expanding alternative alleles...")
    expanded_rows = []
    for _, row in df.iterrows():
        alts = str(row['Alternative']).split(',')
        for alt in alts:
            new_row = row.copy()
            new_row['Alternative'] = alt.strip()
            expanded_rows.append(new_row)
    
    df_expanded = pd.DataFrame(expanded_rows)
    
    # Remove duplicates
    print(f"Removing duplicates: {len(df)} -> {len(df_expanded)} unique variants")
    df_expanded = df_expanded.drop_duplicates(subset=['Chromosome', 'Position', 'Reference', 'Alternative']).reset_index(drop=True)
    
    # Create HGVS notations
    df_expanded['HGVS'] = df_expanded.apply(
        lambda row: f"{row['Chromosome']}:g.{row['Position']}{row['Reference']}>{row['Alternative']}", 
        axis=1
    )
    
    # Check cache hit rate
    uncached_variants = [hgvs for hgvs in df_expanded['HGVS'] if hgvs not in VEP_CACHE]
    cached_count = len(df_expanded) - len(uncached_variants)
    print(f"Cache hit rate: {cached_count}/{len(df_expanded)} ({cached_count/len(df_expanded)*100:.1f}%)")
    
    # Query VEP API with batch requests
    print("Querying VEP API with batch requests...")
    if uncached_variants:
        print(f"Querying VEP for {len(uncached_variants)} uncached variants...")
        
        # Use batch processing
        from tqdm import tqdm
        batch_results = {}
        
        for i in tqdm(range(0, len(uncached_variants), batch_size), desc="VEP Batches"):
            batch = uncached_variants[i:i+batch_size]
            batch_result = query_vep_batch(batch, batch_size)
            batch_results.update(batch_result)
            
            # Save cache periodically
            if i % save_interval == 0:
                save_cache(VEP_CACHE, 'cache/vep_cache.pkl')
    
    # Process VEP results
    print("Processing VEP results...")
    annotation_results = []
    
    for _, row in df_expanded.iterrows():
        hgvs = row['HGVS']
        vep_data = VEP_CACHE.get(hgvs, {})
        
        # Parse VEP data
        parsed_info = parse_snp_info(vep_data)
        
        # Combine original row with annotations
        annotated_row = row.to_dict()
        annotated_row.update(parsed_info)
        annotation_results.append(annotated_row)
    
    # Create final DataFrame
    annotated_df = pd.DataFrame(annotation_results)
    
    # Save results
    if output_csv:
        annotated_df.to_csv(output_csv, index=False)
        print(f"Annotations saved to {output_csv}")
    
    # Save caches
    save_cache(VEP_CACHE, 'cache/vep_cache.pkl')
    save_cache(GWAS_CACHE, 'cache/gwas_cache.pkl')
    save_cache(DBSNP_CACHE, 'cache/dbsnp_cache.pkl')
    
    print(f"OPTIMIZED annotation completed in {time.time() - time.time():.1f} seconds")
    print(f"Cache sizes: VEP={len(VEP_CACHE)}, GWAS={len(GWAS_CACHE)}, dbSNP={len(DBSNP_CACHE)}")
    
    return annotated_df

def annotate_snp_level_progressive(df, output_csv=None, chunk_size=50000):
    """
    Progressive annotation with chunking for very large datasets.
    """
    print("Starting PROGRESSIVE SNP annotation...")
    
    total_rows = len(df)
    annotated_chunks = []
    
    for i in range(0, total_rows, chunk_size):
        chunk = df.iloc[i:i+chunk_size]
        print(f"Processing chunk {i//chunk_size + 1}/{(total_rows-1)//chunk_size + 1}")
        
        # Process chunk with optimized method
        chunk_annotated = annotate_snp_level_optimized(chunk, output_csv=None)
        annotated_chunks.append(chunk_annotated)
        
        # Save intermediate results
        if output_csv:
            intermediate_file = output_csv.replace('.csv', f'_chunk_{i//chunk_size + 1}.csv')
            chunk_annotated.to_csv(intermediate_file, index=False)
    
    # Combine all chunks
    final_df = pd.concat(annotated_chunks, ignore_index=True)
    
    if output_csv:
        final_df.to_csv(output_csv, index=False)
        print(f"Final annotations saved to {output_csv}")
    
    return final_df

def annotate_snp_level(df, output_csv=None, method="optimized", **kwargs):
    """
    Main annotation function with multiple methods.
    
    Args:
        df: Input DataFrame
        output_csv: Output file path
        method: "optimized" or "progressive"
        **kwargs: Additional arguments for specific methods
    """
    if method == "optimized":
        return annotate_snp_level_optimized(df, output_csv, **kwargs)
    elif method == "progressive":
        return annotate_snp_level_progressive(df, output_csv, **kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")

def clear_cache():
    """Clear all caches"""
    global VEP_CACHE, GWAS_CACHE, DBSNP_CACHE
    VEP_CACHE.clear()
    GWAS_CACHE.clear()
    DBSNP_CACHE.clear()
    print("All caches cleared")

def cache_stats():
    """Print cache statistics"""
    print(f"Cache sizes: VEP={len(VEP_CACHE)}, GWAS={len(GWAS_CACHE)}, dbSNP={len(DBSNP_CACHE)}")

def main(input_csv, output_dir, vcf_basename):
    """
    Main function for Step 2: SNP Annotation
    
    Args:
        input_csv (str): Path to input CSV file from Step 1
        output_dir (str): Output directory
        vcf_basename (str): Base name of VCF file
        
    Returns:
        pd.DataFrame: Annotated variants DataFrame
    """
    import os
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load input data
    try:
        df = pd.read_csv(input_csv)
        print(f"Loaded {len(df)} variants for annotation")
    except Exception as e:
        print(f"Error loading input file: {e}")
        return None
    
    # Generate output filename
    output_csv = os.path.join(output_dir, f"final_annotated_output_{vcf_basename}.csv")
    
    # Annotate variants
    try:
        annotated_df = annotate_snp_level(df, output_csv, method="optimized")
        return annotated_df
    except Exception as e:
        print(f"Error annotating variants: {e}")
        return None

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python step2_SNP_annotation.py <input_csv> <output_dir> <vcf_basename>")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_dir = sys.argv[2]
    vcf_basename = sys.argv[3]
    
    result = main(input_csv, output_dir, vcf_basename)
    if result is not None:
        print(f"Step 2 completed successfully. Annotated {len(result)} variants.")
    else:
        print("Step 2 failed.")
        sys.exit(1)

