import pandas as pd
import requests
import concurrent.futures  #Faster parallel API requests
import time  # MINIMAL CHANGE: Added for retry delay
import os

#Global cache for API responses (to avoid duplicate calls)
GENE_LENGTH_CACHE = {}
OPEN_TARGETS_CACHE = {}

def preprocess_data(df):
    """Preprocess input DataFrame: convert scores to numeric and compute initial risk scores."""
    # Create placeholder columns if they don't exist
    if 'clin_sig_score' not in df.columns:
        df['clin_sig_score'] = 0.1  # Default low clinical significance
    else:
        df['clin_sig_score'] = pd.to_numeric(df['clin_sig_score'], errors='coerce').fillna(0)
    
    if 'phenotype_or_disease' not in df.columns:
        df['phenotype_or_disease'] = 0.1  # Default low phenotype score
    else:
        df['phenotype_or_disease'] = pd.to_numeric(df['phenotype_or_disease'], errors='coerce').fillna(0)

    # Map Impact to numeric scores
    impact_scores = {'LOW': 1, 'MODIFIER': 2, 'MODERATE': 3, 'HIGH': 4}
    df['Impact_Score'] = df['Impact'].map(impact_scores).fillna(0)

    # Calculate risk score based on available data
    # If we don't have clinical significance scores, use Impact as primary risk factor
    df['Risk_Score'] = df['clin_sig_score'] + df['phenotype_or_disease'] + (df['Impact_Score'] * 0.1)

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
    # FIXED: Check for invalid gene_id first
    if not gene_id or pd.isna(gene_id) or gene_id == '-' or gene_id == '':
        return None
        
    if gene_id in GENE_LENGTH_CACHE:
        return GENE_LENGTH_CACHE[gene_id]

    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    
    # MINIMAL CHANGE: Simple retry with longer timeout
    for attempt in range(2):  # Try twice
        try:
            response = requests.get(url, headers=headers, timeout=15)  # CHANGED: 5 -> 15 seconds
            if response.ok:
                gene_data = response.json()
                # FIXED: Check if gene_data is None
                if gene_data is None:
                    GENE_LENGTH_CACHE[gene_id] = None
                    return None
                    
                # FIXED: Safely access nested dictionary
                transcripts = gene_data.get("Transcript", []) if gene_data else []
                gene_length = 0
                
                for transcript in transcripts:
                    if transcript is not None:
                        exons = transcript.get("Exon", [])
                        for exon in exons:
                            if exon is not None and "end" in exon and "start" in exon:
                                gene_length += exon["end"] - exon["start"] + 1
                
                GENE_LENGTH_CACHE[gene_id] = gene_length if gene_length > 0 else None
                return GENE_LENGTH_CACHE[gene_id]
        except requests.RequestException:
            if attempt == 0:  # MINIMAL CHANGE: Retry once after delay
                time.sleep(1)
                continue
            GENE_LENGTH_CACHE[gene_id] = None
            return None
    
    GENE_LENGTH_CACHE[gene_id] = None
    return None

def query_open_targets(gene_id):
    """Fetch disease associations from Open Targets API (optimized with caching)."""
    # FIXED: Check for invalid gene_id first
    if not gene_id or pd.isna(gene_id) or gene_id == '-' or gene_id == '':
        return None
        
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
    
    # MINIMAL CHANGE: Simple retry with longer timeout
    for attempt in range(2):  # Try twice
        try:
            response = requests.post(url, json={'query': query, 'variables': {"geneId": gene_id}}, timeout=15)  # CHANGED: 5 -> 15 seconds
            if response.ok:
                data = response.json()
                
                # FIXED: Check if data is None before accessing
                if data is None:
                    OPEN_TARGETS_CACHE[gene_id] = None
                    return None
                
                # FIXED: Safely navigate nested dictionary structure
                data_section = data.get('data') if data else None
                if data_section is None:
                    OPEN_TARGETS_CACHE[gene_id] = None
                    return None
                
                target_section = data_section.get('target') if data_section else None
                if target_section is None:
                    OPEN_TARGETS_CACHE[gene_id] = None
                    return None
                
                associated_diseases = target_section.get('associatedDiseases') if target_section else None
                if associated_diseases is None:
                    OPEN_TARGETS_CACHE[gene_id] = None
                    return None
                
                rows = associated_diseases.get('rows', []) if associated_diseases else []
                
                if rows:
                    # FIXED: Safely process rows
                    processed_rows = []
                    for row in rows:
                        if row is not None and 'disease' in row and 'score' in row:
                            disease = row['disease']
                            if disease is not None and isinstance(disease, dict):
                                processed_row = {
                                    'disease_id': disease.get('id'),
                                    'disease_name': disease.get('name'),
                                    'score': row['score'],
                                    'Gene_ID': gene_id
                                }
                                processed_rows.append(processed_row)
                    
                    if processed_rows:
                        df = pd.DataFrame(processed_rows)
                        OPEN_TARGETS_CACHE[gene_id] = df
                        return df
                
                OPEN_TARGETS_CACHE[gene_id] = None
                return None
                
        except requests.RequestException:
            if attempt == 0:  # MINIMAL CHANGE: Retry once after delay
                time.sleep(1)
                continue
            OPEN_TARGETS_CACHE[gene_id] = None
            return None
    
    OPEN_TARGETS_CACHE[gene_id] = None
    return None

def integrate_gene_information(df_unique, gene_group_sorted):
    """Fetch gene lengths in parallel and merge with gene scores."""
    gene_info = df_unique[['Gene', 'Gene_ID', 'Feature_Type']].drop_duplicates()

    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        gene_lengths = list(executor.map(query_gene_length, gene_info["Gene_ID"]))

    gene_info["Gene_Length_bp"] = gene_lengths
    
    # FIXED: Replace None values with default gene length
    gene_info["Gene_Length_bp"] = gene_info["Gene_Length_bp"].fillna(1000)

    return gene_group_sorted.merge(gene_info, on=['Gene', 'Gene_ID'], how='left')

def calculate_final_risk_scores(gene_group_extended):
    """Compute risk scores using Open Targets API data."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(query_open_targets, gene_group_extended['Gene_ID']))

    # FIXED: Filter out None results more safely
    filtered_results = [df for df in results if df is not None and not df.empty]
    if not filtered_results:
        return None

    all_disease_associations = pd.concat(filtered_results, ignore_index=True)

    merged_df = all_disease_associations.merge(gene_group_extended, on="Gene_ID", how="left")

    merged_df["score"] = pd.to_numeric(merged_df["score"], errors='coerce').fillna(0)
    
    # FIXED: Prevent division by zero
    merged_df["Gene_Length_bp"] = merged_df["Gene_Length_bp"].fillna(1000)
    merged_df = merged_df[merged_df["Gene_Length_bp"] > 0]
    
    merged_df["Calculated_Risk_Score"] = (merged_df["Weighted_Gene_Score"] / merged_df["Gene_Length_bp"]) + merged_df["score"] * 10

    return merged_df

def create_cleaned_risk_assessment(gene_group_extended, disease_associations_df=None):
    """
    Create a CLEANED risk assessment that ONLY includes genes with meaningful disease associations.
    
    CLEANED VERSION: Filters out:
    1. Genes with "No known associations"
    2. Genes with missing names (showing as "-")
    3. Genes with empty or invalid Gene_IDs
    """
    if disease_associations_df is not None and not disease_associations_df.empty:
        # Start with disease associations (already meaningful)
        cleaned_df = disease_associations_df.copy()
        
        # CLEANING STEP 1: Remove genes with missing names
        cleaned_df = cleaned_df[
            (cleaned_df['Gene'] != '-') & 
            (cleaned_df['Gene'].notna()) & 
            (cleaned_df['Gene'] != '') &
            (cleaned_df['Gene'] != 'nan')
        ]
        
        # CLEANING STEP 2: Remove genes with invalid Gene_IDs
        cleaned_df = cleaned_df[
            (cleaned_df['Gene_ID'] != '-') & 
            (cleaned_df['Gene_ID'].notna()) & 
            (cleaned_df['Gene_ID'] != '') &
            (cleaned_df['Gene_ID'] != 'nan')
        ]
        
        # CLEANING STEP 3: Remove genes with very low disease association scores (< 0.1)
        cleaned_df = cleaned_df[cleaned_df['score'] >= 0.1]
        
        # Reorder columns for better readability
        column_order = ['Gene', 'Gene_ID', 'disease_name', 'disease_id', 'score', 'Calculated_Risk_Score',
                       'Num_SNPs', 'Total_Risk_Score', 'Weighted_Gene_Score', 'Feature_Type', 'Gene_Length_bp']
        cleaned_df = cleaned_df[column_order]
        
        # Sort by disease association score (highest first)
        cleaned_df = cleaned_df.sort_values('score', ascending=False).reset_index(drop=True)
        
        print(f"Cleaned risk assessment includes {len(cleaned_df)} meaningful disease associations.")
        print(f"Filtered out genes with missing names, invalid IDs, and low association scores.")
        
        return cleaned_df
        
    else:
        # No disease associations found - return empty DataFrame
        print("No disease associations found. Cleaned risk assessment will be empty.")
        return pd.DataFrame()

def create_comprehensive_summary(gene_group_extended, cleaned_df):
    """
    Create a summary of all genes processed vs. meaningful associations found.
    """
    total_genes = len(gene_group_extended)
    meaningful_genes = len(cleaned_df)
    
    # Count genes with missing names
    missing_names = len(gene_group_extended[
        (gene_group_extended['Gene'] == '-') | 
        (gene_group_extended['Gene'].isna()) | 
        (gene_group_extended['Gene'] == '')
    ])
    
    # Count genes with invalid Gene_IDs
    invalid_ids = len(gene_group_extended[
        (gene_group_extended['Gene_ID'] == '-') | 
        (gene_group_extended['Gene_ID'].isna()) | 
        (gene_group_extended['Gene_ID'] == '')
    ])
    
    summary = {
        'total_genes_processed': total_genes,
        'genes_with_missing_names': missing_names,
        'genes_with_invalid_ids': invalid_ids,
        'meaningful_disease_associations': meaningful_genes,
        'filtering_efficiency': f"{(meaningful_genes/total_genes)*100:.1f}%" if total_genes > 0 else "0%"
    }
    
    print("\\n" + "="*50)
    print("CLEANING SUMMARY:")
    print("="*50)
    print(f"Total genes processed: {summary['total_genes_processed']}")
    print(f"Genes with missing names: {summary['genes_with_missing_names']}")
    print(f"Genes with invalid IDs: {summary['genes_with_invalid_ids']}")
    print(f"Meaningful disease associations: {summary['meaningful_disease_associations']}")
    print(f"Filtering efficiency: {summary['filtering_efficiency']}")
    print("="*50)
    
    return summary

def filter_by_disease(merged_df, disease_name=None, output_basename=""):
    """
    Filter genes associated with specific diseases and save results ONLY when requested.
    
    FIXED: Only generate filtered_by_disease.csv when:
    1. disease_name is specified (user requested disease filtering)
    2. Disease associations are found in the data
    """
    if disease_name is None:
        # No disease filtering requested - don't generate filtered_by_disease.csv
        print("No disease filtering requested. Skipping filtered_by_disease.csv generation.")
        return None
    
    if merged_df is None or merged_df.empty:
        # No disease associations found - don't generate filtered_by_disease.csv
        print("No disease associations found. Skipping filtered_by_disease.csv generation.")
        return None
    
    # Apply disease filtering
    if isinstance(disease_name, list):
        pattern = '|'.join(disease_name)
    else:
        pattern = disease_name
    
    filtered_df = merged_df[merged_df["disease_name"].astype(str).str.contains(pattern, case=False, na=False)]
    
    if filtered_df.empty:
        # No matches found for the specified disease
        print(f"No associations found for disease: {disease_name}")
        print("Skipping filtered_by_disease.csv generation.")
        return None
    
    # Generate the filtered file only when we have matches
    expected_order = ["Gene", "Gene_ID", "disease_name", "disease_id", "Calculated_Risk_Score", "score"]
    filtered_df = filtered_df[expected_order].drop_duplicates()
    
    # Use output basename for unique filenames
    output_file = f"results/filtered_by_disease{output_basename}.csv"
    filtered_df.to_csv(output_file, index=False)
    print(f"Disease-filtered results saved as '{output_file}'")
    print(f"Found {len(filtered_df)} associations for disease: {disease_name}")

    return filtered_df

# ---- MAIN FUNCTION ---- #
def main(input_file, output_dir, vcf_basename, disease_name=None):
    """
    Main function to process input CSV file and return CLEANED genes with meaningful disease associations.
    
    Args:
        input_file (str): Path to input CSV file from Step 2
        output_dir (str): Output directory
        vcf_basename (str): Base name of VCF file
        disease_name (str, optional): Disease name for filtering
        
    Returns:
        pandas.DataFrame: Cleaned risk assessment data
    """
    import pandas as pd
    import os
    
    # Load the input CSV file
    print(f"Loading input file: {input_file}")
    df = pd.read_csv(input_file)
    
    # Create output basename for this analysis
    output_basename = f"_{vcf_basename}" if vcf_basename else ""
    
    # Call the original main function logic
    cleaned_risk_df = main_analysis(df, disease_name, output_basename, output_dir)
    
    return cleaned_risk_df

def main_analysis(df, disease_name=None, output_basename="", output_dir="results"):
    """
    Main function to process input DataFrame and return CLEANED genes with meaningful disease associations.
    
    CLEANED VERSION: Only saves genes with actual disease associations and valid names/IDs.
    """
    print("Starting Optimized DeepCV Risk Assessment...")

    df_unique = preprocess_data(df)
    print(f"Preprocessed {len(df_unique)} unique variants.")

    gene_group_sorted = calculate_gene_scores(df_unique)
    print(f"Calculated gene scores for {len(gene_group_sorted)} genes.")

    gene_group_extended = integrate_gene_information(df_unique, gene_group_sorted)
    print("Integrated gene length information.")

    # ALWAYS try to get disease associations
    print("Fetching disease associations for risk assessment...")
    disease_associations_df = calculate_final_risk_scores(gene_group_extended)
    
    if disease_associations_df is not None:
        print(f"Found {len(disease_associations_df)} disease associations.")
    else:
        print("No disease associations found in Open Targets database.")

    # CLEANED VERSION: Create cleaned risk assessment (only meaningful associations)
    cleaned_risk_df = create_cleaned_risk_assessment(gene_group_extended, disease_associations_df)
    
    # Create summary of cleaning process
    summary = create_comprehensive_summary(gene_group_extended, cleaned_risk_df)
    
    # Save the cleaned risk assessment file (ONLY meaningful associations)
    main_output_file = os.path.join(output_dir, f"final_risk_assessment{output_basename}.csv")
    if not cleaned_risk_df.empty:
        cleaned_risk_df.to_csv(main_output_file, index=False)
        print(f"Cleaned risk assessment saved as '{main_output_file}'")
        print(f"File contains {len(cleaned_risk_df)} meaningful disease associations.")
    else:
        # Create empty file with headers
        empty_df = pd.DataFrame(columns=['Gene', 'Gene_ID', 'disease_name', 'disease_id', 'score', 'Calculated_Risk_Score',
                                       'Num_SNPs', 'Total_Risk_Score', 'Weighted_Gene_Score', 'Feature_Type', 'Gene_Length_bp'])
        empty_df.to_csv(main_output_file, index=False)
        print(f"No meaningful associations found. Empty file saved as '{main_output_file}'")

    # Only generate filtered_by_disease.csv if disease_name is specified
    if disease_name is not None:
        print(f"Disease filtering requested for: {disease_name}")
        filter_by_disease(disease_associations_df, disease_name, output_basename)
    else:
        print("No disease filtering requested. Analysis complete.")
    
    # FIXED: Always return the cleaned risk assessment DataFrame
    return cleaned_risk_df

# Example usage:
# df = pd.read_csv("your_filtered_vcf_data.csv")
# final_results = main(df)


