import pandas as pd
from pathlib import Path
from cyvcf2 import VCF

def parse_vcf(vcf_file, output_csv=None):
    """
    Parse a VCF file to extract detailed variant data, filter variants, and optionally save the results.

    Parameters:
        vcf_file (str): Path to the input VCF file.
        output_csv (str): Path to save filtered CSV. If None, doesn't save to file.

    Returns:
        pd.DataFrame: Filtered DataFrame of variant data.
    """
    vcf_reader = VCF(vcf_file)
    data = []

    for variant in vcf_reader:
        # Extract fields safely, handling missing values
        variant_data = {
            "Chromosome": variant.CHROM,
            "Position": variant.POS,
            "ID": variant.ID if variant.ID else ".",
            "Reference": variant.REF,
            "Alternative": ','.join(variant.ALT) if variant.ALT else ".",
            "Qual": variant.QUAL if variant.QUAL is not None else 0,
            "Filter": variant.FILTER if variant.FILTER else "PASS",
            "Annotation": variant.INFO.get("Annotation", None),
            "Impact": variant.INFO.get("Impact", None),
            "Gene": variant.INFO.get("Gene", None),
            "Gene_ID": variant.INFO.get("Gene_ID", None),
            "Feature_Type": variant.INFO.get("Feature_Type", None),
            "Feature_ID": variant.INFO.get("Feature_ID", None),
            "VAF": variant.INFO.get("VAF", None),
        }

        # Extract genotype details
        try:
            gt = '/'.join(map(str, variant.genotypes[0][:2])) if variant.genotypes else "."
            gq = variant.format('GQ')[0][0] if variant.format('GQ') is not None else 0
            dp = variant.format('DP')[0][0] if variant.format('DP') is not None else 0
        except IndexError:
            gt, gq, dp = ".", 0, 0  # Handle cases where genotype fields are missing

        variant_data.update({
            "GT": gt,
            "GQ": gq,
            "DP": dp,
        })

        data.append(variant_data)

    # Convert to DataFrame
    df = pd.DataFrame(data)

    if df.empty:
        print("No variants found in VCF file.")
        return pd.DataFrame()  # Return empty DataFrame if no data

    # Drop variants with missing mandatory fields
    df.dropna(subset=["Chromosome", "Position", "Reference", "Alternative"], inplace=True)

    print(f"Loaded {len(df)} variants from {vcf_file}")

    # Apply filtering criteria
    filtered_df = df[
        (df['Qual'] > 45) & 
        (df['Filter'] == 'PASS') & 
        (df['DP'] > 40) & 
        (df['GQ'] > 60)
    ]

    print(f"Filtered {len(filtered_df)} variants (passed quality thresholds)")

    # Drop unnecessary columns
    keep_columns = ["Chromosome", "Position", "ID", "Reference", "Alternative", 
                    "Annotation", "Impact", "Gene", "Gene_ID", "Feature_Type", "Feature_ID", "VAF"]
    
    filtered_df = filtered_df[keep_columns]

    # Remove duplicates based on key fields
    filtered_df = filtered_df.drop_duplicates(subset=["Chromosome", "Position", "Reference", "Alternative"]).reset_index(drop=True)

    print(f" {len(filtered_df)} unique variants remain after removing duplicates")

    # Save if output path is provided
    if output_csv:
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        filtered_df.to_csv(output_csv, index=False)
        print(f" Filtered data saved to {output_csv}")

    return filtered_df

def main(vcf_file, output_dir):
    """
    Main function for Step 1: VCF Processing
    
    Args:
        vcf_file (str): Path to input VCF file
        output_dir (str): Output directory
        
    Returns:
        pd.DataFrame: Processed variants DataFrame
    """
    import os
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate output filename
    vcf_basename = Path(vcf_file).stem
    output_csv = os.path.join(output_dir, f"processed_variants_{vcf_basename}.csv")
    
    # Parse VCF file
    try:
        filtered_df = parse_vcf(vcf_file, output_csv)
        return filtered_df
    except Exception as e:
        print(f"Error processing VCF file: {e}")
        return None

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python step1_processing.py <vcf_file> <output_dir>")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    result = main(vcf_file, output_dir)
    if result is not None:
        print(f"Step 1 completed successfully. Processed {len(result)} variants.")
    else:
        print("Step 1 failed.")
        sys.exit(1)

