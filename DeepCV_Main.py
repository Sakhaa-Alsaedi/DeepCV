"""
DeepCV Workflow - Variant Processing, Annotation & Risk Assessment (Demo Version)

===================================================
This script processes VCF files by parsing, filtering, 
annotating SNPs using DeepCV modules, and assessing genetic risk.

Demo Version: ALPHA

Creative & Informative Environment:
   - Designed to streamline genomic data analysis.
   - User-friendly logging and progress tracking.
   - Provides clear output with step-by-step updates.

Requirements:
   - Ensure all dependencies are installed:
     ```bash
     pip install -r requirements.txt
     ```

Contact:
   If you have further questions, please contact:
   sakhaa.alsaedi@kaust.edu.sa
===================================================
"""

import os
import logging
import sys
import time
import subprocess
import pandas as pd
import step1_processing      # for VCF processing
import step2_SNP_annotation  # for SNP annotations
import step3_risk_assessment # for risk scoring

# Define version
DEMO_VERSION = "ALPHA"

# Ensure dependencies are installed
try:
    subprocess.run(["pip", "install", "-r", "requirements.txt"], check=True)
except subprocess.CalledProcessError:
    print("Failed to install dependencies. Please run: pip install -r requirements.txt")

# Create results directory
RESULTS_DIR = "results"
os.makedirs(RESULTS_DIR, exist_ok=True)

def format_time(seconds):
    """Formats time in minutes and seconds for logging."""
    minutes = int(seconds // 60)
    seconds = int(seconds % 60)
    return f"{minutes}m {seconds}s" if minutes > 0 else f"{seconds}s"

def main(input_vcf=None, output_csv="final_annotated_output.csv", disease_name=None, run_step1=True, run_step2=True, run_step3=True):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    print("\n ******************************************** ")
    print(f"*         DeepCV SNP Annotation             *")
    print(f"*       Demo Version: {DEMO_VERSION}                 *")
    print("*                                           *")
    print("*  Variant Processing & Risk Assessment     *")
    print(" ******************************************** \n")


    start_time = time.time()

    try:
        ## --- Step 1: Parse & Filter VCF File ---
        processed_output_path = os.path.join(RESULTS_DIR, "processed_variants.csv")

        if run_step1:
            if input_vcf is None:
                logging.error(f"[{DEMO_VERSION}] No input VCF file specified. Please provide a valid file.")
                sys.exit(1)

            logging.info(f"[{DEMO_VERSION}] Starting VCF parsing and filtering...")
            step1_start = time.time()
            processed_df = step1_processing.parse_vcf(input_vcf)
            step1_time = time.time() - step1_start
            logging.info(f"[{DEMO_VERSION}] VCF parsing complete. {len(processed_df)} variants after filtering. Time taken: {format_time(step1_time)}")

            # Save processed output
            processed_df.to_csv(processed_output_path, index=False)
            logging.info(f"[{DEMO_VERSION}] Processed VCF saved to: {processed_output_path}")

        elif os.path.exists(processed_output_path):
            logging.info(f"[{DEMO_VERSION}] Skipping Step 1 - Loading existing processed variants")
            processed_df = pd.read_csv(processed_output_path)
        else:
            logging.error(f"[{DEMO_VERSION}] Missing required file: {processed_output_path}. Run Step 1 first.")
            sys.exit(1)

        ## --- Step 2: Annotate SNPs ---
        annotated_output_path = os.path.join(RESULTS_DIR, "final_annotated_output.csv")

        if run_step2:
            logging.info(f"[{DEMO_VERSION}] Starting SNP annotation...")
            step2_start = time.time()
            annotated_df = step2_SNP_annotation.annotate_snp_level(processed_df)
            step2_time = time.time() - step2_start
            logging.info(f"[{DEMO_VERSION}] SNP annotation complete. Annotated {len(annotated_df)} variants. Time taken: {format_time(step2_time)}")

            # Save annotated output
            annotated_df.to_csv(annotated_output_path, index=False)
            logging.info(f"[{DEMO_VERSION}] SNP annotation saved to: {annotated_output_path}")

        elif os.path.exists(annotated_output_path):
            logging.info(f"[{DEMO_VERSION}] Skipping Step 2 - Loading existing annotated SNPs")
            annotated_df = pd.read_csv(annotated_output_path)
        else:
            logging.error(f"[{DEMO_VERSION}] Missing required file: {annotated_output_path}. Run Step 2 first.")
            sys.exit(1)

        ## --- Step 3: Perform Risk Assessment ---
        final_risk_path = os.path.join(RESULTS_DIR, "final_risk_assessment.csv")

        if run_step3:
            logging.info(f"[{DEMO_VERSION}] Starting genetic risk assessment...")

            step3_start = time.time()
            risk_assessment_df = step3_risk_assessment.main(annotated_df, disease_name)
            step3_time = time.time() - step3_start

            if risk_assessment_df is not None:
                risk_assessment_df.to_csv(final_risk_path, index=False)
                logging.info(f"[{DEMO_VERSION}] Risk assessment results saved to: {final_risk_path}")
                logging.info(f"[{DEMO_VERSION}] Processed {len(risk_assessment_df)} gene-disease associations. Time taken: {format_time(step3_time)}")
            else:
                logging.warning(f"[{DEMO_VERSION}] Risk assessment did not return any results.")

        else:
            logging.info(f"[{DEMO_VERSION}] Skipping Step 3 - Risk assessment was not executed.")

        ## --- Total Time ---
        total_time = time.time() - start_time
        logging.info(f"[{DEMO_VERSION}] Total processing time: {format_time(total_time)}")

    except Exception as e:
        logging.error(f"[{DEMO_VERSION}] An error occurred: {e}")
        logging.error("If the issue persists, please contact: sakhaa.alsaedi@kaust.edu.sa")
        sys.exit(1)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Process VCF file, annotate SNPs, and assess genetic risk.")
    parser.add_argument("--input_vcf", type=str, help="Path to input VCF file")
    parser.add_argument("--output_csv", type=str, default="results/final_annotated_output.csv", help="Path to final SNP annotation output file")
    parser.add_argument("--disease_name", type=str, nargs="+", default=None, help="Filter results for one or more diseases (optional)")
    
    # Add flags for skipping specific steps
    parser.add_argument("--skip_step1", action="store_true", help="Skip VCF parsing and filtering (Step 1)")
    parser.add_argument("--skip_step2", action="store_true", help="Skip SNP annotation (Step 2)")
    parser.add_argument("--skip_step3", action="store_true", help="Skip genetic risk assessment (Step 3)")

    args = parser.parse_args()

    main(
        input_vcf=args.input_vcf, 
        output_csv=args.output_csv, 
        disease_name=args.disease_name, 
        run_step1=not args.skip_step1,
        run_step2=not args.skip_step2,
        run_step3=not args.skip_step3
    )
