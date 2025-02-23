# DeepCV (Alpha version, 10 Feb 2025) - Genetic Risk Scoring Workflow  

![Image](https://github.com/user-attachments/assets/6733cd05-ea8f-4194-a533-12be4e781858)


## 🔹 Overview
## Deep Common Variant Analysis Tool for Risk Assessment. Decode your genome profile, assess disease risk, protect and proactive health with DeepCV!

The **DeepCV Alpha Version** is a computational workflow for **variant annotation and genetic risk scoring using statistical methods and LLM-Powered multiomics risk scoring system for each gene**
*This pipeline has been tested using **synthetic data generated orgnal form from `NA12878.ann.vcf.gz`** and **GWAS catalog database**.


### ✅ Purpose  
- Computes **genetic risk scores** for **each gene associated with a disease**.  
- Provides **a filterable gene list** for disease-specific risk analysis.

## 📂 DeepCV Workflow - Files & Descriptions

### **📌 How the Data is Used in DeepCV**
1️⃣ **The user provides a VCF file** (`.vcf` or `.vcf.gz`) as input.  
2️⃣ **DeepCV extracts relevant SNPs and mutations**, removing non-informative variants.  
3️⃣ **Annotation is performed** using clinical databases and functional impact predictors.  
4️⃣ **Risk scores are computed**, generating a list of **genes associated with disease risk**.  
5️⃣ **Results are stored in structured output files**, such as `final_risk_assessment.csv`.  

Below is an overview of the **DeepCV workflow** and the corresponding scripts.

| **Filename**                   | **Description**                                                 |
|--------------------------------|-----------------------------------------------------------------|
| `DeepCV_Main.py`              | The main script that runs the entire pipeline.                 |
| `environment.yml`             | Conda environment configuration file for dependency management. |
| `setup_deepcv.py`             | Automates setup, environment creation, and installation.       |
| `step1_processing.py`         | Parses and filters VCF files for variant processing.           |
| `step2_SNP_annotation.py`     | Annotates SNPs using functional impact scores and clinical data. |
| `step3_risk_assessment.py`    | Computes genetic risk scores based on annotated SNPs.         |

---

# Data Overview

### Input  
- **Annotated VCF file** (processed using **Variant Effect Predictor (VEP)**).

DeepCV processes the following **synthetic VCF files** to simulate real-world variant annotation and risk assessment.

| **Dataset Name**        | **File Format** | **Description**                                  |
|------------------------|--------------|------------------------------------------------|
| `Demo.vcf`    | `.vcf`       | Synthetic VCF file for hypertension (HTN) and T2D genetic risk assessment. |
| `Demo.vcf.gz` | `.vcf.gz`    | Compressed version of `Demo.vcf` for optimized storage. |

### Output Files & Directories  

- **List of genes with genetic risk scores** and their association with diseases saved as CSV file.
- **Filtered gene lists with risk scores** based on disease-specific queries.

DeepCV saves all results inside the results/ directory:

After execution, DeepCV generates several files:

| **Filename**                   | **Description**                                             |
|--------------------------------|-------------------------------------------------------------|
| Step1 output: `processed_variants.csv` | filtered and processed vcf2csv version of `Demovcf`  |
| Step2 output: `final_annotated_output.csv`   | Annotated SNPs with functional & clinical significance     |
| Step3 output: `final_risk_assessment.csv`    | Computed genetic risk scores for all genes                 |
| Step3 output: `filtered_by_disease.csv`      | Disease-specific filtered results                          |


# DeepCV Alpha - Genetic Risk Scoring Workflow  
The workflow consists of four main steps: VCF processing and filtering, annotation, risk score assessment, and risk graph construction.
The repository provides code for steps 1 to 3.Below is a detailed description of each step and how it will be executed using DeepCV:

# **DeepCV: Command Guide & Usage Instructions**

- To run DeepCV Alpha, use the following command:

## 1. Setting Up DeepCV (First-Time Users)
Ensure dependencies are installed before running DeepCV:

- create a Conda environment:

```bash
conda env create -f environment.yml
conda activate DeepCV
```

## 2. Executing DeepCV Alpha Workflow Steps

🔹 Full Pipeline with All Diseases
To run **all steps** (VCF processing, SNP annotation, and risk assessment) in one command:

```bash
python DeepCV_Main.py --input_vcf ./path/input.vcf 
```

## 🔹 Filtering Genes by Disease  

To filter the **genetic risk scores** based on a specific disease, use:
```bash
python DeepCV_Main.py input.vcf --disease_name "Your Disease Name"

```
for example: 

```bash
python DeepCV_Main.py --input_vcf ./path/input.vcf --disease_name "Type 2 Diabetes Mellitus"
```

---

📩 **For more details, read our upcoming research paper!**  
📩 **Contact:** [sakhaa.alsaedi@kaust.edu.sa](mailto:sakhaa.alsaedi@kaust.edu.sa)  



