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
| `requirements.txt`            | List of required Python packages for installation via pip.     |
| `setup_deepcv.py`             | Automates setup, environment creation, and installation.       |
| `step1_processing.py`         | Parses and filters VCF files for variant processing.           |
| `step2_SNP_annotation.py`     | Annotates SNPs using functional impact scores and clinical data. |
| `step3_risk_assessment.py`    | Computes genetic risk scores based on annotated SNPs.         |

---

### ✅ Input  
- **Annotated VCF file** (processed using **Variant Effect Predictor (VEP)**).
- ### **🔹 Data Characteristics**
- **VCF files** contain **SNPs and structural variants** relevant to **hypertension and T2D**.
- **Variants were generated based on GWAS risk loci**, ensuring realistic mutation patterns.
- **The `.vcf.gz` format is used** for better performance in large-scale genomic processing.

- ## 📊 DeepCV Data Overview

DeepCV processes the following **synthetic VCF files** to simulate real-world variant annotation and risk assessment.

| **Dataset Name**        | **File Format** | **Description**                                  |
|------------------------|--------------|------------------------------------------------|
| `HTN_Synthetic.vcf`    | `.vcf`       | Synthetic VCF file for **hypertension (HTN)** genetic risk assessment. |
| `T2D_Synthetic.vcf`    | `.vcf`       | Synthetic VCF file for **type 2 diabetes (T2D)** genetic risk assessment. |
| `HTN_Synthetic.vcf.gz` | `.vcf.gz`    | Compressed version of `HTN_Synthetic.vcf` for optimized storage. |
| `T2D_Synthetic.vcf.gz` | `.vcf.gz`    | Compressed version of `T2D_Synthetic.vcf` for optimized storage. |

### ✅ Output  
- **List of genes with genetic risk scores** and their association with diseases saved as CSV file.
- **Filtered gene lists** based on disease-specific queries.
- 

# DeepCV Alpha - Genetic Risk Scoring Workflow  
The workflow has 4 main steps: vcf processing and filtring, anntionon, risk score assment, and risk graph constriction. The repostir provide code form step 1 to 3. 
in the follong disted of each step and how it wwill be exuted using DeepCV : To run **DeepCV Alpha**, use the following command:

![Image](https://github.com/user-attachments/assets/6733cd05-ea8f-4194-a533-12be4e781858)

# **DeepCV: Command Guide & Usage Instructions 🚀**

DeepCV is a bioinformatics pipeline that processes **VCF files**, performs **SNP annotation**, and conducts **genetic risk assessment**.

This guide provides a **step-by-step catalog** of all DeepCV commands.

---

##Setting Up DeepCV (First-Time Users)
Ensure dependencies are installed before running DeepCV:

```bash
pip install -r requirements.txt
```
Or create a Conda environment:

```bash
conda env create -f environment.yml
conda activate DeepCV
```

 

#  Example Workflows
## **1️⃣ Running the Full Pipeline**
To run **all steps** (VCF processing, SNP annotation, and risk assessment) in one command:

```bash
python DeepCV_Main.py --input_vcf input.vcf --output_path "deepcv_results.csv"
```
🔹 Full Pipeline with Single Disease
```bash
python DeepCV_Main.py --input_vcf sample.vcf --output_path "diabetes_risk.csv" --disease_name "Type 2 Diabetes Mellitus"
```

This command will:
- Process input.vcf
- Annotate SNPs using Ensembl VEP API
- Perform genetic risk scoring
-Save results to deepcv_results.csv

🔹 Annotate & Assess Multiple Diseases
```bash
python DeepCV_Main.py --skip_step1 --disease_name "Cardiovascular Disease Hypertension"

```
🔹 Only Compute Risk Scores for Annotated Data
```bash
python DeepCV_Main.py --skip_step1 --skip_step2 --output_path "risk_assessment.csv"

```


## **2️⃣ Step-by-Step Execution
🔹 Step 1: Process & Filter VCF File
Run only the VCF processing & filtering step:

```bash
python DeepCV_Main.py --input_vcf input.vcf --skip_step2 --skip_step3
```
✅ This will:
- Extract variant information from `input.vcf`
- Filter SNPs based on quality score & depth
- Save processed variants to: `results/processed_variants.csv`

🔹 Step 2: Annotate SNPs
If you already have a processed VCF file and want to run only SNP annotation:
```bash
python DeepCV_Main.py --skip_step1 --skip_step3
```
✅ This will:
Annotate variants with ClinVar, dbSNP, and GWAS Catalog
Save annotated results to: `results/final_annotated_output.csv`

🔹 Step 3: Perform Genetic Risk Assessment
If you have an annotated file and only want to compute genetic risk scores:
```batch
python DeepCV_Main.py --skip_step1 --skip_step2
```
✅ This will:
- Use Open Targets API to find disease associations
- Compute risk scores for genes
- Save results to: `results/final_risk_assessment.csv`

``

3️⃣ Specifying Disease for Risk Assessment
By default, DeepCV analyzes all disease associations, but you can filter for specific diseases.

🔸 Single Disease Analysis
```bash
python DeepCV_Main.py --input_vcf input.vcf --disease_name "Type 2 Diabetes Mellitus"
```
✅ This will: Only keep genes related to Type 2 Diabetes Mellitus in the final output.

🔸 Multiple Disease Analysis
To analyze multiple diseases, use space-separated disease names inside quotes:
```bash
python DeepCV_Main.py --disease_name "Type 2 Diabetes Mellitus Hypertension"

```
4️⃣ Output Files & Directories
DeepCV saves all results inside the results/ directory:

After execution, DeepCV generates several files:

| **Filename**                   | **Description**                                             |
|--------------------------------|-------------------------------------------------------------|
| `final_annotated_output.csv`   | Annotated SNPs with functional & clinical significance     |
| `final_risk_assessment.csv`    | Computed genetic risk scores for all genes                 |
| `filtered_by_disease.csv`      | Disease-specific filtered results                          |
| `filtered_by_disease_cleaned.csv` | Duplicate-removed version of `filtered_by_disease.csv`  |


---


## 🔹 DeepCV Alpha Workflow Steps  
1️⃣ **Input Processing Filtering:**  (**`DeepCV_Main.py`** ) **`step1_processing.py`** extracts relevant variants from **VCF files**.
   - Reads **VEP-annotated VCF files**.  
   - Extracts **genomic variants** and their **functional impact**.

```bash
python DeepCV_Main.py input.vcf
```

2️⃣ **Variant Annotation**    **`step2_SNP_annotation.py`** adds **functional & clinical annotations**.  
   - Integrates **GWAS catalog mutations**.  
   - Assigns **clinical significance scores** to each variant.  
   - Computes **weighted genetic risk scores** for each gene.

```bash
python DeepCV_Main.py input.vcf --output_path "deepcv_results.csv"
```

3️⃣ **Risk Calculation:**  **`step3_risk_assessment.py`** computes **risk scores** based on **gene-disease associations**.  
   - Produces a **list of genes with associated risk scores**.  
   - Allows **filtering genes** based on disease-specific associations.
  
4️⃣ The final risk assessment results are saved in **CSV files** for further analysis.

## 🔹 Filtering Genes by Disease  

To filter the **genetic risk scores** based on a specific disease, use:
```bash
python DeepCV_Main.py input.vcf --disease_name "Your Disease Name"

```

---

📩 **For more details, read our upcoming research paper!**  
📩 **Contact:** [sakhaa.alsaedi@kaust.edu.sa](mailto:sakhaa.alsaedi@kaust.edu.sa)  



