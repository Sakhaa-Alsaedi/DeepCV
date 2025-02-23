# 🧬 DeepCV Alpha - Genetic Risk Scoring Workflow  

## 🔹 Overview
The **DeepCV Alpha Version** is a validated workflow for **computational variant annotation and genetic risk scoring**. This pipeline has been tested using **synthetic data from `NICP.vcf`** and **GWAS catalog mutations**.

### ✅ Purpose  
- Computes **genetic risk scores** for **each gene associated with a disease**.  
- Provides **a filterable gene list** for disease-specific risk analysis.  

### ✅ Input  
- **Annotated VCF file** (processed using **Variant Effect Predictor (VEP)**).  

### ✅ Output  
- **List of genes with genetic risk scores** and their association with diseases.  
- **Filtered gene lists** based on disease-specific queries.

---

## 🔹 DeepCV Alpha Workflow Steps  
1️⃣ **Input Processing:**  
   - Reads **VEP-annotated VCF files**.  
   - Extracts **genomic variants** and their **functional impact**.

2️⃣ **Variant Annotation & Risk Calculation:**  
   - Integrates **GWAS catalog mutations**.  
   - Assigns **clinical significance scores** to each variant.  
   - Computes **weighted genetic risk scores** for each gene.  

3️⃣ **Filtering & Output Generation:**  
   - Produces a **list of genes with associated risk scores**.  
   - Allows **filtering genes** based on disease-specific associations.  

---

## 🔹 Example Execution  
To run **DeepCV Alpha**, use the following command:

```bash
python DeepCV_Main.py input.vcf --output_path "deepcv_results.csv"
```

## 📂 DeepCV Workflow - Files & Descriptions

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

📌 **How the Workflow Works:**
1️⃣ **`DeepCV_Main.py`** calls **Step 1, 2, and 3** in sequence.  
2️⃣ **`step1_processing.py`** extracts relevant variants from **VCF files**.  
3️⃣ **`step2_SNP_annotation.py`** adds **functional & clinical annotations**.  
4️⃣ **`step3_risk_assessment.py`** computes **risk scores** based on **gene-disease associations**.  
5️⃣ The final risk assessment results are saved in **CSV files** for further analysis.

🚀 **Run DeepCV with:**
```bash
python DeepCV_Main.py input.vcf
```

### 📁 Output Files

After execution, DeepCV generates several files:

| **Filename**                   | **Description**                                             |
|--------------------------------|-------------------------------------------------------------|
| `final_annotated_output.csv`   | Annotated SNPs with functional & clinical significance     |
| `final_risk_assessment.csv`    | Computed genetic risk scores for all genes                 |
| `filtered_by_disease.csv`      | Disease-specific filtered results                          |
| `filtered_by_disease_cleaned.csv` | Duplicate-removed version of `filtered_by_disease.csv`  |

