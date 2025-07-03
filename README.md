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

Below is an overview of the **DeepCV workflow** and the corresponding scripts.

| **Filename**                   | **Description**                                                 |
|--------------------------------|-----------------------------------------------------------------|
| `DeepCV_Main.py`              | The main script that runs the entire pipeline.                 |
| `environment.yml`             | Conda environment configuration file for dependency management. |
| `setup_deepcv.py`             | Automates setup, environment creation, and installation.       |
| `step1_processing.py`         | Parses and filters VCF files for variant processing.           |
| `step2_SNP_annotation.py`     | Annotates SNPs using functional impact scores and clinical data. |
| `step3_risk_assessment.py`    | Computes genetic risk scores based on annotated SNPs.         |
| `step4_riskGraph.py`    | Generates genetic risk graphs based gene-disease networks.         |

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
| Step4 output: `Risk_Graph.png`      | Risk graph results                          |



## DeepCV Alpha - Genetic Risk Scoring Workflow  
The workflow consists of four main steps: VCF processing and filtering, annotation, risk score assessment, and risk graph construction.
The repository provides code for steps 1 to 4.Below is a detailed description of each step and how it will be executed using DeepCV:

![Image](https://github.com/user-attachments/assets/ce343b15-477a-4d99-9bad-6acf0692988e)

### Key Features

- ** VCF Parsing & Quality Control** - Comprehensive variant filtering and quality assessment
- ** SNP Annotation** - VEP API integration with intelligent caching system
- **⚕ Risk Assessment** - Disease-gene association analysis with risk scoring
- ** Network Visualization** - STRING protein-protein interaction networks
- ** Comprehensive Reporting** - Step-by-step visualizations and PDF reports
- ** Performance Optimized** - Caching system for faster repeated analyses

### **DeepCV: Command Guide & Usage Instructions**

- To run DeepCV Alpha, use the following command:

## 1. Setting Up DeepCV (First-Time Users)
Ensure dependencies are installed before running DeepCV:

- create a Conda environment:

```bash
conda env create -f environment.yml
conda activate DeepCV
```
- or run create setup_deepcv.py:

```bash
python setup_deepcv.py --input_vcf ./path/input.vcf
conda activate DeepCV
```

## 2. Executing DeepCV Alpha Workflow Steps

🔹 Full Pipeline with All Diseases
To run **all steps** (VCF processing, SNP annotation, and risk assessment) in one command:

```bash
# Run complete pipeline with all features
python DeepCV_Main_Fixed.py --input_vcf your_file.vcf --output_dir results

# Example with Demo data
python DeepCV_Main.py --input_vcf Demo.vcf --output_dir results
```

### Advanced Usage

```bash
# Custom STRING confidence threshold
python DeepCV_Main.py --input_vcf your_file.vcf --output_dir results --string_confidence 0.7

# Skip specific components
python DeepCV_Main.py --input_vcf your_file.vcf --output_dir results --skip_pdf
python DeepCV_Main.py --input_vcf your_file.vcf --output_dir results --skip_step_viz
python DeepCV_Main.py --input_vcf your_file.vcf --output_dir results --skip_visualization
```

## 📋 Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input_vcf` | Input VCF file path (required) | - |
| `--output_dir` | Output directory | `results` |
| `--string_confidence` | STRING confidence threshold (0.0-1.0) | `0.7` |
| `--disease_name` | Filter by specific disease name | `None` |
| `--skip_pdf` | Skip PDF report generation | `False` |
| `--skip_step_viz` | Skip step-by-step visualizations | `False` |
| `--skip_visualization` | Skip network analysis | `False` |
| `--help` | Show help message | - |

## 📊 Pipeline Steps

### Step 1: VCF Parsing & Quality Control
- Loads and filters VCF files based on quality thresholds
- Removes duplicate variants
- Generates quality distribution plots

### Step 2: SNP Annotation
- Annotates variants using VEP (Variant Effect Predictor) API
- Intelligent caching system for improved performance
- Adds functional consequence and impact information

### Step 3: Risk Assessment
- Calculates gene-based risk scores
- Fetches disease associations from Open Targets
- Generates comprehensive gene-disease association matrix

### Step 4: Network Visualization
- Creates protein-protein interaction networks using STRING database
- Generates comprehensive disease-gene-PPI networks
- Produces high-quality publication-ready visualizations

## 📁 Output Structure

```
results/
├── 📊 Core Analysis Files
│   ├── processed_variants_[filename].csv          # Step 1: Filtered variants
│   ├── final_annotated_output_[filename].csv      # Step 2: Annotated variants
│   └── final_risk_assessment_[filename].csv       # Step 3: Disease associations
├── 📈 Step-by-Step Visualizations
│   ├── step1_quality_distributions_[filename].png
│   ├── step1_variant_types_[filename].png
│   ├── step2_annotation_summary_[filename].png
│   ├── step2_impact_distribution_[filename].png
│   ├── step2_annotation_coverage_[filename].png
│   ├── step3_risk_distribution_[filename].png
│   ├── step3_top_genes_[filename].png
│   ├── step3_disease_analysis_[filename].png
│   └── step3_gene_disease_heatmap_[filename].png
├── 🕸️ Network Analysis
│   ├── comprehensive_network.png
│   ├── disease_gene_ppi_network.png
│   └── ppi_source_target_interactions.csv
└── 📄 Comprehensive Report
    └── DeepCV_Analysis_Report_[filename].pdf
```

## 🔬 Example Analysis Results

Using the provided `Demo.vcf` file:

### Input Data
- **Total variants in VCF**: 7,169
- **Quality-filtered variants**: 1,436
- **Processing time**: ~3 minutes

### Key Findings
- **Disease associations identified**: 16,791
- **Unique genes with disease links**: 913
- **Unique diseases associated**: 2,836
- **Average risk score**: 0.355
- **Maximum risk score**: 0.897

### Top Risk Genes Identified
1. **FBN1** (score: 0.897) - Marfan syndrome
2. **SOS1** (score: 0.871) - Noonan syndrome
3. **SLC22A5** (score: 0.867) - Carnitine deficiency
4. **SPAST** (score: 0.867) - Spastic paraplegia
5. **GALC** (score: 0.859) - Krabbe disease

## 🛠️ Installation

### Option 1: Clone Repository
```bash
git clone https://github.com/Sakhaa-Alsaedi/DeepCV.git
cd DeepCV
pip install -r requirements.txt
```

### Option 2: Download Files
Download the following core files:
- `DeepCV_Main_Fixed.py` (main pipeline)
- `step1_processing_fixed.py`
- `step2_SNP_annotation_fixed.py`
- `step3_risk_assessment.py`
- `step4_riskGraph_updated.py`
- `deepcv_step_visualizations.py`
- `deepcv_pdf_report.py`

## 📋 Requirements

### Python Dependencies
```
pandas>=1.3.0
numpy>=1.21.0
requests>=2.25.0
matplotlib>=3.3.0
seaborn>=0.11.0
plotly>=5.0.0
networkx>=3.0.0
reportlab>=3.6.0
tqdm>=4.62.0
joblib>=1.1.0
cyvcf2>=0.30.0
```

### External APIs
- **VEP API**: For variant annotation (no API key required)
- **STRING API**: For protein-protein interactions (no API key required)
- **Open Targets API**: For disease associations (no API key required)

## 🎯 Performance Features

- **Intelligent Caching**: VEP annotations are cached for faster repeated analyses
- **Batch Processing**: API requests are batched for optimal performance
- **Parallel Processing**: Multi-threaded operations where applicable
- **Memory Optimization**: Efficient data handling for large VCF files

## 📊 Visualization Gallery

### Quality Control Plots
- Variant quality score distributions
- Variant type classifications
- Chromosome distribution analysis

### Annotation Analysis
- Impact severity distributions
- Functional consequence analysis
- Annotation coverage heatmaps

### Risk Assessment
- Gene risk score distributions
- Top risk genes analysis
- Disease category analysis
- Gene-disease association heatmaps

### Network Analysis
- Comprehensive gene-disease-PPI networks
- STRING protein interaction networks
- High-confidence interaction filtering


## 📚 Documentation

### File Descriptions
- `DeepCV_Main.py`: Main pipeline orchestrator
- `step1_processing.py`: VCF parsing and quality control
- `step2_SNP_annotation.py`: Variant annotation using VEP
- `step3_risk_assessment.py`: Disease risk assessment
- `step4_riskGraph.py`: Network visualization
- `deepcv_step_visualizations.py`: Step-by-step plotting
- `deepcv_pdf_report.py`: PDF report generation

### API Documentation
- [VEP REST API](https://rest.ensembl.org/)
- [STRING API](https://string-db.org/help/api/)
- [Open Targets API](https://platform-docs.opentargets.org/)

## 🤝 Contributing

We welcome contributions! Please feel free to submit issues, feature requests, or pull requests.

### Development Setup
```bash
git clone https://github.com/Sakhaa-Alsaedi/DeepCV.git
cd DeepCV
pip install -r requirements.txt
# Make your changes
# Test with Demo.vcf
python DeepCV_Main.py --input_vcf ./Data/Demo.vcf --output_dir test_results
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.



📩 **For more details, read our upcoming research paper!**  
📩 **Contact:** [sakhaa.alsaedi@kaust.edu.sa](mailto:sakhaa.alsaedi@kaust.edu.sa)  
