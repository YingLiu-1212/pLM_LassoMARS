# Drug Response Prediction and Clinical Analysis Pipeline

A comprehensive machine learning framework for predicting amino acid mutation effects on universal stress response using high-throughput drug screening, protein language models (pLM), and amino acid topology attributes.

## Overview

This repository contains code for training and validating machine learning models that predict how amino acid mutations affect cellular responses to targeted cancer therapies. The framework integrates structural biology, pLM embeddings, and drug screen data to enable accurate prediction of drug response mechanisms.

## Key Features

- **Multi-modal Feature Integration**: Combines protein language model embeddings, structural topology features, and biophysical properties
- **Ensemble Modeling**: Uses Lasso-MARS ensemble approach for robust prediction
- **Cross-Validation**: Extensive validation across multiple drug classes and clinical datasets
- **Clinical Translation**: Applications to TCGA and MSK cancer genomics data
- **Structural Insights**: Analysis of mutation-phosphorylation site spatial relationships

## Model Architecture

### Input Features
- **Protein Embeddings**: 1024-dimensional embeddings from protein language models (ProtT5)
- **Structural Topology**: Contact flexibility (CF) and local density (LD) metrics from AlphaFold2 structures
- **Biophysical Properties**: Amino acid mass, isoelectric point, hydrophobicity, and structural angles

### Machine Learning Approach
- **Feature Selection**: LASSO regularization for informative feature identification
- **Non-linear Modeling**: Multivariate Adaptive Regression Splines (MARS) for complex relationships
- **Ensemble Prediction**: Combined predictions from multiple drug-specific models

## Installation

### 1.1 Operating Systems
- **Linux** 
- Tested on: Rocky Linux 8.7 

### 1.2 Software Dependencies

#### R Environment
- **R version**: 4.1.3

#### R Packages with Version Requirements:
```r
# Core packages (with tested versions)
dplyr (>= 1.1.4)
tidyr (>= 1.3.1)
ggplot2 (>= 3.5.1)
stringr (>= 1.5.0)
matrixStats (>= 1.4.1)
VennDiagram (>= 1.7.3)
pROC (>= 1.18.5)  # for ROC analysis
survival (>= 3.7.0)  # for survival analysis
survminer (>= 0.4.9)  # for survival plots
earth (>= 5.3.4)  # for MARS modeling
glmnet (>= 4.1.8)  # for LASSO regression
```

#### Additional Requirements:
- **Disk Space**: Minimum 120GB (for ProtT5 embeddings)
- **Memory**: 64GB RAM recommended
- **Processor**: Multi-core CPU (8+ cores recommended)

### R Dependencies
```r
install.packages(c(
  "dplyr", "tidyr", "stringr", "earth", "glmnet", "pROC", 
  "Rfast", "survival", "survminer", "ggplot2", "ggpubr",
  "RColorBrewer", "VennDiagram", "bio3d"
))
```

## Project Structure

### Core Analysis Pipeline
```
├── 0.similar_protein_by_cosine.R      # Protein similarity analysis
├── 1.0.training_prepare.R             # Data preprocessing and feature engineering
├── 2.0.model_construct.R              # Model training and evaluation
├── 3.0.cross_drug_validation.R        # Cross-drug validation
├── 4.1.TCGA_mut_proc.R               # TCGA clinical data processing
├── 4.2.TCGA_predict.R                # TCGA prediction and survival analysis
├── 4.3.MSK_mut_proc.R                # MSK clinical data processing  
├── 4.4.MSK_predict.R                 # MSK prediction and survival analysis
├── 5.1.plot_Predict_tcga_msk_PTM.R   # Phosphorylation site analysis
├── 5.2.plot_cosmic_AB_PTM.R          # COSMIC resistance mutation analysis
├── 5.3.clinical_core_surv_Phos_PDB.R # Clinical survival analysis
└── function_pLM_LassoMARS.R          # Core modeling functions
```

## Installation Guide

### Step 1: Install R and Required Packages
```bash

# Install required packages
Rscript -e "install.packages(c('dplyr', 'tidyr', 'ggplot2', 'stringr', 'matrixStats', 'VennDiagram', 'RColorBrewer', 'pROC', 'survival', 'survminer', 'earth', 'glmnet'), dependencies=TRUE)"
```

### Step 2: Clone Repository and Set Up Directory Structure
```bash
# Clone repository
git clone https://github.com/YingLiu-1212/pLM_LassoMARS.git

```
### Installation Time
- **Data download**: 1-2 hours (depending on internet speed and file sizes)
- **Total setup**: 1.5-3 hours for complete setup

## Usage

### Step 1: Data Preparation
```bash
# Process training data for specific drugs
Rscript 1.0.training_prepare.R
```
#### Expected output:
- training_sites_abem.txt, training_sites_bini.txt, training_sites_olap.txt
- KO_phenotype_[drug].txt files
- training_site_CF_[drug].txt files
- mars_train_[drug]_info.txt files


### Step 2: Model Training
```bash
# Train Lasso-MARS ensemble models
Rscript 2.0.model_construct.R
```

#### Model Files (c.training)
1. **RData files**: Saved models pLM_LassoMARS.RData
2. **CSV files**: Model coefficients for interpretability

#### Visualization Output (o.output_figures/)
1. **Model component plots**: Show MARS basis functions and relationships
2. **Combined ROC curve**: Compares performance across all three models
3. **AUC values**: Quantitative performance metrics


#### Feature Space
Models are trained on:
- Protein embeddings (1024 dimensions)
- Amino acid embeddings (1024 dimensions)
- Structural topology features
- Amino acid biophysical properties (nAngels, Mass, IP, HM)
- Knockout phenotype information

#### Performance Metrics
- **AUC (Area Under ROC Curve)**: Primary performance metric
- **ROC curves**: Visualization of true positive vs. false positive rates
- **Model coefficients**: Interpretability of feature importance

### Step 3: Cross-Drug Validation
```bash
# Validate models on independent drug datasets
Rscript 3.0.cross_drug_validation.R
```
#### Supporting Data Files (a.data/)
**`AF2_WG_CF.RData`** - AlphaFold2 structural features and confidence metrics
**`embeddings_per_protein.txt`** - Protein embeddings from protein language models
**`cross_drug_eb_per_residue_all.txt`** - Amino acid residue embeddings
**`aa_property.txt`** - Amino acid biophysical properties
**`Prot_scope_common.csv`** - Core protein subset for focused analysis
**`pLM_LassoMARS.RData`** - Pre-trained machine learning models

#### Processing Files (d.validation_cross_drug/)
1. **`cross_drug_screen_phenotype.csv`** - Processed drug resistance phenotypes
2. **`cross_drug_edit_AA_canonical.txt`** - Canonical amino acid editing information
3. **`cross_drug_edit_AA_phenotype_AF2CF.txt`** - Phenotype data with AlphaFold2 structural features
4. **`cross_drug_full_phenotype_info.txt`** - Complete phenotype-structure dataset for all proteins
5. **`cross_drug_Core_phenotype_info.txt`** - Phenotype-structure dataset for core protein subset

#### Model Input Files (d.validation_cross_drug/)
6. **`Drug_[drugname]_full_test_info.txt`** - Test data information for full protein sets
7. **`Drug_[drugname]_full_eb_all.txt`** - Embedding data for full protein sets
8. **`Drug_[drugname]_Core_test_info.txt`** - Test data information for core protein subsets
9. **`Drug_[drugname]_Core_eb_all.txt`** - Embedding data for core protein subsets


### Step 4: Clinical Prediction
```bash
# Apply models to clinical datasets
Rscript 4.2.TCGA_predict.R
Rscript 4.4.MSK_predict.R
```
#### Expected output:
- `tcga_core73_prediction_info.txt` - Site-level prediction results
- `tcga_patient_stress_prediction.txt` - Patient-level prediction results
- `Surv_OS_StressResponse_tcga.pdf` - Overall survival curves
- `Surv_OS_StressResponse_tcga_type.pdf` - Survival curves stratified by cancer type
- `msk_core73_prediction_info.txt` - Site-level prediction results
- `msk_patient_stress_prediction.txt` - Patient-level prediction results
- `Surv_OS_StressResponse_msk.pdf` - Overall survival curves
- `Surv_OS_StressResponse_msk_type.pdf` - Survival curves stratified by cancer type


### Step 5: COSMIC and Screening Data Phosphorylation Analysis 

```bash
Rscript 5.1.plot_Predict_tcga_msk_PTM.R
Rscript 5.2.plot_cosmic_AB_PTM.R
Rscript 5.3.clinical_core_surv_Phos_PDB.R  # Run clinical survival analysis by phosphorylation proximity script
```
#### Expected output:
- `PDB_AA_Unires_clinical_prediction_boxplot.pdf` - Comparative boxplot visualization
- `Surv_OS_PhosDist3_clinical_COSMIC_core.pdf` - Survival curves based on phosphorylation distance


## Reproduction Instructions
- 1. Download complete dataset
- 2. Extract all compressed directory
- 3. Run scripts in order

### Expected Runtime
- **Data preparation**: ~2 hours
- **Script execution**: ~12 hours


## License
This project is licensed under the MIT License - see the LICENSE file for details.