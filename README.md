# Universal Stress Response Prediction

A comprehensive machine learning framework for predicting amino acid mutation effects on universal stress response using high-throughput drug screening, protein language models, and amino acid topology attributes.

## Overview

This repository contains code for training and validating machine learning models that predict how amino acid mutations affect cellular responses to targeted cancer therapies. The framework integrates structural biology, deep learning embeddings, and functional genomics data to enable accurate prediction of drug resistance mechanisms.

## Key Features

- **Multi-modal Feature Integration**: Combines protein language model embeddings, structural topology features, and biophysical properties
- **Ensemble Modeling**: Uses Lasso-MARS ensemble approach for robust prediction
- **Cross-Validation**: Extensive validation across multiple drug classes and clinical datasets
- **Clinical Translation**: Applications to TCGA and MSK cancer genomics data
- **Structural Insights**: Analysis of mutation-phosphorylation site spatial relationships

## Model Architecture

### Input Features
- **Protein Embeddings**: 1024-dimensional embeddings from protein language models
- **Structural Topology**: Contact frequency (CF) and local density (LD) metrics from AlphaFold2 predictions
- **Biophysical Properties**: Amino acid mass, isoelectric point, hydrophobicity, and structural angles

### Machine Learning Approach
- **Feature Selection**: LASSO regularization for informative feature identification
- **Non-linear Modeling**: Multivariate Adaptive Regression Splines (MARS) for complex relationships
- **Ensemble Prediction**: Combined predictions from multiple drug-specific models

## Installation

### Prerequisites
- R (version 4.0+)
- Python (version 3.8+)
- PyMOL (for structural visualization)

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

### Data Requirements

#### Input Data
- **Drug Screening Data**: Base editor-based functional genomics screens
- **Protein Structures**: AlphaFold2 predictions 
- **Clinical Mutations**: TCGA and MSK sequencing data
- **Protein Embeddings**: ProtT5 model outputs
- **Amino Acid Properties**: Biophysical and chemical characteristics

#### Processed Features
- **Topology Metrics**: CF10, CF10QS, LD15, LD15RK, CF10QS_FPI
- **Embedding Vectors**: 1024-dimensional protein and residue embeddings
- **Clinical Annotations**: Survival outcomes, mutation types, treatment history

## Usage

### Step 1: Data Preparation
```bash
# Process training data for specific drugs
Rscript 1.0.training_prepare.R
```

### Step 2: Model Training
```bash
# Train Lasso-MARS ensemble models
Rscript 2.0.model_construct.R
```

### Step 3: Cross-Drug Validation
```bash
# Validate models on independent drug datasets
Rscript 3.0.cross_drug_validation.R
```

### Step 4: Clinical Prediction
```bash
# Apply models to clinical datasets
Rscript 4.2.TCGA_predict.R
Rscript 4.4.MSK_predict.R
```
