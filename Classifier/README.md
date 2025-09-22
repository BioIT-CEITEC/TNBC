#  TNBC Proteomic Classification Pipeline

This repository contains a complete R-based pipeline for classifying triple-negative breast cancer (TNBC) subtypes using proteomic data. It includes dimensionality reduction, feature selection, model training, hyperparameter tuning, and external validation.

All steps are modular and reproducible using the provided R scripts. No prior knowledge of R or proteomics is required—just follow the instructions below.

**Repository:** [BioIT-CEITEC/TNBC](https://github.com/BioIT-CEITEC/TNBC.git)  
**Data folder:** [ TNBC Proteomic Data](https://drive.google.com/drive/u/0/folders/1BMSASG2U_DUofrV2QQpVWJmudhL63LXg)

---

##  Repository Contents

| Script Name                          | Purpose                                                                 |
|-------------------------------------|-------------------------------------------------------------------------|
| `TNBC_separability.R`               | Visualizes sample clustering and separability using t-SNE, UMAP, and LDA |
| `Feature_selection_TNBC.R`          | Selects the most informative proteins using multiple feature selection methods |
| `Bayes_interference_tuning.R`       | Optimizes Random Forest hyperparameters using Bayesian optimization     |
| `Random_forest_TNBC.R`              | Trains and evaluates the final Random Forest classification model       |
| `final_model_TNBC.rds`              | Pre-trained Random Forest model          |

---

##  Step-by-Step Workflow

### 1. Sample Visualization and Separability (`TNBC_separability.R`)

**Goal**: Assess whether patient samples naturally cluster based on protein expression.

**Methods**:
- t-SNE and UMAP for dimensionality reduction and visualization.
- LDA to evaluate linear separability of predefined groups.

---

### 2. Feature Selection (`Feature selection TNBC.R`)

**Goal**: Identify the most informative proteins for classification.

**Methods**:
- Remove highly correlated and low-variance proteins.
- Apply Boruta and Recursive Feature Elimination (RFE) using Random Forest.

---

### 3. Hyperparameter Optimization (`Bayes_interference_tuning.R`)

**Goal**: Tune Random Forest parameters to maximize classification accuracy.

**Methods**:
- Random Search over 100 combinations.
- Bayesian Optimization to iteratively refine parameter selection.

---

### 4. Model Training and Evaluation (`Random_forest_TNBC.R`)

**Goal**: Train a Random Forest classifier and evaluate its performance.

**Methods**:
- Hierarchical clustering for group consistency.
- Leave-One-Out Cross-Validation for robust accuracy estimation.



##  Script Inputs

Each script begins with a clearly defined `# INPUTS` block specifying required file paths and objects. You must update these paths to match your local setup or use the provided data folder.

**Example:**
```r
#### INPUTS ----------
dir_im <- "TNBC/Classifier"
model_im <- "final_model_TNBC.rds"
