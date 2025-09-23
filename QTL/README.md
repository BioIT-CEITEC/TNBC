# Matrix pQTL Analysis: CNV and SNP-Based Cis/Trans pQTL Mapping

This repository contains two R scripts for performing cis/trans pQTL analysis using the [MatrixEQTL](https://cran.r-project.org/web/packages/MatrixEQTL/index.html) package. The analysis integrates either **copy number variation (CNV)** or **single nucleotide polymorphism (SNP)** data with **log₂-transformed protein intensities**, stratified by patient groups.

---

##  Data Access

All required input files are available in the shared folder:  
🔗 [Google Drive: QTL Data](https://drive.google.com/drive/u/0/folders/1LGKPmp7TodcTyKqCAYKAAHWOxuf9mDX5)

---

##  Script Overview

### `CNV_CisTrans_eQTL.R`
- **Genotype input**: CNV data
- **Expression input**: Protein intensities (log₂-transformed)
- **Output folder**: `Res/CNV_CisTrans_LOG`

### `SNPS_CisTrans_eQTL_fdr05.R`
- **Genotype input**: SNP data
- **Expression input**: Protein intensities (log₂-transformed)
- **Output folder**: `Res/SNPS_CisTrans_LOG_fdr05`

Both scripts:
- Filter genotype data to retain only variants with ≥5% representation across all three states (0, 1, 2)
- Log₂-transform protein intensities
- Match sample IDs across genotype and expression matrices
- Run MatrixEQTL for each patient group (and for all samples combined)
- Export summary tables and full significant results

---

##  Required Inputs


| File | Description |
|------|-------------|
| `CNV_tab.txt` / `SNP_may25_noNA.txt` | Raw genotype matrix (CNV or SNP) |
| `CNV_location.txt` / `SNP_location.txt` | Genomic coordinates for variants |
| `prot_CNV.txt` / `prot2_log2.txt` | Protein intensity matrix |
| `prot_location_cnv.txt` / `prot_location_snps.txt` | Genomic coordinates for proteins |
| `Pacient_names.Bouchal.xlsx` | Sample-to-patient mapping |
| `full_clust_tab.tsv` | Patient group assignments |

---

##  Outputs

Each script generates the following per group (including `"all"`):

### 🔹 QQ Plots
- `group_<name>_qqplot.png`  
  Visualizes p-value distribution for cis/trans associations.

### 🔹 Summary Tables
- `group_<name>_summary.csv`  
  Aggregated counts of significant associations per protein and variant type.

### 🔹 Full Significant Results
- `group_<name>_signif.csv`  
  All significant cis/trans eQTLs with FDR below threshold.

All outputs are saved in:
- `Res/CNV_CisTrans_LOG/` for CNV script
- `Res/SNPS_CisTrans_LOG_fdr05/` for SNP script

---

##  How to Run

1. Open R or RStudio.
2. Set your working directory to the project folder.
3. Run either script:
   ```r
   source("CNV_CisTrans_eQTL.R")
   # or
   source("SNPS_CisTrans_eQTL_fdr05.R")
