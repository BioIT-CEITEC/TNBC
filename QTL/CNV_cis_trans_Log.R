# CIS, TRANS 
# eQTL analysis by patient group using MatrixEQTL 
# -----------------------------------------------
rm(list=ls())

# Load required libraries
library(MatrixEQTL)
library(readxl)
library(dplyr)
library(data.table)

#### INPUTS -------
base.dir        <- "QTL"
cnv_file        <- file.path(base.dir, "CNV_tab.txt")
prot_file       <- file.path(base.dir, "prot_CNV.txt")
cnvpos_file     <- file.path(base.dir, "CNV_location.txt")
protpos_file    <- file.path(base.dir, "prot_location_cnv.txt")
patient_map_xls <- "Pacient_names.xlsx"
group_tab_file  <- "full_clust_tab.tsv"
output_dir      <- file.path(base.dir, "Res/CNV_CisTrans_LOG")

# Set working directory and create output folder
setwd(base.dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#------------------------------------------------
# 1. Read in genotype and expression data
#------------------------------------------------
CNV <- fread(cnv_file)

# Define filtering function
filter_cnv <- function(row) {
  counts <- table(factor(row, levels = c(0, 1, 2)))
  percentages <- counts / length(row)
  all(percentages > 0.05)
}

# Apply row-wise filtering
filtered_cnv <- CNV[apply(CNV, 1, filter_cnv), ]

# Save filtered CNV data
filtered_cnv_file <- file.path(base.dir, "CNV_tab_filtered.txt")
write.table(filtered_cnv, filtered_cnv_file, sep = "\t", row.names = FALSE, quote = TRUE)

snps <- SlicedData$new()
snps$fileDelimiter      <- "\t"
snps$fileOmitCharacters <- "NA"
snps$fileSkipRows       <- 1
snps$fileSkipColumns    <- 1
snps$fileSliceSize      <- 2000
snps$LoadFile(filtered_cnv_file)

# Log2 transform protein data
prot <- fread(prot_file)
prot[, 2:ncol(prot)] <- log2(prot[, 2:ncol(prot)] + 1)
prot_log_file <- file.path(base.dir, "prot_CNV_log2.txt")
write.table(prot, prot_log_file, sep = "\t", row.names = FALSE, quote = TRUE)

gene <- SlicedData$new()
gene$fileDelimiter      <- "\t"
gene$fileOmitCharacters <- "NA"
gene$fileSkipRows       <- 1
gene$fileSkipColumns    <- 1
gene$fileSliceSize      <- 2000
gene$LoadFile(prot_log_file)

cvrt <- SlicedData$new()

# Read positions
cnvpos <- read.table(cnvpos_file, header = TRUE, stringsAsFactors = FALSE) %>%
  filter(snpid %in% filtered_cnv$V1)
protpos <- read.table(protpos_file, header = TRUE, stringsAsFactors = FALSE)

# ------------------------------------------------
# 2. Filter out genes (proteins) without position
# ------------------------------------------------
common_genes <- intersect(rownames(gene), protpos$id)
print(paste("Common genes found:", length(common_genes)))
gene$RowReorder(match(common_genes, rownames(gene)))
protpos <- protpos %>% filter(id %in% common_genes)
stopifnot(all(rownames(gene) %in% protpos$id))

# Filter SNPs to only those with known positions
common_snps <- intersect(rownames(snps), cnvpos$snpid)
print(paste("Common snps found:", length(common_snps)))
snps$RowReorder(match(common_snps, rownames(snps)))
cnvpos <- cnvpos %>% filter(snpid %in% common_snps)
stopifnot(nrow(snps) == nrow(cnvpos))
stopifnot(all(rownames(snps) %in% cnvpos$snp))

#------------------------------------------------
# 3. Read patient-to-sample mapping
#------------------------------------------------
patient_map <- read_excel(patient_map_xls) %>%
  select(ID = `sample name`, sample = PROTEIN, tumor) %>%
  mutate(
    sample = sprintf("P_%03d", sample),
    ID     = if_else(is.na(ID), sub("_T", "", tumor), sub("_N", "", ID))
  ) %>%
  select(ID, sample)

#------------------------------------------------
# 4. Read group assignments
#------------------------------------------------
group_tab <- fread(group_tab_file) %>%
  select(sample, PROT_clust) %>%
  left_join(patient_map, by = "sample") %>%
  select(ID, PROT_clust) %>%
  distinct() %>%
  rstatix::drop_na(PROT_clust)

#------------------------------------------------
# 5. Check for sample mismatches
#------------------------------------------------
cn_snps  <- gsub('"', '', snps$columnNames)
cn_genes <- gsub('"', '', gene$columnNames)
samples_all <- group_tab$ID
not_in_snps  <- setdiff(samples_all, cn_snps)
not_in_genes <- setdiff(samples_all, cn_genes)
cat(length(not_in_snps)," Samples not found in SNP columns:\n", sort(not_in_snps), sep = " ", fill = TRUE)
cat(length(not_in_genes),",Samples not found in gene columns:\n", sort(not_in_genes), sep = " ", fill = TRUE)

#------------------------------------------------
# 6. eQTL engine settings
#------------------------------------------------
useModel          <- modelLINEAR
pvOutputThreshold <- 1e-2
errorCovariance   <- numeric()
cisDist           <- 1e6
fdr_threshold     <- 0.01

#------------------------------------------------
# 7. Loop over groups, run eQTL, export results
#------------------------------------------------
cn <- gsub('"','', snps$columnNames)
all_idx <- seq_along(cn)

for (grp in c("all", unique(group_tab$PROT_clust))) {
  cat("\nRunning group:", grp, "\n")
  
  if (grp=="all") {
    snp_idx <- all_idx
  } else {
    samp    <- group_tab %>% filter(PROT_clust==grp) %>% pull(ID)
    snp_idx <- which(cn %in% samp)
  }
  # now *always* clone & subsample:
  # determine SNP columns matching samp, preserving original order
  snps_grp    <- snps$Clone()
  snps_grp$ColumnSubsample(snp_idx)
  grp_snps    <- gsub('"', '', snps_grp$columnNames)
  
  gene_idx    <- match(grp_snps, cn_genes)
  gene_grp    <- gene$Clone()
  gene_grp$ColumnSubsample(gene_idx)
  
  all(sort(snps_grp$columnNames)==sort(gene_grp$columnNames))
  
  # Check matching sample names
  if (!all(sort(snps_grp$columnNames) == sort(gene_grp$columnNames))) {
    stop("Mismatch between SNP and gene column names for group ", grp)
  }
  cat(grp,"script sees", length(snps_grp$columnNames), "samples\n")
  
  # Run MatrixEQTL
  me_gr <- Matrix_eQTL_main(
    snps                = snps_grp,
    gene                = gene_grp,
    cvrt                = cvrt,
    useModel            = useModel,
    errorCovariance     = errorCovariance,
    verbose             = TRUE,
    snpspos             = cnvpos,
    genepos             = protpos,
    cisDist             = cisDist,
    pvOutputThreshold   = 1e-4,
    pvOutputThreshold.cis = 1e-4,
    output_file_name    = tempfile(),
    output_file_name.cis = tempfile(),
    pvalue.hist         = ifelse(length(snps_grp$columnNames) >5,"qqplot",FALSE),
    min.pv.by.genesnp   = FALSE,
    noFDRsaveMemory     = FALSE
  )
  
  # Filter results
  png(file.path("Res/CNV_CisTrans_LOG", paste0("group_", grp, "_qqplot.png")), width = 800, height = 600)
  if (length(snps_grp$columnNames) >5 ){plot(me_gr, pch = 16, cex = 0.7)}
  dev.off()
  
  # --- 8. Filter, label and save cis/trans results ---
  # Set FDR threshold
  fdr_threshold <- 0.01
  
  # Extract cis and trans eQTLs
  cis_eqtls   <- me_gr$cis$eqtls %>%filter(FDR < fdr_threshold)%>% mutate(Type = "cis")
  trans_eqtls <- me_gr$trans$eqtls %>%filter(FDR < fdr_threshold)%>% mutate(Type = "trans")
  
  # Combine and filter significant
  eqtls_all <- bind_rows(cis_eqtls, trans_eqtls) %>%
    arrange(FDR) %>%
    filter(FDR < fdr_threshold)
  
  
  print(head(eqtls_all, 10))
  # Save to environment
  assign(paste0("eQTLs_", grp), eqtls_all)
  
  # --- Summary table ---
  sum_fn <- file.path("Res/CNV_CisTrans_LOG", paste0("group_", grp, "_summary.csv"))
  eqtls_all %>%
    group_by(gene, Type) %>%
    summarise(n = n(), CNVs = paste(sort(snps), collapse = ";"), .groups = "drop") %>%
    rename(Protein = gene) %>%
    arrange(desc(n)) %>%
    write.csv(sum_fn, row.names = FALSE)
  
  # --- Full significant results ---
  full_fn <- file.path("Res/CNV_CisTrans_LOG", paste0("group_", grp, "_signif.csv"))
  eqtls_all %>%
    arrange(FDR) %>%
    rename(Protein = gene, CNV = snps) %>%
    write.table(full_fn, sep = ";", row.names = FALSE, quote = FALSE)
  
  message("Group ", grp, ": done in ", me_gr$time.in.sec, "s — ",
          nrow(eqtls_all), " significant eQTLs (", sum(eqtls_all$Type == "cis"), " cis, ",
          sum(eqtls_all$Type == "trans"), " trans)")
}

#------------------------------------------------
# 8. Print group summary
#------------------------------------------------
for (grp in c("all", unique(group_tab$PROT_clust))) {
  cat("\nGroup:", grp, "with n =",
      if (grp == "all") length(all_idx) else length(which(group_tab$PROT_clust == grp)), "samples\n")
  se_gr <- get(paste("eQTLs", grp, sep = "_"))
  cat("Significant eQTLs:", nrow(se_gr), "\n")
}

