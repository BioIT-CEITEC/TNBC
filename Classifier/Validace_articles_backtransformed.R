# ARTICLES eval - 
# 1) SKONTROLOVAŤ, ČI JE 45 ZVOLENÝCH PROTEÍNOV DOSTUPNÝCH V DÁTACH 
# 2) VYROBIŤ SKUPINY PRE RANDOM FOREST  
# 3) SPOČÍTAŤ PROT_fc SKUPIN VS. VŠETCI PRE VŠETKY DOSTUPNÉ PROTEÍNOVÉ EXPRESIE

# LOAD LIBRARIES
library(dplyr)
library(data.table)
library(randomForest)
library(readxl)
library(tools)
library(tidyr)
library(ggplot2)
library(patchwork)

# SET WORKING DIRECTORY AND DATE
setwd("~/CEITEC/20241219_Bouchal_TNBC") 
wd <- getwd()
datum <- Sys.Date()

# LOAD TRAINING DATA
df <- readRDS("~/CEITEC/20241219_Bouchal_TNBC/Res/Filtered_Features_2025-02-27.rds")
data <- fread("~/CEITEC/20241219_Bouchal_TNBC/Katka_files/results_2025/full_clust_tab.tsv")
head(data)
head(df)

# IDENTIFY PROTEINS FROM TRAINING DATA (USING GENE IDs)
proteiny <- colnames(df)[grep("ENSG", colnames(df))]
length(proteiny)
proteiny[proteiny %in% data$gene_id]
length(proteiny[proteiny %in% data$gene_id])

protein <- data %>% 
  filter(gene_id %in% proteiny) %>% 
  select(gene_id, gene_name) %>% 
  distinct()

# READ IN ALL ARTICLES
folder_path <- "C:/Users/FNOadmin/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Articles"
files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
data_list <- lapply(files, read_excel)
names(data_list) <- sub("_proteinmatrix", "", file_path_sans_ext(basename(files)))

# FUNCTION TO DETERMINE BEST CANDIDATE COLUMN FOR PROTEIN IDENTIFIERS
get_best_identifier_column <- function(df, protein) {
  best_col <- NULL
  best_overlap <- 0
  prot_ids <- tolower(protein$gene_id)
  prot_names <- tolower(protein$gene_name)
  
  for (col in names(df)) {
    vals <- tolower(as.character(df[[col]]))
    overlap_ids <- intersect(vals, prot_ids)
    overlap_names <- intersect(vals, prot_names)
    overlap <- length(unique(c(overlap_ids, overlap_names)))
    
    if (overlap > best_overlap) {
      best_overlap <- overlap
      best_col <- col
    }
  }
  return(list(column = best_col, overlap = best_overlap))
}

# ITERATE THROUGH EACH ARTICLE TO CHECK IF IT CONTAINS ALL 45 PROTEINS
articles_with_all_proteins <- c()
for (article_name in names(data_list)) {
  df_article <- data_list[[article_name]]
  candidate <- get_best_identifier_column(df_article, protein)
  best_col <- candidate$column
  overlap_count <- candidate$overlap
  
  cat("Article:", article_name, "\n")
  if (!is.null(best_col)) {
    cat("  Best identifier column:", best_col, "\n")
    cat("  Overlap count:", overlap_count, "\n")
    if (overlap_count == nrow(protein)) {
      articles_with_all_proteins <- c(articles_with_all_proteins, article_name)
    }
  } else {
    cat("  No candidate identifier column found.\n")
  }
  cat("\n")
}
cat("Articles that contain all 45 proteins:\n")
print(articles_with_all_proteins)

# SELECT EXTERNAL DATASET (THE WINNER, e.g. "Anurag_CancerDiscov_2022")
df_test <- read_excel("Anurag_CancerDiscov_2022_matrix_gene_level.xlsx")
head(df_test)
dim(df_test)
names(df_test)
df_test %>% select(id:entry_name) %>% distinct() %>% dim()
uniqueN(df_test$id)
uniqueN(df_test$entry_name)
uniqueN(df_test$geneSymbol)

# REMOVE DUPLICATE ROWS PER PROTEIN (KEEP ROW WITH FEWEST MISSING VALUES)
df_test <- df_test %>%
  mutate(na_count = rowSums(is.na(select(., ends_with("Baseline"))))) %>%
  group_by(geneSymbol) %>%
  slice_min(order_by = na_count, with_ties = FALSE) %>% 
  ungroup() %>%
  select(-na_count)

### ------------- NEW SECTION: BACK-TRANSFORMATION TO MATCH MODEL INPUT -------------
# THE EXTERNAL DATA IS PROVIDED AS LOG2-TRANSFORMED TMT RATIOS.
# WE NEED TO RECONSTRUCT "RAW-LIKE" INTENSITIES SO THAT WE CAN CALCULATE FOLD CHANGES FROM RAW DATA.

# STEP A: IDENTIFY INTENSITY COLUMNS
intensity_cols <- grep("Baseline$", names(df_test), value = TRUE)
df_test <- df_test %>% mutate(across(ends_with("Baseline"), as.numeric))

# STEP B: BACK-TRANSFORM FROM LOG2 TO LINEAR SPACE
# (i.e., TMT Ratio in linear scale = 2^(log2 value))
df_linear <- df_test %>%
  mutate(across(all_of(intensity_cols), ~ 2^(.x)))

# STEP C: APPROXIMATE THE POOLED REFERENCE
# (Since the pooled reference is not explicitly given, we use the median
#  of the linear values for each protein as a proxy.)
df_linear <- df_linear %>%
  rowwise() %>%
  mutate(pooled_ref = median(c_across(all_of(intensity_cols)), na.rm = TRUE)) %>%
  ungroup()

# STEP D: RECONSTRUCT RAW-LIKE INTENSITIES
# Multiply each sample’s linear TMT ratio by the approximated pooled reference.
df_raw_like <- df_linear %>%
  mutate(across(all_of(intensity_cols), ~ .x * pooled_ref))

# STEP E: CONVERT THE RAW-LIKE INTENSITIES BACK TO LOG2 SCALE
df_log2_raw <- df_raw_like %>%
  mutate(across(all_of(intensity_cols), ~ log2(.x)))
### ------------- END BACK-TRANSFORMATION SECTION -------------


### ------------- NEW SECTION: IMPUTATION ON LOG2-TRANSFORMED DATA -------------
# IMPUTATION IS BEST PERFORMED ON THE LOG2 DATA, AS THE DISTRIBUTION IS MORE NORMALIZED.
# COMPUTE THE 5TH PERCENTILE ACROSS ALL INTENSITY VALUES (LOG2 SPACE)
low_threshold <- df_log2_raw %>%
  select(all_of(intensity_cols)) %>%
  unlist() %>%
  na.omit() %>%
  quantile(probs = 0.05)

# CREATE A POOL OF LOW-ABUNDANCE VALUES (IN LOG2 SPACE)
low_abundance_pool <- df_log2_raw %>%
  select(all_of(intensity_cols)) %>%
  unlist() %>%
  na.omit() %>%
  .[. <= low_threshold]

# SET SEED FOR REPRODUCIBILITY
set.seed(42)

# IMPUTE MISSING VALUES BY SAMPLING FROM THE LOW-ABUNDANCE POOL
df_imputed <- df_log2_raw %>%
  mutate(across(all_of(intensity_cols), 
                ~ ifelse(is.na(.), sample(low_abundance_pool, size = sum(is.na(.)), replace = TRUE), .)))
### ------------- END IMPUTATION SECTION -------------

# OPTIONAL: COMPARE MISSING VALUES BEFORE & AFTER IMPUTATION
missing_before <- df_log2_raw %>% 
  select(all_of(intensity_cols)) %>% summarise(across(everything(), ~ sum(is.na(.))))
missing_after <- df_imputed %>% 
  select(all_of(intensity_cols)) %>% summarise(across(everything(), ~ sum(is.na(.))))

# (Density/Boxplots/Violin plots can be added as needed for quality check)

# FILTER EXTERNAL DATA TO THE 45 PROTEINS OF INTEREST (USING GENE SYMBOL) AND JOIN WITH 'protein' INFO
df_test_final <- df_imputed %>% 
  filter(geneSymbol %in% protein$gene_name) %>% 
  inner_join(protein, by = c("geneSymbol" = "gene_name"))
dim(df_test_final)
head(df_test_final)

# PREPARE LONG-FORMAT DATA AND COMPUTE FOLD CHANGE (FC)
# CALCULATE FC AS THE DIFFERENCE BETWEEN EACH VALUE AND THE MEAN FOR THAT GENE
dft <- as.data.table(df_test_final %>% 
                       select(-c(id:accession_numbers, entry_name)) %>%
                       pivot_longer(
                         cols = ends_with("Baseline"),  # PATIENT/SAMPLE COLUMNS
                         names_to = "sample",           # NEW COLUMN FOR SAMPLE NAMES
                         values_to = "prot_log"         # LOG2 INTENSITY VALUES
                       ) %>% 
                       mutate(prot_log = as.numeric(prot_log),
                              sample = stringr::str_remove(sample, "_?Baseline")) %>% 
                       select(sample, geneSymbol, gene_id, prot_log))
head(dft)

# COMPUTE PROT_fc FOR EACH GENE BASED ON GENE_ID (MODEL IS TRAINED ON GENE_ID)
dft[, PROT_fc := prot_log - mean(prot_log, na.rm = TRUE), by = gene_id]
summary(dft$PROT_fc)
sum(is.na(dft$PROT_fc)) / nrow(dft)
sum(is.na(dft$prot_log))

# LOAD THE PRE-TRAINED MODEL (TRAINED ON GENE_ID AS FEATURES)
model <- readRDS("~/CEITEC/20241219_Bouchal_TNBC/final_model_TNBC.rds")

# ONLY CONSIDER CASES WHERE ALL OBSERVATIONS ARE AVAILABLE
# CREATE A MATRIX FOR THE MODEL: ROWS = SAMPLES, COLUMNS = GENES (BASED ON GENE_ID)
x <- dcast(dft, gene_id ~ sample, value.var = "PROT_fc")
gene_ids <- x$gene_id
x_matrix <- as.matrix(t(x[, -1, with = FALSE]))
colnames(x_matrix) <- gene_ids

# PREDICT CLUSTERS/GROUPS USING THE MODEL
PROT_clust <- data.frame(sample = names(predict(model, x_matrix)),
                         PROT_clust = predict(model, x_matrix))
table(predict(model, x_matrix))

# ASSIGN GROUPS BACK TO THE ORIGINAL DATASET
fdf <- as.data.table(df_imputed %>% 
                       select(-c(id:accession_numbers, entry_name)) %>%
                       pivot_longer(
                         cols = ends_with("Baseline"),  # PATIENT COLUMNS
                         names_to = "sample",           # NEW COLUMN FOR SAMPLE NAMES
                         values_to = "prot_log"         # LOG2 INTENSITY VALUES
                       ) %>% 
                       mutate(prot_log = as.numeric(prot_log),
                              sample = stringr::str_remove(sample, "_?Baseline")) %>% 
                       select(sample, geneSymbol, prot_log))
rownames(fdf) <- paste(fdf$sample,fdf$geneSymbol, sep = "_")
fdf <- merge(fdf, as.data.frame(PROT_clust), by = "sample", all.x = TRUE)

# OPTIONAL: COMPUTE GROUP-SPECIFIC FOLD CHANGES FOR EACH GENE
lapply(unique(fdf$geneSymbol), function(x) {
  lapply(unique(fdf$PROT_clust), function(y) {
    fdf[geneSymbol == x & PROT_clust == y, PROT_fc_group := 
          mean(fdf[geneSymbol == x & PROT_clust == y, prot_log]) - 
          mean(fdf[geneSymbol == x & PROT_clust != y, prot_log])]
  })
})

head(fdf)
fdf %>% filter(geneSymbol == "A1BG") %>% arrange(PROT_clust) %>% head()

# SAVE FINAL EVALUATION TABLE TO CSV
write.csv(fdf, "~/CEITEC/20241219_Bouchal_TNBC/Res/Anurag_eval_trans.csv", row.names = FALSE)
