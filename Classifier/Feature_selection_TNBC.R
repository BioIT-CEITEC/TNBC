 ####### Feature selecition Boucal TNBC---
# ✅ Removes highly correlated proteins (r > 0.9) to prevent redundancy.
# ✅ Filters out low-variance features (bottom 10%) to keep informative features.
# ✅ Uses Boruta to eliminate unimportant features while handling feature interactions.
# ✅ Uses RFE with RF to select the most predictive proteins.
# ✅ Saves final feature matrix for further analysis in Random Forest or other models.



library(dplyr)
library(data.table)
library(data.table)
library(dplyr)
library(randomForest)
library(Boruta)
library(caret)

#### IMPUTS ----------
dir_im<-"TNBC/Classifier"
data_im <- ("full_clust_tab.tsv")


setwd(dir_im)
wd<-getwd()
datum<-Sys.Date()


# Load data   
#data <- fread("all_rna_prot_clust_tab.tsv")
data<-fread(data_im)
head(data)

cor_gene_tpm <- data[,.(corelation = cor(RNA_log,prot_log)),by = .(gene_id,gene_name)]
cor_gene_tpm <- cor_gene_tpm[!is.na(corelation)]

setorder(cor_gene_tpm,-corelation,na.last = T)

gene_number <- 1223

data <- data[gene_id %in% cor_gene_tpm[1:gene_number]$gene_id]

data[,PROT_fc := prot_log - mean(prot_log), .(gene_id)]
data <- data[!(sample %in% c("P_010","P_078","P_037","P_079","P_093","P_092"))]

head(data)


# Reshape the data to a matrix where rows are patients and columns are proteins
prot_matrix  <- dcast(data, sample ~ gene_id, value.var = "PROT_fc")
samples_labs <- prot_matrix $sample
prot_matrix  <- as.matrix(prot_matrix [, -1, with = FALSE])  # Remove sample column
rownames(prot_matrix ) <- samples_labs


# --- 1. Correlation Filtering ---
cor_matrix <- cor(prot_matrix, use = "pairwise.complete.obs")
high_cor_features <- findCorrelation(cor_matrix, cutoff = 0.9, names = TRUE)
prot_matrix <- prot_matrix[, !colnames(prot_matrix) %in% high_cor_features]  # Drop correlated proteins
cat("Removed", length(high_cor_features), "highly correlated proteins.\n")

# --- 2. Variance Thresholding ---
protein_variance <- apply(prot_matrix, 2, var)
low_variance_threshold <- quantile(protein_variance, 0.1)  # 10th percentile
low_var_features <- names(protein_variance[protein_variance < low_variance_threshold])
prot_matrix <- prot_matrix[, !colnames(prot_matrix) %in% low_var_features]
cat("Removed", length(low_var_features), "low variance proteins.\n")

# --- 3. Boruta Feature Selection ---
# Prepare class labels (PROT_clust)
y <- as.factor(data %>% distinct(sample, PROT_clust) %>% pull(PROT_clust))

# Merge labels with filtered features
datas2 <- merge(as.data.frame(data %>% distinct(sample, PROT_clust)), prot_matrix, by.y = "row.names", by.x = "sample") %>% dplyr::select(-sample)

set.seed(123)
boruta_result <- Boruta(PROT_clust ~ ., data = datas2, doTrace = 2, maxRuns = 100)
selected_features_boruta <- getSelectedAttributes(boruta_result, withTentative = TRUE)

cat("Boruta selected", length(selected_features_boruta), "features.\n")
prot_matrix <- prot_matrix[, colnames(prot_matrix) %in% selected_features_boruta]

# --- 4. Recursive Feature Elimination (RFE) ---
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(123)
rfe_result <- rfe(x = prot_matrix, y = y, sizes = seq(5, ncol(prot_matrix), by = 5), rfeControl = control)
selected_features_rfe <- predictors(rfe_result)

cat("RFE selected", length(selected_features_rfe), "features.\n")

# Final feature matrix
final_prot_matrix <- prot_matrix[, colnames(prot_matrix) %in% selected_features_rfe]

# Save results
saveRDS(final_prot_matrix, paste0("Filtered_Features", ".rds"))

cat("Final feature selection completed. Features saved for further modeling.\n")


