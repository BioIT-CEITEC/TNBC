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

setwd("~/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC")
wd<-getwd()
datum<-Sys.Date()


# Load data   
#data <- fread("all_rna_prot_clust_tab.tsv")
data<-fread("~/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Katka_files/results_2025/full_clust_tab.tsv")
head(data)

cor_gene_tpm <- data[,.(corelation = cor(RNA_log,prot_log)),by = .(gene_id,gene_name)]
cor_gene_tpm <- cor_gene_tpm[!is.na(corelation)]

setorder(cor_gene_tpm,-corelation,na.last = T)

gene_number <- 1223

data <- data[gene_id %in% cor_gene_tpm[1:gene_number]$gene_id]


#data<-data %>% dplyr::select(sample,PROT_clust,protein_name:gene_name,prot_raw,prot_log)

# Proteinové fc su absolutne mimo - vypočítať znova !

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
saveRDS(final_prot_matrix, paste0("Res/Filtered_Features_", datum, ".rds"))

cat("Final feature selection completed. Features saved for further modeling.\n")


###### Alternative with UMAP for dim reduction
# Prepare data for UMAP
library(umap)
x  <- dcast(data, sample ~ gene_id, value.var = "PROT_fc")
samples_labs <- x $sample
x  <- as.matrix(x [, -1, with = FALSE])  # Remove sample column
rownames(x ) <- samples_labs

# Normalize data
x <- scale(x)  

# Run UMAP
umap_model <- umap(x, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_data <- as.data.frame(umap_model$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")

# Add classification labels
umap_data$sample <- rownames(x)  # Assign correct sample names

# Merge with PROT_clust labels
umap_data <- merge(umap_data, unique(data[, .(sample, PROT_clust)]), by = "sample")
umap_data$PROT_clust <- as.factor(umap_data$PROT_clust)

# Plot UMAP
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = PROT_clust)) +
  geom_point(size = 3) +
  ggtitle("UMAP Projection of Protein Intensities") +
  theme_minimal()

# Train-test split
set.seed(123)
train_index <- createDataPartition(umap_data$PROT_clust, p = 0.8, list = FALSE)
train_data <- umap_data[train_index, ]
test_data <- umap_data[-train_index, ]

# Define hyperparameter grid for tuning
rf_grid <- expand.grid(mtry = c(2, 5, 10))  # RF: Number of features per split
# Define hyperparameter grid for tuning
xgb_grid <- expand.grid(
  nrounds = seq(100, 500, 100),   # Number of boosting rounds
  max_depth = c(3, 6, 9),         # Tree depth
  eta = c(0.01, 0.1, 0.3),        # Learning rate
  gamma = c(0, 0.1, 1),           # Minimum loss reduction
  colsample_bytree = c(0.5, 0.8, 1),  # Feature fraction per tree
  min_child_weight = c(1, 5, 10), # Minimum sum of instance weights in a node
  subsample = c(0.5, 0.75, 1)     # Row subsampling
)
# Train RF model with hyperparameter tuning
set.seed(123)
rf_tuned <- train(PROT_clust ~ ., data = train_data, method = "rf",
                  trControl = trainControl(method = "cv", number = 5),
                  tuneGrid = rf_grid)

# Train XGBoost model with hyperparameter tuning
set.seed(123)
# Train tuned XGBoost model
xgb_tuned <- train(
  PROT_clust ~ ., 
  data = train_data, 
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 5), # 5-fold cross-validation
  tuneGrid = xgb_grid
)
# Evaluate models
rf_pred <- predict(rf_tuned, test_data)
xgb_pred <- predict(xgb_tuned, test_data)

rf_acc <- mean(rf_pred == test_data$PROT_clust)
xgb_acc <- mean(xgb_pred == test_data$PROT_clust)

print(paste("Random Forest Accuracy:", rf_acc))
print(paste("XGBoost Accuracy:", xgb_acc))

# Compare feature importance
rf_importance <- varImp(rf_tuned)
xgb_importance <- varImp(xgb_tuned)

plot(rf_importance, main = "RF Feature Importance")
plot(xgb_importance, main = "XGB Feature Importance")