# ==================== Code 1: Random Search with Best Iteration Logging ====================
# Load required libraries
library(data.table)
library(dplyr)
library(randomForest)
library(caret)
library(openxlsx)
library(rstatix)
library(data.table)
library(dplyr)
library(randomForest)
library(caret)
library(rBayesianOptimization)

### IMPUTS ----------
dir_im<-"TNBC/Classifier"
df_im <- "Filtered_Features.rds"
data_im <- "full_clust_tab.tsv"


setwd(dir_im)
wd<-getwd()
datum<-Sys.Date()


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
# ==================== Bayesian Optimization with Improved Handling ====================


# Read filtered features and filter data accordingly
df <- as.data.frame(readRDS(df_im))
data <- data %>% filter(gene_id %in% colnames(df))

# Reshape data: rows = samples, columns = proteins
x <- dcast(data, gene_id ~ sample, value.var = "PROT_fc")
gene_id <- x$gene_id
x <- as.matrix(t(x[, -1, with = FALSE]))
colnames(x) <- gene_id
y <- as.data.frame(data %>% distinct(sample, PROT_clust))
rownames(y) <- y$sample
y <- as.data.frame(y %>% dplyr::select(-sample))
datas2 <- merge(y, x, by.y = "row.names", by.x = "row.names") %>% dplyr::select(-Row.names)
datas2$PROT_clust <- as.factor(datas2$PROT_clust)

# Total number of proteins (features)
total <- dim(x)[2]

# Define search bounds for hyperparameters
lower_mtry <- floor(sqrt(total))
upper_mtry <- floor(total / 2)
lower_nodesize <- 1
upper_nodesize <- 10

# Define the objective function for Bayesian Optimization
rf_cv_bayes <- function(mtry, nodesize) {
  mtry <- as.integer(round(mtry))
  nodesize <- as.integer(round(nodesize))
  n_iter <- 10  # number of iterations to average over
  acc_vec <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    # For each unique class, sample one index
    unique_classes <- unique(datas2$PROT_clust)
    inds <- sapply(unique_classes, function(cl) {
      idx <- which(datas2[, 1] == cl)
      if (length(idx) > 0) sample(idx, 1) else NA
    })
    # If any class is missing in the hold-out, assign 0 and continue
    if (any(is.na(inds)) || length(inds) < length(unique_classes)) {
      acc_vec[i] <- 0
      next
    }
    
    tryCatch({
      # Compute class weights for training set
      train_labels <- factor(datas2[-inds, 1])
      tbl <- table(train_labels)
      class_weights <- as.vector(max(tbl) / tbl)
      names(class_weights) <- names(tbl)
      
      model <- randomForest::randomForest(
        x = datas2[-inds, -1],
        y = factor(datas2[-inds, 1]),
        data = datas2[-inds, ],
        ntree = 500,  # fewer trees for speed during optimization
        mtry = mtry,
        nodesize = nodesize,
        classwt = class_weights,
        importance = FALSE,
        replace = TRUE
      )
      
      preds <- predict(model, datas2[inds, -1])
      actual <- factor(datas2[inds, 1])
      pred <- factor(preds, levels = levels(actual))
      cm <- caret::confusionMatrix(pred, actual)
      if ("Balanced Accuracy" %in% colnames(cm$byClass)) {
        acc_vec[i] <- mean(cm$byClass[,"Balanced Accuracy"], na.rm = TRUE)
      } else {
        acc_vec[i] <- mean(cm$byClass[,"Sensitivity"], na.rm = TRUE)
      }
    }, error = function(e) {
      acc_vec[i] <- 0
    })
  }
  
  MeanBA <- mean(acc_vec, na.rm = TRUE)
  if (is.nan(MeanBA)) MeanBA <- 0
  return(list(Score = MeanBA))
}

set.seed(100)
opt_result <- BayesianOptimization(
  FUN = rf_cv_bayes,
  bounds = list(mtry = c(lower_mtry, upper_mtry),
                nodesize = c(lower_nodesize, upper_nodesize)),
  init_points = 10,
  n_iter = 20,
  acq = "ucb",
  kappa = 2.576,
  eps = 0.0,
  verbose = TRUE
)

print(opt_result)
saveRDS(opt_result,"opt_result.rds")
# ==================== End Code 2 ====================
