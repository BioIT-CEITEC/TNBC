# ==================== Code 1: Random Search with Best Iteration Logging ====================
# Load required libraries
library(data.table)
library(dplyr)
library(randomForest)
library(caret)
library(openxlsx)
library(rstatix)

# (Assuming 'data' is already loaded in your workspace)

# Read filtered features and filter data accordingly
df <- readRDS("~/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Res/Filtered_Features_2025-02-27.rds")
data <- data %>% filter(gene_id %in% colnames(df))

# Create Excel workbook and formatting style
datum <- Sys.Date()
wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
wb <- createWorkbook(Sys.setenv("R_ZIPCMD" = "c:/RBuildTools/3.5/bin/zip.exe"))
hs1 <- createStyle(halign = "CENTER", textDecoration = "Bold", border = "Bottom", fontColour = "black")
addWorksheet(wb, "TNBC")

# Reshape the data: rows = samples, columns = proteins
x <- dcast(data, gene_id ~ sample, value.var = "PROT_fc")
gene_id <- x$gene_id
x <- as.matrix(t(x[, -1, with = FALSE]))
colnames(x) <- gene_id
y <- as.data.frame(data %>% distinct(sample, PROT_clust))
rownames(y) <- y$sample
y <- as.data.frame(y %>% dplyr::select(-sample))

# Number of groups and total proteins
numg <- uniqueN(y$PROT_clust)
total <- dim(x)[2]

# Merge labels with data so that datas2 has sample labels in the first column
datas2 <- merge(y, x, by.y = "row.names", by.x = "row.names") %>% dplyr::select(-Row.names)
datas2$PROT_clust <- as.factor(datas2$PROT_clust)

#### Set number of trees and iterations ----
nr <- 100   # number of iterations (random search runs)
nt <- 1000  # number of trees in RF

# Define candidate hyperparameter values
mtry_candidates <- unique(c(floor(sqrt(total)), floor(total/3), floor(total/2)))
nodesize_candidates <- c(1, 5, 10)

# Initialize matrices to log feature importance and performance
imp <- matrix(0, nrow = dim(x)[2], ncol = nr)
acc <- matrix(0, nrow = dim(x)[2], ncol = nr)
for (i in paste0("g", 1:numg)) {
  assign(i, matrix(0, nrow = dim(x)[2], ncol = nr))
}
ta <- 0  # to accumulate confusion matrices

# Data frame to log hyperparameters and performance per iteration
hp_results <- data.frame(iteration = 1:nr, mtry = NA, nodesize = NA, balanced_acc = NA)
performance <- numeric(nr)

set.seed(100)
for (j in 1:nr) {
  # Stratified hold-out: select one sample per group
  inds <- matrix(0, nrow = numg, ncol = 1)
  if (min(table(y$PROT_clust)) == 1) {
    for (i in as.numeric(names(table(y$PROT_clust)[table(y$PROT_clust) != 1]))) {
      ind <- sample(which(datas2[,1] == unique(datas2[,1])[i]), 1)
      inds[i] <- ind  
    }
  } else {
    for (i in 1:numg) {
      ind <- sample(which(datas2[,1] == unique(datas2[,1])[i]), 1)
      inds[i] <- ind
    }
  }
  
  # Randomly choose hyperparameters
  chosen_mtry <- sample(mtry_candidates, 1)
  chosen_nodesize <- sample(nodesize_candidates, 1)
  hp_results$mtry[j] <- chosen_mtry
  hp_results$nodesize[j] <- chosen_nodesize
  
  # Compute class weights from training set (datas2[-inds,])
  train_labels <- factor(datas2[-inds, 1])
  tbl <- table(train_labels)
  class_weights <- as.vector(max(tbl) / tbl)
  names(class_weights) <- names(tbl)
  
  # Train Random Forest with chosen hyperparameters and class weights
  borut <- randomForest::randomForest(
    x = datas2[-inds, -1],
    y = factor(datas2[-inds, 1]),
    data = datas2[-inds, ],
    ntree = nt,
    mtry = chosen_mtry,
    nodesize = chosen_nodesize,
    classwt = class_weights,
    importance = TRUE,
    replace = TRUE
  )
  
  # Log importance measures
  imp[, j] <- importance(borut)[, "MeanDecreaseGini"]
  acc[, j] <- importance(borut)[, "MeanDecreaseAccuracy"]
  
  # For each group, store group-specific importance
  for (i in 1:numg) {
    matrix_name <- paste0("g", i)
    eval(parse(text = paste0(matrix_name, "[,j] <- importance(borut)[,", i, "]")))
  }
  
  # Predict on held-out samples and compute balanced accuracy
  ptest <- predict(borut, datas2[inds, -1])
  actual <- factor(datas2[inds, 1])
  pred <- factor(ptest, levels = levels(actual))
  conf_mat <- caret::confusionMatrix(pred, actual)
  # If "Balanced Accuracy" is available, use it; otherwise, average sensitivity
  if ("Balanced Accuracy" %in% colnames(conf_mat$byClass)) {
    bal_acc <- mean(conf_mat$byClass[,"Balanced Accuracy"], na.rm = TRUE)
  } else {
    bal_acc <- mean(conf_mat$byClass[,"Sensitivity"], na.rm = TRUE)
  }
  performance[j] <- bal_acc
  hp_results$balanced_acc[j] <- bal_acc
  
  # Accumulate confusion matrix
  ta <- table(datas2[inds, 1], ptest) + ta
}

# Average importance over iterations
means <- cbind(rowMeans(imp))
for (i in 1:numg) {
  assign(paste0("g", i), cbind(rowMeans(get(paste0("g", i)))))
}
rownames(means) <- rownames(imp) <- colnames(datas2[, 2:ncol(datas2)])
el <- as.data.frame(t(imp)) %>% rstatix::get_summary_stats(type = "mean_se")
dt <- data.frame(el[c(1, 3:4)])
colnames(dt) <- c("MDI", "Means", "SE")
imps <- dt[order(-dt$Means), ]

# Identify best iteration (highest balanced accuracy)
best_iter <- which.max(performance)
best_hp <- hp_results[best_iter, ]
print(paste("Best iteration:", best_iter))
print(best_hp)
print(summary(hp_results))

# Write results to Excel
for (i in 1:numg) {
  g <- get(paste0("g", i))
  rownames(g) <- colnames(datas2[, 2:ncol(datas2)])
  g <- arrange(as.data.frame(g), desc(V1))
  colnames(g) <- levels(datas2$PROT_clust)[i]
  addWorksheet(wb, levels(datas2$PROT_clust)[i])
  writeData(wb, levels(datas2$PROT_clust)[i], head(g, 500),
            startRow = 1, startCol = 1, colNames = TRUE, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, levels(datas2$PROT_clust)[i], "Pacienti v skupine",
            startRow = 1, startCol = 4, colNames = TRUE, rowNames = TRUE, headerStyle = hs1)
  writeData(wb, levels(datas2$PROT_clust)[i], names(y[y == i]),
            startRow = 2, startCol = 4, colNames = TRUE, rowNames = TRUE, headerStyle = hs1)
  assign(paste0("g", i), g)
}

writeData(wb, "TNBC", head(imps, 500), startRow = 1, startCol = 1,
          colNames = TRUE, rowNames = FALSE, headerStyle = hs1)
writeData(wb, "TNBC", "Confusion matrix", startRow = 1, startCol = 5,
          colNames = TRUE, rowNames = TRUE, headerStyle = hs1)
writeData(wb, "TNBC", ta, startRow = 2, startCol = 5, colNames = TRUE, rowNames = TRUE)

print(head(as.data.frame(head(g3))))
saveWorkbook(wb, paste0("Res/Random_forests_PreselectedFeatures", datum, ".xlsx"), overwrite = TRUE)
# ==================== End Code 1 ====================
# ==================== Code 2: Bayesian Optimization with Improved Handling ====================
# Load required libraries
library(data.table)
library(dplyr)
library(randomForest)
library(caret)
library(rBayesianOptimization)

# Read filtered features and filter data accordingly
df <- readRDS("~/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Res/Filtered_Features_2025-02-27.rds")
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
