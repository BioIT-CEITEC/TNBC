### Separability 
# Load libraries
library(dplyr)
library(data.table)

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


#data<-data %>% dplyr::select(sample,PROT_clust,protein_name:gene_name,prot_raw,prot_log)


data[,PROT_fc := prot_log - mean(prot_log), .(gene_id)]
data <- data[!(sample %in% c("P_010","P_078","P_037","P_079","P_093","P_092"))]

head(data)

#install.packages(c("Rtsne", "umap", "ggplot2", "MASS"))  # Install if not already installed
library(Rtsne)
library(umap)
library(ggplot2)
library(MASS)  # For LDA


# Reshape the data to a matrix where rows are patients and columns are proteins
data_matrix <- dcast(data, sample ~ gene_id, value.var = "PROT_fc")
samples_labs <- data_matrix$sample
data_matrix <- as.matrix(data_matrix[, -1, with = FALSE])  # Remove sample column
rownames(data_matrix) <- samples_labs
# Extract labels (cluster assignments)
labels <- as.factor(data$PROT_clust[match(rownames(data_matrix), data$sample)]) 
labels


### t-SNE-----
#t-SNE is a non-linear dimensionality reduction technique that helps visualize high-dimensional data in 2D.
set.seed(42)  # For reproducibility
tsne_result <- Rtsne(data_matrix, perplexity = 30, theta = 0.5, dims = 2, check_duplicates = FALSE)

# Convert to data frame
tsne_df <- data.frame(X = tsne_result$Y[,1], Y = tsne_result$Y[,2], Cluster = labels)

# Plot t-SNE results
ggplot(tsne_df, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  ggtitle("t-SNE Visualization of Protein Clusters")
ggsave("Pics/tsne_plot.png", width = 10, height = 7, units = "in")


### UMAP-----

set.seed(42)
umap_result <- umap(data_matrix)

# 2. Build the plotting data.frame
umap_df <- data.frame(
  X       = umap_result$layout[, 1],
  Y       = umap_result$layout[, 2],
  Cluster = labels
)

# 3. Define your 7‐cluster color palette
cluster_colors <- c(
  "1" = "#A6CEE3",
  "2" = "#1F78B4",
  "3" = "#B2DF8A",
  "4" = "#33A02C",
  "5" = "#FB9A99",
  "6" = "#E31A1C",
  "7" = "#FDBF6F"
)

# 4. Ensure the Cluster column is a factor with matching levels
umap_df$Cluster <- factor(umap_df$Cluster, levels = names(cluster_colors))

# 5. Plot and save
p_umap <- ggplot(umap_df, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(name = "Cluster", values = cluster_colors) +
  theme_minimal() +
  ggtitle("UMAP Visualization of Protein Clusters")

# Display the plot
print(p_umap)

# Save to file
ggsave("Pics/umap_plot.png", plot = p_umap,
       width = 10, height = 7, units = "in",bg="white") 
### LDA-----
lda_result <- lda(Cluster ~ ., data = as.data.frame(cbind(Cluster = labels, data_matrix)))

# Project data onto first 2 LDA components
lda_values <- predict(lda_result)$x
lda_df <- data.frame(X = lda_values[,1], Y = lda_values[,2], Cluster = labels)

# Plot LDA results
ggplot(lda_df, aes(x = X, y = Y, color = Cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  ggtitle("LDA Visualization of Protein Clusters")
ggsave("Pics/lda_plot.png", width = 10, height = 7, units = "in")
