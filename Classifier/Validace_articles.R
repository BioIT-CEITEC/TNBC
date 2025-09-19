# ARTICLES eval - 
# 1) skontrolovať či je 45 zvolených poroteinov dostupných v dátach 
# vyrobiť skupiny pre random forest  
# spočítať PROT_fc skupin vs všetci pre všetky dostupné proteinové expresie

# Load libraries
library(dplyr)
library(data.table)
library(randomForest)

setwd("~/CEITEC/20241219_Bouchal_TNBC") #/OneDrive/Documents/
wd<-getwd()
datum<-Sys.Date()

df <- readRDS("~/CEITEC/20241219_Bouchal_TNBC/Res/Filtered_Features_2025-02-27.rds")
data<-fread("~/CEITEC/20241219_Bouchal_TNBC/Katka_files/results_2025/full_clust_tab.tsv")#/OneDrive/Documents/
head(data)
head(df)

proteiny<-colnames(df)[grep("ENSG",colnames(df))]
length(proteiny)
proteiny[proteiny %in% data$gene_id]
length(proteiny[proteiny %in% data$gene_id])

protein<-data %>% filter(gene_id %in% proteiny) %>% select(gene_id,gene_name) %>% distinct()

# read in all articles 
library(readxl)
library(tools)  # for file_path_sans_ext()

# Define the folder path (use forward slashes or double backslashes)
folder_path <- "C:/Users/FNOadmin/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Articles"

# List all .xlsx files in the folder (with full paths)
files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Read each file into a list; names are the file names without the extension
data_list <- lapply(files, read_excel)
names(data_list) <-articles<- sub("_proteinmatrix","",file_path_sans_ext(basename(files)))

for (i in articles){
  print(i)
  #print(head(data_list[[i]]))
  print(names(data_list[[i]]))
}

library(readxl)
library(tools)


# Assume 'data_list' is your named list of article data frames
# For example, it could have been created like this:
# folder_path <- "C:/Users/FNOadmin/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Articles"
# files <- list.files(path = folder_path, pattern = "\\.xlsx$", full.names = TRUE)
# data_list <- lapply(files, read_excel)
# names(data_list) <- file_path_sans_ext(basename(files))

# Function to determine the best candidate column for protein identifiers in a data frame.
get_best_identifier_column <- function(df, protein) {
  best_col <- NULL
  best_overlap <- 0
  
  # Lower-case the protein identifiers from your list for comparison
  prot_ids <- tolower(protein$gene_id)
  prot_names <- tolower(protein$gene_name)
  
  # Iterate over each column in the article data frame
  for (col in names(df)) {
    # Convert the column values to character and lower-case
    vals <- tolower(as.character(df[[col]]))
    
    # Compute overlap with both gene_id and gene_name from protein
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

# Now, iterate through each article to check if it contains all 45 proteins.
articles_with_all_proteins <- c()

for (article_name in names(data_list)) {
  df <- data_list[[article_name]]
  
  candidate <- get_best_identifier_column(df, protein)
  best_col <- candidate$column
  overlap_count <- candidate$overlap
  
  cat("Article:", article_name, "\n")
  if (!is.null(best_col)) {
    cat("  Best identifier column:", best_col, "\n")
    cat("  Overlap count:", overlap_count, "\n")
    
    # If the overlap equals 45 (all proteins found), mark this article.
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

# jjhu mame víťaza 
#Anurag_CancerDiscov_2022
df_test<-data_list[["Anurag_CancerDiscov_2022"]]
df_test<-read_excel("Anurag_CancerDiscov_2022_matrix_gene_level.xlsx")

head(df_test)
dim(df_test)
names(df_test)
df_test %>% select(id:entry_name)%>% distinct%>%dim()
uniqueN(df_test$id)
uniqueN(df_test$entry_name)
uniqueN(df_test$geneSymbol)


library(dplyr)
library(tidyr)

#nemože byť viac riadkov pre jeden protein nameraných - vymažem ten kde je viac riadkov missing
df_test %>% group_by(geneSymbol)%>% summarise(n=n()) %>% arrange(desc(n))

df_test <- df_test %>%
  mutate(na_count = rowSums(is.na(select(., ends_with("Baseline"))))) %>%
  group_by(geneSymbol) %>%
  slice_min(order_by = na_count, with_ties = FALSE) %>% 
  ungroup() %>%
  select(-na_count)

#### IMputation based on global imputing
# Step 1: Identify intensity columns (assuming they end with "Baseline")
intensity_cols <- grep("Baseline$", names(df_test), value = TRUE)
str(df_test)

df_test <- df_test %>%
  mutate(across(ends_with("Baseline"), as.numeric))
# Step 2: Identify low-abundance values (5th percentile of log2 values)
low_threshold <- df_test %>%
  select(all_of(intensity_cols)) %>%
  unlist() %>%
  na.omit() %>%
  quantile(probs = 0.05)  # 5th percentile in log2 space

# Extract all values below this threshold (creating a low-abundance pool)
low_abundance_pool <- df_test %>%
  select(all_of(intensity_cols)) %>%
  unlist() %>%
  na.omit() %>%
  .[. <= low_threshold]  # Keep values in log2 space

# Step 3: Impute missing values by randomly sampling from low-abundance pool
set.seed(42)  # For reproducibility

df_imputed <- df_test %>%
  mutate(across(all_of(intensity_cols), 
                ~ ifelse(is.na(.), sample(low_abundance_pool, size = sum(is.na(.)), replace = TRUE), .)))

# Check the first few rows
head(df_imputed)
sum(is.na(df_imputed))
summary(df_imputed)
summary(df_test)

# IMPUTATION COMPARISON ------

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots

# Step 1: Count missing values before and after
missing_before <- df_test %>%
  select(ends_with("Baseline")) %>%
  summarise(across(everything(), ~ sum(is.na(.)), .names = "missing_{.col}")) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "missing_count")

missing_after <- df_imputed %>%
  select(ends_with("Baseline")) %>%
  summarise(across(everything(), ~ sum(is.na(.)), .names = "missing_{.col}")) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "missing_count")

# Step 2: Compare the density distribution before & after imputation
df_long_before <- df_test %>%
  pivot_longer(cols = ends_with("Baseline"), names_to = "sample", values_to = "prot_log") %>%
  mutate(status = ifelse(is.na(prot_log), "Missing", "Observed"))

df_long_after <- df_imputed %>%
  pivot_longer(cols = ends_with("Baseline"), names_to = "sample", values_to = "prot_log") %>%
  mutate(status = "Imputed")

df_density <- bind_rows(df_long_before, df_long_after)

p_density <- ggplot(df_density %>% filter(geneSymbol %in% protein$gene_name), aes(x = prot_log, fill = status)) +
  geom_density(alpha = 0.5) +
  facet_grid(~status)+
  labs(title = "Density of Protein Intensities Before vs After Imputation",
       x = "Log2 Intensity", y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("Missing" = "red", "Observed" = "blue", "Imputed" = "green"))
p_density
# Step 3: Boxplot of protein intensities (before vs after) for the 45 proteins
df_45_before <- df_test %>%
  filter(geneSymbol %in% protein$gene_name) %>%
  pivot_longer(cols = ends_with("Baseline"), names_to = "sample", values_to = "prot_log") %>%
  mutate(status = as.character(ifelse(is.na(prot_log), "Missing", "Observed")))

df_45_after <- df_imputed %>%
  filter(geneSymbol %in% protein$gene_name) %>%
  pivot_longer(cols = ends_with("Baseline"), names_to = "sample", values_to = "prot_log") %>%
  mutate(status = "Imputed")

df_45_combined <- bind_rows(df_45_before, df_45_after)

p_boxplot <- ggplot(df_45_combined, aes(x = geneSymbol, y = prot_log, fill = status)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  labs(title = "Boxplot of 45 Proteins Before vs After Imputation",
       x = "Gene ID", y = "Log2 Intensity") +
  theme_minimal() +
  scale_fill_manual(values = c("Missing" = "red", "Observed" = "blue", "Imputed" = "green"))
p_boxplot

# Step 4: Violin plot to compare imputed vs. non-imputed values
p_violin <- ggplot(df_45_combined, aes(x = status, y = prot_log, fill = status)) +
  geom_violin(alpha = 0.5) +
  labs(title = "Violin Plot of 45 Proteins Before vs After Imputation",
       x = "Status", y = "Log2 Intensity") +
  theme_minimal() +
  scale_fill_manual(values = c("Missing" = "red", "Observed" = "blue", "Imputed" = "green"))
p_violin

  df_test<-df_imputed %>% filter(geneSymbol%in% protein$gene_name) %>% inner_join(protein, by = c("geneSymbol" = "gene_name")) 
dim(df_test )
head(df_test)

library(tidyr)

dft <- as.data.table(df_test %>% select(-c(id:accession_numbers,entry_name))%>%
  pivot_longer(
    cols = ends_with("Baseline"),  # patient columns
    names_to = "sample",           # new column for sample names
    values_to = "prot_log"         # log2 intensity values
  )%>% 
  mutate(prot_log= as.numeric(prot_log),
         sample = stringr::str_remove(sample, "_?Baseline"))%>% 
  select(sample,geneSymbol,gene_id,prot_log) )

head(dft)
head(data)
str(dft)
uniqueN(dft$sample)
dft%>% filter(is.na(as.numeric(prot_log)))%>% group_by(geneSymbol)%>% summarise(n=n())

uniqueN(dft$sample)

library(data.table)
dft[,PROT_fc := prot_log - mean(prot_log, na.rm = TRUE), .(gene_id)]
summary(dft$PROT_fc)
sum(is.na(dft$PROT_fc)/nrow(dft))
sum(is.na(dft$prot_log))

model <- readRDS("~/CEITEC/20241219_Bouchal_TNBC/final_model_TNBC.rds")

#Only cinsider when all observation available
x<-dcast(dft, gene_id ~ sample, value.var = "PROT_fc")
gene_id <- x$gene_id
x <- (as.matrix(t(x[, -1, with = FALSE])))
colnames(x) <- gene_id

PROT_clust<-data.frame(sample=names(predict(model,x)),PROT_clust=predict(model,x))
table(predict(model,x))


# priradenie skupín k pôvodnéu datasetu
fdf<- as.data.table(df_imputed %>% select(-c(id:accession_numbers,entry_name))%>%
                       pivot_longer(
                         cols = ends_with("Baseline"),  # patient columns
                         names_to = "sample",           # new column for sample names
                         values_to = "prot_log"         # log2 intensity values
                       )%>% 
                       mutate(prot_log= as.numeric(prot_log),
                              sample = stringr::str_remove(sample, "_?Baseline"))%>% 
                       select(sample,geneSymbol,prot_log) )
rownames(fdf)<-fdf$sample
dim(fdf)
head(fdf)

fdf<-data.table::merge.data.table(fdf,as.data.frame(PROT_clust))

lapply(unique(fdf$geneSymbol), function(x) {
  dt_list <- lapply(unique(fdf$PROT_clust), function(y) {
    # Compute log2FC using already calculated log values
    fdf[geneSymbol == x & PROT_clust == y, PROT_fc := 
          mean(fdf[geneSymbol == x & PROT_clust == y, prot_log]) - 
          mean(fdf[geneSymbol == x & PROT_clust != y, prot_log])]
  })
})


head(fdf)
fdf %>% filter(geneSymbol=="A1BG") %>%arrange(PROT_clust)%>% head()

# ulož csv
write.csv(fdf,"~/CEITEC/20241219_Bouchal_TNBC/Res/Anurag_eval.csv", row.names = FALSE)
