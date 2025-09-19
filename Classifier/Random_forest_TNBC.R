### Klasifikátor - based on PROTEINY ?
 rm(list=ls())
# Load libraries
library(dplyr)
library(data.table)
 library(randomForest)
 ### Autosaving in xlsx
 library(openxlsx)

setwd("~/CEITEC/20241219_Bouchal_TNBC") #/OneDrive/Documents/
wd<-getwd()
datum<-Sys.Date()


# Load data   
#data <- fread("all_rna_prot_clust_tab.tsv")
data<-fread("~/CEITEC/20241219_Bouchal_TNBC/Katka_files/results_2025/full_clust_tab.tsv")#/OneDrive/Documents/
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


# Load necessary libraries
library(data.table)


#Initialize result list
res <- list()
left_out<-unique(data$sample)[1]
# Perform leave-one-out  in hierarchical clustering
for (left_out in unique(data$sample)) {
  cat("Leaving out:", left_out, "\n")

  # Training set: all samples except the one left out
  training_data <- data[sample != left_out]

  # Test set: the left-out sample
  test_data <- data[sample == left_out]

  #training_data <- data
  
  # Prepare matrix for clustering
  train_matrix <- dcast(training_data, gene_id ~ sample, value.var = "PROT_fc")
  rownames(train_matrix) <- train_matrix$gene_id
  train_matrix <- as.matrix(train_matrix[, -1, with = FALSE])
  
  # Perform hierarchical clustering on the training set
  hc <- hclust(dist(t(train_matrix)), method = "ward.D")
  
  # Assign cluster labels to the training set
  cluster_labels <- cutree(hc, k = 7)  # Adjust `k` as needed
  

  # Chceck the results agains Previous clust in Katka file - if they cluster 
   kontorla <-merge(as.data.frame(data%>%dplyr::select(sample,PROT_clust)%>%distinct(sample,.keep_all = T)),
        as.data.frame(cluster_labels),by.x="sample",by.y=0,all.x = TRUE) %>% 
    distinct(sample,.keep_all = TRUE) %>%
   # {table(.$cluster_labels, .$PROT_clust)} # alternatíva k tabulke
    count(cluster_labels, PROT_clust)  
   #print(kontorla)  
  if (nrow(na.omit(kontorla)) == 7) {
    print("OK - clustrovanie konzistetné!")
  } else {
    print("Nekonzistetné clustrovanie")
    print(kontorla)  
  }
}
table(cluster_labels)

# RF ------
data_save<-data

{
  wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
  wb <- createWorkbook(Sys.setenv("R_ZIPCMD" = "c:/RBuildTools/3.5/bin/zip.exe"))
  
  #formaring style
  hs1 <- createStyle( halign = "CENTER", textDecoration = "Bold",border = "Bottom", fontColour = "black")
  #excel result worksheet
  addWorksheet(wb, "TNBC")
  
  
  #Redesign the data in a way that there will be individual lines on rows
  #and proteins  in columns. 
  
  
  #Only cinsider when all observation available
  x<-dcast(data, gene_id ~ sample, value.var = "PROT_fc")
  gene_id <- x$gene_id
  x <- (as.matrix(t(x[, -1, with = FALSE])))
  colnames(x) <- gene_id
  y<-as.data.frame(data %>% distinct(sample,PROT_clust))
   rownames(y)<-y$sample
   sum(rownames(x) %in% y$sample)
   y<-as.data.frame(y %>% dplyr::select(-sample))
 #check if the columns are the same
  
 
  
   
  #number of group
  numg<-uniqueN(y$PROT_clust)
  
  #total number of proteins
  total<-dim(x)[2]
  
  datas2<-merge(y,x,by.y = "row.names",by.x = "row.names")%>% dplyr::select(-Row.names)
  y<-as.factor(t(y$PROT_clust))
  datas2[1:6,1:6]
  datas2$PROT_clust<-as.factor(datas2$PROT_clust)
  #### Number of tree and repetition ----
  nr = 100
  nt = 1000
  
  
  #implementation of random forest classifier with nr repeats.
  imp<-matrix(0,(dim(x)[2]),nr)
  acc<-matrix(0,(dim(x)[2]),nr)
  #for each group separately
  for (i in paste0("g",1:numg)){
    assign(i,matrix(0,(dim(x)[2]),nr)) 
  }
  
  ta<-0
  set.seed(100)
  for (j in 1:nr){
    inds<-matrix(0,numg,1)
    #stratified hold-out cross validation (1 obs. per group)
    if (min(table(y))==1){ # v prípade že je v skupine iba jeden pacient nebudeme mať testovaciu sadu, iba trenovaciu - TA bude zla
      for (i in as.numeric(names(table(y)[table(y)!=1]))){ 
        
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) # v prípade že ide len o jednu sample a to je inteager  x >= 1, sampling via sample takes place from 1:x
        inds[i]<-ind  
      }
    }else {
      for (i in 1:numg){ 
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) 
        inds[i]<-ind
      }
    }
    
    
    borut<-randomForest::randomForest(datas2[-inds,-1],factor(datas2[-inds,1]),data=datas2[-inds,],ntree=nt,importance = T,replace = T)
    head(importance(borut))
    imp[,j]<-importance(borut)[,"MeanDecreaseGini"]
    acc[,j]<-importance(borut)[,"MeanDecreaseAccuracy"]
    
    #importance for each group
    # Assign importance for each group
    for (i in 1:numg) {
      matrix_name <- paste0("g", i)
      eval(parse(text = paste0(matrix_name, "[,j] <- importance(borut)[,", i, "]")))
    }
    
    #alternative - does not work
    
    #purrr::map(1:numg, ~ assign(paste0("g", .x, "[,j]"), importance(borut)[,.x]))
    
    
    
    ptest<-predict(borut,datas2[inds,-1])
    ta<-table(datas2[inds,1],ptest)+ta
  }
  means<-cbind(rowMeans(imp))
  
  for (i in 1:numg) {
    assign(paste0("g", i), cbind(rowMeans(get(paste0("g", i)))))
  }
  
  # g1<-cbind(rowMeans(g1))
  
  
  rownames(means)<-rownames(imp)<-colnames(datas2[,2:(dim(datas2)[2])])
  #calculate rows variances from MDI (Mean Decrease in Impurity)
  el<-as.data.frame(t(imp)) %>%  rstatix::get_summary_stats(type="mean_se")
  
  dt<-data.frame(el[c(1,3:4)])
  colnames(dt)<-c("MDI","Means","SE")
  imps<-dt[order(-dt$Means),]
  
  
  #Write in Excel ----
  for (i in 1:numg){
    g<-get(paste0("g",i))
    rownames(g)<-colnames(datas2[,2:(dim(datas2)[2])])
    #calculate rows variances from MDI (Mean Decrease in Impurity)
    g<-arrange(as.data.frame(g),desc(V1))
    colnames(g)<-levels(datas2$PROT_clust)[i]
    addWorksheet(wb,levels(datas2$PROT_clust)[i])
    writeData(wb,levels(datas2$PROT_clust)[i],head(g,500), startRow = 1, startCol = 1, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],"Pacienti v skupine", startRow = 1, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],names(y[y==i]), startRow = 2, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    
    assign(paste0("g",i),g)
  }

  # Save the results
  writeData(wb,"TNBC",head(imps,500), startRow = 1, startCol = 1, colNames = T, rowNames = F, headerStyle = hs1)  
  writeData(wb,"TNBC","Confusion matrix", startRow = 1, startCol = 5, colNames = T, rowNames = T, headerStyle = hs1) 
  writeData(wb,"TNBC",ta, startRow = 2, startCol = 5, colNames = T, rowNames = T) 

  print(head(as.data.frame(head(g3))))
  # Save the workbook
  saveWorkbook(wb,paste0("Res/Random_forests_",datum,".xlsx"),overwrite = T)
  }
########################################################################
# Selected features based on MDi.R---------
imps%>% head(100)

#select features with MDI >= 2
features<-imps$MDI[imps$Means>=0.2]
length(features)
#select first 50
# features<-imps$MDI[1:50]

data %>% filter(gene_id %in% features) %>% distinct(gene_id) %>% nrow()
data<-data_save %>% filter(gene_id %in% colnames(df))

{
  wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
  wb <- createWorkbook(Sys.setenv("R_ZIPCMD" = "c:/RBuildTools/3.5/bin/zip.exe"))
  
  #formaring style
  hs1 <- createStyle( halign = "CENTER", textDecoration = "Bold",border = "Bottom", fontColour = "black")
  #excel result worksheet
  addWorksheet(wb, "TNBC")
  
  
  #Redesign the data in a way that there will be individual lines on rows
  #and proteins  in columns. 
  
  
  #Only cinsider when all observation available
  x<-dcast(data, gene_id ~ sample, value.var = "PROT_fc")
  gene_id <- x$gene_id
  x <- (as.matrix(t(x[, -1, with = FALSE])))
  colnames(x) <- gene_id
  y<-as.data.frame(data %>% distinct(sample,PROT_clust))
  rownames(y)<-y$sample
  sum(rownames(x) %in% y$sample)
  y<-as.data.frame(y %>% dplyr::select(-sample))
  #check if the columns are the same
  
  
  
  
  #number of group
  numg<-uniqueN(y$PROT_clust)
  
  #total number of proteins
  total<-dim(x)[2]
  
  datas2<-merge(y,x,by.y = "row.names",by.x = "row.names")%>% dplyr::select(-Row.names)
  y<-as.factor(t(y$PROT_clust))
  datas2[1:6,1:6]
  datas2$PROT_clust<-as.factor(datas2$PROT_clust)
  #### Number of tree and repetition ----
  nr = 100
  nt = 1000
  
  
  #implementation of random forest classifier with nr repeats.
  imp<-matrix(0,(dim(x)[2]),nr)
  acc<-matrix(0,(dim(x)[2]),nr)
  #for each group separately
  for (i in paste0("g",1:numg)){
    assign(i,matrix(0,(dim(x)[2]),nr)) 
  }
  
  ta<-0
  set.seed(100)
  for (j in 1:nr){
    inds<-matrix(0,numg,1)
    #stratified hold-out cross validation (1 obs. per group)
    if (min(table(y))==1){ # v prípade že je v skupine iba jeden pacient nebudeme mať testovaciu sadu, iba trenovaciu - TA bude zla
      for (i in as.numeric(names(table(y)[table(y)!=1]))){ 
        
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) # v prípade že ide len o jednu sample a to je inteager  x >= 1, sampling via sample takes place from 1:x
        inds[i]<-ind  
      }
    }else {
      for (i in 1:numg){ 
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) 
        inds[i]<-ind
      }
    }
    
    
    borut<-randomForest::randomForest(datas2[-inds,-1],factor(datas2[-inds,1]),data=datas2[-inds,],ntree=nt,importance = T,replace = T)
    head(importance(borut))
    imp[,j]<-importance(borut)[,"MeanDecreaseGini"]
    acc[,j]<-importance(borut)[,"MeanDecreaseAccuracy"]
    
    #importance for each group
    # Assign importance for each group
    for (i in 1:numg) {
      matrix_name <- paste0("g", i)
      eval(parse(text = paste0(matrix_name, "[,j] <- importance(borut)[,", i, "]")))
    }
    
    #alternative - does not work
    
    #purrr::map(1:numg, ~ assign(paste0("g", .x, "[,j]"), importance(borut)[,.x]))
    
    
    
    ptest<-predict(borut,datas2[inds,-1])
    ta<-table(datas2[inds,1],ptest)+ta
  }
  means<-cbind(rowMeans(imp))
  
  for (i in 1:numg) {
    assign(paste0("g", i), cbind(rowMeans(get(paste0("g", i)))))
  }
  
  # g1<-cbind(rowMeans(g1))
  
  
  rownames(means)<-rownames(imp)<-colnames(datas2[,2:(dim(datas2)[2])])
  #calculate rows variances from MDI (Mean Decrease in Impurity)
  el<-as.data.frame(t(imp)) %>%  rstatix::get_summary_stats(type="mean_se")
  
  dt<-data.frame(el[c(1,3:4)])
  colnames(dt)<-c("MDI","Means","SE")
  imps<-dt[order(-dt$Means),]
  
  
  #Write in Excel ----
  for (i in 1:numg){
    g<-get(paste0("g",i))
    rownames(g)<-colnames(datas2[,2:(dim(datas2)[2])])
    #calculate rows variances from MDI (Mean Decrease in Impurity)
    g<-arrange(as.data.frame(g),desc(V1))
    colnames(g)<-levels(datas2$PROT_clust)[i]
    addWorksheet(wb,levels(datas2$PROT_clust)[i])
    writeData(wb,levels(datas2$PROT_clust)[i],head(g,500), startRow = 1, startCol = 1, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],"Pacienti v skupine", startRow = 1, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],names(y[y==i]), startRow = 2, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    
    assign(paste0("g",i),g)
  }
  
  # Save the results
  writeData(wb,"TNBC",head(imps,500), startRow = 1, startCol = 1, colNames = T, rowNames = F, headerStyle = hs1)  
  writeData(wb,"TNBC","Confusion matrix", startRow = 1, startCol = 5, colNames = T, rowNames = T, headerStyle = hs1) 
  writeData(wb,"TNBC",ta, startRow = 2, startCol = 5, colNames = T, rowNames = T) 
  
  print(head(as.data.frame(head(g3))))
  # Save the workbook
  saveWorkbook(wb,paste0("Res/Random_forests_imps_Feature",datum,".xlsx"),overwrite = T)
}

########################################################################
# Preselected features based on Feature selection TNBC.R---------

df <- readRDS("~/CEITEC/20241219_Bouchal_TNBC/Res/Filtered_Features_2025-02-27.rds")
head(df)
data %>% filter(gene_id %in% colnames(df)) %>% distinct(gene_id) %>% nrow()
data<-data %>% filter(gene_id %in% colnames(df))

{
  wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
  wb <- createWorkbook(Sys.setenv("R_ZIPCMD" = "c:/RBuildTools/3.5/bin/zip.exe"))
  
  #formaring style
  hs1 <- createStyle( halign = "CENTER", textDecoration = "Bold",border = "Bottom", fontColour = "black")
  #excel result worksheet
  addWorksheet(wb, "TNBC")
  
  
  #Redesign the data in a way that there will be individual lines on rows
  #and proteins  in columns. 
  
  
  #Only cinsider when all observation available
  x<-dcast(data, gene_id ~ sample, value.var = "PROT_fc")
  gene_id <- x$gene_id
  x <- (as.matrix(t(x[, -1, with = FALSE])))
  colnames(x) <- gene_id
  y<-as.data.frame(data %>% distinct(sample,PROT_clust))
  rownames(y)<-y$sample
  sum(rownames(x) %in% y$sample)
  y<-as.data.frame(y %>% dplyr::select(-sample))
  #check if the columns are the same
  
  
  
  
  #number of group
  numg<-uniqueN(y$PROT_clust)
  
  #total number of proteins
  total<-dim(x)[2]
  
  datas2<-merge(y,x,by.y = "row.names",by.x = "row.names")%>% dplyr::select(-Row.names)
  y<-as.factor(t(y$PROT_clust))
  datas2[1:6,1:6]
  datas2$PROT_clust<-as.factor(datas2$PROT_clust)
  #### Number of tree and repetition ----
  nr = 100
  nt = 1000
  
  
  #implementation of random forest classifier with nr repeats.
  imp<-matrix(0,(dim(x)[2]),nr)
  acc<-matrix(0,(dim(x)[2]),nr)
  #for each group separately
  for (i in paste0("g",1:numg)){
    assign(i,matrix(0,(dim(x)[2]),nr)) 
  }
  
  ta<-0
  set.seed(100)
  for (j in 1:nr){
    inds<-matrix(0,numg,1)
    #stratified hold-out cross validation (1 obs. per group)
    if (min(table(y))==1){ # v prípade že je v skupine iba jeden pacient nebudeme mať testovaciu sadu, iba trenovaciu - TA bude zla
      for (i in as.numeric(names(table(y)[table(y)!=1]))){ 
        
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) # v prípade že ide len o jednu sample a to je inteager  x >= 1, sampling via sample takes place from 1:x
        inds[i]<-ind  
      }
    }else {
      for (i in 1:numg){ 
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) 
        inds[i]<-ind
      }
    }
    
    
    borut<-randomForest::randomForest(datas2[-inds,-1],factor(datas2[-inds,1]),data=datas2[-inds,],ntree=nt,importance = T,replace = T)
    head(importance(borut))
    imp[,j]<-importance(borut)[,"MeanDecreaseGini"]
    acc[,j]<-importance(borut)[,"MeanDecreaseAccuracy"]
    
    #importance for each group
    # Assign importance for each group
    for (i in 1:numg) {
      matrix_name <- paste0("g", i)
      eval(parse(text = paste0(matrix_name, "[,j] <- importance(borut)[,", i, "]")))
    }
    

    ptest<-predict(borut,datas2[inds,-1])
    ta<-table(datas2[inds,1],ptest)+ta
  }
  means<-cbind(rowMeans(imp))
  
  for (i in 1:numg) {
    assign(paste0("g", i), cbind(rowMeans(get(paste0("g", i)))))
  }
  
  # g1<-cbind(rowMeans(g1))
  
  
  rownames(means)<-rownames(imp)<-colnames(datas2[,2:(dim(datas2)[2])])
  #calculate rows variances from MDI (Mean Decrease in Impurity)
  el<-as.data.frame(t(imp)) %>%  rstatix::get_summary_stats(type="mean_se")
  
  dt<-data.frame(el[c(1,3:4)])
  colnames(dt)<-c("MDI","Means","SE")
  imps<-dt[order(-dt$Means),]
  
  
  #Write in Excel ----
  for (i in 1:numg){
    g<-get(paste0("g",i))
    rownames(g)<-colnames(datas2[,2:(dim(datas2)[2])])
    #calculate rows variances from MDI (Mean Decrease in Impurity)
    g<-arrange(as.data.frame(g),desc(V1))
    colnames(g)<-levels(datas2$PROT_clust)[i]
    addWorksheet(wb,levels(datas2$PROT_clust)[i])
    writeData(wb,levels(datas2$PROT_clust)[i],head(g,500), startRow = 1, startCol = 1, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],"Pacienti v skupine", startRow = 1, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],names(y[y==i]), startRow = 2, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    
    assign(paste0("g",i),g)
  }
  
  # Save the results
  writeData(wb,"TNBC",head(imps,500), startRow = 1, startCol = 1, colNames = T, rowNames = F, headerStyle = hs1)  
  writeData(wb,"TNBC","Confusion matrix", startRow = 1, startCol = 5, colNames = T, rowNames = T, headerStyle = hs1) 
  writeData(wb,"TNBC",ta, startRow = 2, startCol = 5, colNames = T, rowNames = T) 
  
  print(head(as.data.frame(head(g3))))
  # Save the workbook
  saveWorkbook(wb,paste0("Res/Random_forests_PreselectedFeatures",datum,".xlsx"),overwrite = T)
}


# Varianta pre Bouchala - vybraných 9 proteinov s biologickým významom ---------
# Load data

X9_protein <- read_excel("9_protein.xlsx")
head(X9_protein)

data<-fread("~/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Katka_files/results_2025/full_clust_tab.tsv")
data[,PROT_fc := prot_log - mean(prot_log), .(gene_id)]
data <- data[!(sample %in% c("P_010","P_078","P_037","P_079","P_093","P_092"))]

head(data)

# filter data to only include the 9 preselected proteins
data <-
  data %>% filter(gene_name %in% X9_protein$gene)

#kontrola
head(data)  
data %>% distinct(gene_name) %>% nrow()
data %>% distinct(sample) %>% nrow()
data %>% distinct(PROT_clust) %>% nrow()

#RF -
  library(randomForest)
### Autosaving in xlsx
library(openxlsx)
{
  wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
  wb <- createWorkbook(Sys.setenv("R_ZIPCMD" = "c:/RBuildTools/3.5/bin/zip.exe"))
  
  #formaring style
  hs1 <- createStyle( halign = "CENTER", textDecoration = "Bold",border = "Bottom", fontColour = "black")
  #excel result worksheet
  addWorksheet(wb, "TNBC")
  
  
  #Redesign the data in a way that there will be individual lines on rows
  #and proteins  in columns. 
  
  
  #Only cinsider when all observation available
  x<-dcast(data, gene_id ~ sample, value.var = "PROT_fc")
  gene_id <- x$gene_id
  x <- (as.matrix(t(x[, -1, with = FALSE])))
  colnames(x) <- gene_id
  y<-as.data.frame(data %>% distinct(sample,PROT_clust))
  rownames(y)<-y$sample
  sum(rownames(x) %in% y$sample)
  y<-as.data.frame(y %>% dplyr::select(-sample))
  #check if the columns are the same
  
  
  
  
  #number of group
  numg<-uniqueN(y$PROT_clust)
  
  #total number of proteins
  total<-dim(x)[2]
  
  datas2<-merge(y,x,by.y = "row.names",by.x = "row.names")%>% dplyr::select(-Row.names)
  y<-as.factor(t(y$PROT_clust))
  datas2[1:6,1:6]
  datas2$PROT_clust<-as.factor(datas2$PROT_clust)
  #### Number of tree and repetition ----
  nr = 100
  nt = 1000
  
  
  #implementation of random forest classifier with nr repeats.
  imp<-matrix(0,(dim(x)[2]),nr)
  acc<-matrix(0,(dim(x)[2]),nr)
  #for each group separately
  for (i in paste0("g",1:numg)){
    assign(i,matrix(0,(dim(x)[2]),nr)) 
  }
  
  ta<-0
  set.seed(100)
  for (j in 1:nr){
    inds<-matrix(0,numg,1)
    #stratified hold-out cross validation (1 obs. per group)
    if (min(table(y))==1){ # v prípade že je v skupine iba jeden pacient nebudeme mať testovaciu sadu, iba trenovaciu - TA bude zla
      for (i in as.numeric(names(table(y)[table(y)!=1]))){ 
        
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) # v prípade že ide len o jednu sample a to je inteager  x >= 1, sampling via sample takes place from 1:x
        inds[i]<-ind  
      }
    }else {
      for (i in 1:numg){ 
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) 
        inds[i]<-ind
      }
    }
    
    
    borut<-randomForest::randomForest(datas2[-inds,-1],factor(datas2[-inds,1]),data=datas2[-inds,],ntree=nt,importance = T,replace = T)
    head(importance(borut))
    imp[,j]<-importance(borut)[,"MeanDecreaseGini"]
    acc[,j]<-importance(borut)[,"MeanDecreaseAccuracy"]
    
    #importance for each group
    # Assign importance for each group
    for (i in 1:numg) {
      matrix_name <- paste0("g", i)
      eval(parse(text = paste0(matrix_name, "[,j] <- importance(borut)[,", i, "]")))
    }
    
    #alternative - does not work
    
    #purrr::map(1:numg, ~ assign(paste0("g", .x, "[,j]"), importance(borut)[,.x]))
    
    
    
    ptest<-predict(borut,datas2[inds,-1])
    ta<-table(datas2[inds,1],ptest)+ta
  }
  means<-cbind(rowMeans(imp))
  
  for (i in 1:numg) {
    assign(paste0("g", i), cbind(rowMeans(get(paste0("g", i)))))
  }
  
  # g1<-cbind(rowMeans(g1))
  
  
  rownames(means)<-rownames(imp)<-colnames(datas2[,2:(dim(datas2)[2])])
  #calculate rows variances from MDI (Mean Decrease in Impurity)
  el<-as.data.frame(t(imp)) %>%  rstatix::get_summary_stats(type="mean_se")
  
  dt<-data.frame(el[c(1,3:4)])
  colnames(dt)<-c("MDI","Means","SE")
  imps<-dt[order(-dt$Means),]
  
  
  #Write in Excel ----
  for (i in 1:numg){
    g<-get(paste0("g",i))
    rownames(g)<-colnames(datas2[,2:(dim(datas2)[2])])
    #calculate rows variances from MDI (Mean Decrease in Impurity)
    g<-arrange(as.data.frame(g),desc(V1))
    colnames(g)<-levels(datas2$PROT_clust)[i]
    addWorksheet(wb,levels(datas2$PROT_clust)[i])
    writeData(wb,levels(datas2$PROT_clust)[i],head(g,500), startRow = 1, startCol = 1, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],"Pacienti v skupine", startRow = 1, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],names(y[y==i]), startRow = 2, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    
    assign(paste0("g",i),g)
  }
  
  # Save the results
  writeData(wb,"TNBC",head(imps,500), startRow = 1, startCol = 1, colNames = T, rowNames = F, headerStyle = hs1)  
  writeData(wb,"TNBC","Confusion matrix", startRow = 1, startCol = 5, colNames = T, rowNames = T, headerStyle = hs1) 
  writeData(wb,"TNBC",ta, startRow = 2, startCol = 5, colNames = T, rowNames = T) 
  
  print(head(as.data.frame(head(g3))))
  # Save the workbook
  saveWorkbook(wb,paste0("Res/Random_forests_9p_",datum,".xlsx"),overwrite = T)
}


# CRosss tab variant ----
# Pacient vs zvysok, ten z 9 proteinov ktory bude mat max fold change ten = clustru.
# Tzn ma z bouchal subor zaujima iba cluster = protein

#nesedia mi hodnoty ?
data[, PROT_group_fc_2 := prot_log - mean(prot_log), by = .(gene_id, PROT_clust )]%>% head()
data[, PROT_fc_2 := prot_log - mean(prot_log), by = .(gene_id )]%>% head()
data%>% dplyr::select(gene_id,PROT_clust,PROT_fc,PROT_fc_2,PROT_group_fc,PROT_group_fc_2)%>% head()

library(readxl)
clusters <- read_excel("241205_9_Best_targets.xlsx")
head(clusters)


crs_df<-data%>% dplyr::filter(data$protein_name %in% clusters$UniprotID)
crs_df%>% arrange(sample,desc(PROT_fc))%>% head(9)

crs<-crs_df %>%
  group_by(sample) %>%
  filter(PROT_fc == max(PROT_fc)) %>%
  select(gene_id,protein_name, PROT_clust) %>% full_join(clusters%>% dplyr::select(UniprotID,cluster), by = c("protein_name" = "UniprotID"))%>%
  mutate(Bouchal_clust=word(cluster,2,sep="_"))
crs%>% filter(Bouchal_clust==7)
table(crs$PROT_clust,crs$Bouchal_clust)


# ploty pre týchto 9 proteinov a 7 skupin based on Blouchal_clust a based on 
df_long <- crs_df %>%
  pivot_longer(
    cols = starts_with("Protein"),   # adjust if your column names differ
    names_to = "Protein", 
    values_to = "Intensity"
  )

# Create boxplots
ggplot(crs_df, aes(x =protein_name  , y = PROT_fc, fill = protein_name)) +
  geom_boxplot() +
  facet_wrap(~ PROT_clust, scales = "free_y") +  # one facet per protein intensity column
  theme_bw() +
  labs(
    title = "Protein Intensities by Group",
    x = "Group",
    y = "Protein Intensity"
  )


mean_cl_normal <- function(x) {
  m <- mean(x)
  se <- sd(x) / sqrt(length(x))
  data.frame(y = m, ymin = m - 1.96 * se, ymax = m + 1.96 * se)
}



ggplot(crs_df, aes(x = protein_name, y = PROT_fc, fill = protein_name)) +
  # Create bars for the mean PROT_fc per protein_name
  stat_summary(fun = mean, geom = "bar", color = "black", position = "dodge") +
  # Add error bars for the 95% CI computed by mean_cl_normal()
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = "dodge") +
  facet_wrap(~ PROT_clust, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  labs(
    title = "Protein FC by cluster group",
    x = "",
    y = "Protein Fold Change (PROT_fc)",
    fill= "Protein"
  )


# Vybraných 9 proteinov s biologickým významom + tie selected by importance ---------
# Load data
library(readxl)
X9_protein <- read_excel("9_protein.xlsx")
head(X9_protein)

data<-fread("~/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/Katka_files/results_2025/full_clust_tab.tsv")
data[,PROT_fc := prot_log - mean(prot_log), .(gene_id)]
data <- data[!(sample %in% c("P_010","P_078","P_037","P_079","P_093","P_092"))]

head(data)

# filter data to only include the 9 preselected proteins
data <-
  data %>% filter(gene_name %in% X9_protein$gene| gene_id %in% features)

#kontrola
head(data)  
data %>% distinct(gene_name) %>% nrow()
data %>% distinct(sample) %>% nrow()
data %>% distinct(PROT_clust) %>% nrow()

#RF -
library(randomForest)
### Autosaving in xlsx
library(openxlsx)
{
  wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
  wb <- createWorkbook(Sys.setenv("R_ZIPCMD" = "c:/RBuildTools/3.5/bin/zip.exe"))
  
  #formaring style
  hs1 <- createStyle( halign = "CENTER", textDecoration = "Bold",border = "Bottom", fontColour = "black")
  #excel result worksheet
  addWorksheet(wb, "TNBC")
  
  
  #Redesign the data in a way that there will be individual lines on rows
  #and proteins  in columns. 
  
  
  #Only cinsider when all observation available
  x<-dcast(data, gene_id ~ sample, value.var = "PROT_fc")
  gene_id <- x$gene_id
  x <- (as.matrix(t(x[, -1, with = FALSE])))
  colnames(x) <- gene_id
  y<-as.data.frame(data %>% distinct(sample,PROT_clust))
  rownames(y)<-y$sample
  sum(rownames(x) %in% y$sample)
  y<-as.data.frame(y %>% dplyr::select(-sample))
  #check if the columns are the same
  
  
  
  
  #number of group
  numg<-uniqueN(y$PROT_clust)
  
  #total number of proteins
  total<-dim(x)[2]
  
  datas2<-merge(y,x,by.y = "row.names",by.x = "row.names")%>% dplyr::select(-Row.names)
  y<-as.factor(t(y$PROT_clust))
  datas2[1:6,1:6]
  datas2$PROT_clust<-as.factor(datas2$PROT_clust)
  #### Number of tree and repetition ----
  nr = 100
  nt = 1000
  
  
  #implementation of random forest classifier with nr repeats.
  imp<-matrix(0,(dim(x)[2]),nr)
  acc<-matrix(0,(dim(x)[2]),nr)
  #for each group separately
  for (i in paste0("g",1:numg)){
    assign(i,matrix(0,(dim(x)[2]),nr)) 
  }
  
  ta<-0
  set.seed(100)
  for (j in 1:nr){
    inds<-matrix(0,numg,1)
    #stratified hold-out cross validation (1 obs. per group)
    if (min(table(y))==1){ # v prípade že je v skupine iba jeden pacient nebudeme mať testovaciu sadu, iba trenovaciu - TA bude zla
      for (i in as.numeric(names(table(y)[table(y)!=1]))){ 
        
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) # v prípade že ide len o jednu sample a to je inteager  x >= 1, sampling via sample takes place from 1:x
        inds[i]<-ind  
      }
    }else {
      for (i in 1:numg){ 
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) 
        inds[i]<-ind
      }
    }
    
    
    borut<-randomForest::randomForest(datas2[-inds,-1],factor(datas2[-inds,1]),data=datas2[-inds,],ntree=nt,importance = T,replace = T)
    head(importance(borut))
    imp[,j]<-importance(borut)[,"MeanDecreaseGini"]
    acc[,j]<-importance(borut)[,"MeanDecreaseAccuracy"]
    
    #importance for each group
    # Assign importance for each group
    for (i in 1:numg) {
      matrix_name <- paste0("g", i)
      eval(parse(text = paste0(matrix_name, "[,j] <- importance(borut)[,", i, "]")))
    }
    
    #alternative - does not work
    
    #purrr::map(1:numg, ~ assign(paste0("g", .x, "[,j]"), importance(borut)[,.x]))
    
    
    
    ptest<-predict(borut,datas2[inds,-1])
    ta<-table(datas2[inds,1],ptest)+ta
  }
  means<-cbind(rowMeans(imp))
  
  for (i in 1:numg) {
    assign(paste0("g", i), cbind(rowMeans(get(paste0("g", i)))))
  }
  
  # g1<-cbind(rowMeans(g1))
  
  
  rownames(means)<-rownames(imp)<-colnames(datas2[,2:(dim(datas2)[2])])
  #calculate rows variances from MDI (Mean Decrease in Impurity)
  el<-as.data.frame(t(imp)) %>%  rstatix::get_summary_stats(type="mean_se")
  
  dt<-data.frame(el[c(1,3:4)])
  colnames(dt)<-c("MDI","Means","SE")
  imps<-dt[order(-dt$Means),]
  
  
  #Write in Excel ----
  for (i in 1:numg){
    g<-get(paste0("g",i))
    rownames(g)<-colnames(datas2[,2:(dim(datas2)[2])])
    #calculate rows variances from MDI (Mean Decrease in Impurity)
    g<-arrange(as.data.frame(g),desc(V1))
    colnames(g)<-levels(datas2$PROT_clust)[i]
    addWorksheet(wb,levels(datas2$PROT_clust)[i])
    writeData(wb,levels(datas2$PROT_clust)[i],head(g,500), startRow = 1, startCol = 1, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],"Pacienti v skupine", startRow = 1, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],names(y[y==i]), startRow = 2, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    
    assign(paste0("g",i),g)
  }
  
  # Save the results
  writeData(wb,"TNBC",head(imps,500), startRow = 1, startCol = 1, colNames = T, rowNames = F, headerStyle = hs1)  
  writeData(wb,"TNBC","Confusion matrix", startRow = 1, startCol = 5, colNames = T, rowNames = T, headerStyle = hs1) 
  writeData(wb,"TNBC",ta, startRow = 2, startCol = 5, colNames = T, rowNames = T) 
  
  print(head(as.data.frame(head(g3))))
  # Save the workbook
  saveWorkbook(wb,paste0("Res/Random_forests_9p_plusFS",datum,".xlsx"),overwrite = T)
}


# Výsledný model s tuningom podla bayesa ----
# Parametre z r scripu - Bayes_intereference_tuning
opt_result <- readRDS("~/OneDrive/Documents/CEITEC/20241219_Bouchal_TNBC/opt_result.rds")

# filter data to only include the 45 preselected proteins
data <-
  data %>% filter(gene_id %in% features)

#kontrola
head(data)  
data %>% distinct(gene_name) %>% nrow()
data %>% distinct(sample) %>% nrow()
data %>% distinct(PROT_clust) %>% nrow()

#RF -
library(randomForest)
### Autosaving in xlsx
library(openxlsx)
{
  wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
  wb <- createWorkbook(Sys.setenv("R_ZIPCMD" = "c:/RBuildTools/3.5/bin/zip.exe"))
  
  #formaring style
  hs1 <- createStyle( halign = "CENTER", textDecoration = "Bold",border = "Bottom", fontColour = "black")
  #excel result worksheet
  addWorksheet(wb, "TNBC")
  
  
  #Redesign the data in a way that there will be individual lines on rows
  #and proteins  in columns. 
  
  
  #Only cinsider when all observation available
  x<-dcast(data, gene_id ~ sample, value.var = "PROT_fc")
  gene_id <- x$gene_id
  x <- (as.matrix(t(x[, -1, with = FALSE])))
  colnames(x) <- gene_id
  y<-as.data.frame(data %>% distinct(sample,PROT_clust))
  rownames(y)<-y$sample
  sum(rownames(x) %in% y$sample)
  y<-as.data.frame(y %>% dplyr::select(-sample))
  #check if the columns are the same
  
  
  
  
  #number of group
  numg<-uniqueN(y$PROT_clust)
  
  #total number of proteins
  total<-dim(x)[2]
  
  datas2<-merge(y,x,by.y = "row.names",by.x = "row.names")%>% dplyr::select(-Row.names)
  y<-as.factor(t(y$PROT_clust))
  datas2[1:6,1:6]
  datas2$PROT_clust<-as.factor(datas2$PROT_clust)
  #### Number of tree and repetition ----
  nr = 100
  nt = 1000
  
  
  #implementation of random forest classifier with nr repeats.
  imp<-matrix(0,(dim(x)[2]),nr)
  acc<-matrix(0,(dim(x)[2]),nr)
  #for each group separately
  for (i in paste0("g",1:numg)){
    assign(i,matrix(0,(dim(x)[2]),nr)) 
  }
  
  ta<-0
  set.seed(100)
  for (j in 1:nr){
    inds<-matrix(0,numg,1)
    #stratified hold-out cross validation (1 obs. per group)
    if (min(table(y))==1){ # v prípade že je v skupine iba jeden pacient nebudeme mať testovaciu sadu, iba trenovaciu - TA bude zla
      for (i in as.numeric(names(table(y)[table(y)!=1]))){ 
        
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) # v prípade že ide len o jednu sample a to je inteager  x >= 1, sampling via sample takes place from 1:x
        inds[i]<-ind  
      }
    }else {
      for (i in 1:numg){ 
        ind<-sample(which(datas2[,1]==unique(datas2[,1])[i]),1) 
        inds[i]<-ind
      }
    }
    
    # ---- Handling Imbalanced Groups ----
    # Compute class weights from training set (datas2[-inds,])
    train_labels <- factor(datas2[-inds, 1])
    tbl <- table(train_labels)
    # Higher weight for underrepresented classes (weight = max(freq)/freq)
    class_weights <- as.vector(max(tbl) / tbl)
    names(class_weights) <- names(tbl)
    
    borut<-randomForest::randomForest(datas2[-inds,-1],factor(datas2[-inds,1]),data=datas2[-inds,],ntree=nt,importance = T,replace = T)
    
    borut <- randomForest::randomForest(
      x = datas2[-inds, -1],
      y = factor(datas2[-inds, 1]),
      data = datas2[-inds, ],
      ntree = nt,  # fewer trees for speed during optimization
      mtry = opt_result$Best_Par["mtry"],
      nodesize = opt_result$Best_Par["nodesize"],
      classwt = class_weights,
      importance = T,
      replace = TRUE
    )
    
    head(importance(borut))
    imp[,j]<-importance(borut)[,"MeanDecreaseGini"]
    acc[,j]<-importance(borut)[,"MeanDecreaseAccuracy"]
    
    #importance for each group
    # Assign importance for each group
    for (i in 1:numg) {
      matrix_name <- paste0("g", i)
      eval(parse(text = paste0(matrix_name, "[,j] <- importance(borut)[,", i, "]")))
    }
    
    
    
    ptest<-predict(borut,datas2[inds,-1])
    ta<-table(datas2[inds,1],ptest)+ta
  }
  means<-cbind(rowMeans(imp))
  
  for (i in 1:numg) {
    assign(paste0("g", i), cbind(rowMeans(get(paste0("g", i)))))
  }
  
  # g1<-cbind(rowMeans(g1))
  
  
  rownames(means)<-rownames(imp)<-colnames(datas2[,2:(dim(datas2)[2])])
  #calculate rows variances from MDI (Mean Decrease in Impurity)
  el<-as.data.frame(t(imp)) %>%  rstatix::get_summary_stats(type="mean_se")
  
  dt<-data.frame(el[c(1,3:4)])
  colnames(dt)<-c("MDI","Means","SE")
  imps<-dt[order(-dt$Means),]
  
  
  #Write in Excel ----
  for (i in 1:numg){
    g<-get(paste0("g",i))
    rownames(g)<-colnames(datas2[,2:(dim(datas2)[2])])
    #calculate rows variances from MDI (Mean Decrease in Impurity)
    g<-arrange(as.data.frame(g),desc(V1))
    colnames(g)<-levels(datas2$PROT_clust)[i]
    addWorksheet(wb,levels(datas2$PROT_clust)[i])
    writeData(wb,levels(datas2$PROT_clust)[i],head(g,500), startRow = 1, startCol = 1, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],"Pacienti v skupine", startRow = 1, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    writeData(wb,levels(datas2$PROT_clust)[i],names(y[y==i]), startRow = 2, startCol = 4, colNames = T, rowNames = T, headerStyle = hs1)  
    
    assign(paste0("g",i),g)
  }
  
  # Save the results
  writeData(wb,"TNBC",head(imps,500), startRow = 1, startCol = 1, colNames = T, rowNames = F, headerStyle = hs1)  
  writeData(wb,"TNBC","Confusion matrix", startRow = 1, startCol = 5, colNames = T, rowNames = T, headerStyle = hs1) 
  writeData(wb,"TNBC",ta, startRow = 2, startCol = 5, colNames = T, rowNames = T) 
  
  print(head(as.data.frame(head(g3))))
  # Save the workbook
  saveWorkbook(wb,paste0("Res/Random_forests_FS_tuned2",datum,".xlsx"),overwrite = T)
}

# FInal model based on this----------

df <- readRDS("~/CEITEC/20241219_Bouchal_TNBC/Res/Filtered_Features_2025-02-27.rds")
opt_result <- readRDS("~/CEITEC/20241219_Bouchal_TNBC/opt_result.rds")

#kontrola
head(data)  
data %>% distinct(gene_name) %>% nrow()
data %>% distinct(sample) %>% nrow()
data %>% distinct(PROT_clust) %>% nrow()

#RF -
library(randomForest)
### Autosaving in xlsx
library(openxlsx)

  
  
  #Redesign the data in a way that there will be individual lines on rows
  #and proteins  in columns. 
  
  
  #Only cinsider when all observation available
  x<-dcast(data, gene_id ~ sample, value.var = "PROT_fc")
  gene_id <- x$gene_id
  x <- (as.matrix(t(x[, -1, with = FALSE])))
  colnames(x) <- gene_id
  y<-as.data.frame(data %>% distinct(sample,PROT_clust))
  rownames(y)<-y$sample
  sum(rownames(x) %in% y$sample)
  y<-as.data.frame(y %>% dplyr::select(-sample))
  #check if the columns are the same
  
  
  
  
  #number of group
  numg<-uniqueN(y$PROT_clust)
  
  #total number of proteins
  total<-dim(x)[2]
  
  datas2<-merge(y,x,by.y = "row.names",by.x = "row.names")%>% dplyr::select(-Row.names)
  y<-as.factor(t(y$PROT_clust))
  datas2[1:6,1:6]
  datas2$PROT_clust<-as.factor(datas2$PROT_clust)
  #### Number of tree and repetition ----
  nt = 1000
  
  set.seed(100)

    
    # ---- Handling Imbalanced Groups ----
    # Compute class weights from training set (datas2[-inds,])
    train_labels <- factor(datas2[, 1])
    tbl <- table(train_labels)
    # Higher weight for underrepresented classes (weight = max(freq)/freq)
    class_weights <- as.vector(max(tbl) / tbl)
    names(class_weights) <- names(tbl)
    
    final_model <- randomForest::randomForest(
      x = datas2[, -1],
      y = factor(datas2[, 1]),
      data = datas2[, ],
      ntree = nt,  # fewer trees for speed during optimization
      mtry =  opt_result$Best_Par["mtry"], #nemusí byť harcodované on si to zapamatá
      nodesize = opt_result$Best_Par["nodesize"],
      classwt = class_weights,
      importance = T,
      replace = TRUE
    )
    
    saveRDS(final_model, "final_model_TNBC.rds")
    
    # Install and load the pmml package if not already installed
    install.packages("r2pmml")
    library(r2pmml)
    
    # Export the final model to a PMML file
    r2pmml(final_model, file = "final_model_TNBC.pmml")
    
# CRosss tab variant ----
# Pacient vs zvysok, ten z 9 proteinov ktory bude mat max fold change ten = clustru.
# Tzn ma z bouchal subor zaujima iba cluster = protein


library(readxl)
clusters <- read_excel("241205_9_Best_targets.xlsx")
head(clusters)


crs_df<-data%>% dplyr::filter(data$protein_name %in% clusters$UniprotID)
crs_df%>% arrange(sample,desc(PROT_fc))%>% head(9)

crs<-crs_df %>%
  group_by(sample) %>%
  filter(PROT_fc == max(PROT_fc)) %>%
  select(gene_id,protein_name, PROT_clust) %>% full_join(clusters%>% dplyr::select(UniprotID,cluster), by = c("protein_name" = "UniprotID"))%>%
  mutate(Bouchal_clust=stringr::word(cluster,2,sep="_"))
crs%>% filter(Bouchal_clust==7)
table(crs$PROT_clust,crs$Bouchal_clust)


# ploty pre týchto 9 proteinov a 7 skupin based on Blouchal_clust a based on 
df_long <- crs_df %>%
  tidyr::pivot_longer(
    cols = starts_with("Protein"),   # adjust if your column names differ
    names_to = "Protein", 
    values_to = "Intensity"
  )

# Create boxplots
ggplot(crs_df, aes(x =protein_name  , y = PROT_fc, fill = protein_name)) +
  geom_boxplot() +
  facet_wrap(~ PROT_clust, scales = "free_y") +  # one facet per protein intensity column
  theme_bw() +
  labs(
    title = "Protein Intensities by Group",
    x = "Group",
    y = "Protein Intensity"
  )


mean_cl_normal <- function(x) {
  m <- mean(x)
  se <- sd(x) / sqrt(length(x))
  data.frame(y = m, ymin = m - 1.96 * se, ymax = m + 1.96 * se)
}



ggplot(crs_df, aes(x = protein_name, y = PROT_fc, fill = protein_name)) +
  # Create bars for the mean PROT_fc per protein_name
  stat_summary(fun = mean, geom = "bar", color = "black", position = "dodge") +
  # Add error bars for the 95% CI computed by mean_cl_normal()
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2, position = "dodge") +
  facet_wrap(~ PROT_clust, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  labs(
    title = "Protein FC by cluster group",
    x = "",
    y = "Protein Fold Change (PROT_fc)",
    fill= "Protein"
  )

ggsave("Pics/Prote9_byCluster.jpg")
