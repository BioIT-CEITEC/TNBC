### Klasifikátor - based on PROTEINY ?
 rm(list=ls())
# Load libraries
library(dplyr)
library(data.table)
 library(randomForest)
 ### Autosaving in xlsx
 library(openxlsx)

 #### IMPUTS ----------
 dir_im<-"TNBC/Classifier"
 data_im <- "full_clust_tab.tsv"
 df_im <- "Filtered_Features.rds"
 param_im<-"opt_result.rds"
 
setwd(dir_im) 
wd<-getwd()
datum<-Sys.Date()


# Load data   
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

df <- readRDS(df_im)
opt_result <- readRDS(param_im)

data %>% filter(gene_id %in% colnames(df)) %>% distinct(gene_id) %>% nrow()
data<-data %>% filter(gene_id %in% colnames(df))

data %>% distinct(gene_name) %>% nrow()
data %>% distinct(sample) %>% nrow()
data %>% distinct(PROT_clust) %>% nrow()

# RF ------
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
  set.seed(100)
  
  
  
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
  saveWorkbook(wb,paste0("Random_forests_",datum,".xlsx"),overwrite = T)
}

# Setting up final model
  
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
      mtry =  opt_result$Best_Par["mtry"], 
      nodesize = opt_result$Best_Par["nodesize"],
      classwt = class_weights,
      importance = T,
      replace = TRUE
    )
    
    saveRDS(final_model, "final_model_TNBC.rds")
  
