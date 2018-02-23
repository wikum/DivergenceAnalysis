

sink("log/6_1.log.txt", split=TRUE)

tryCatch({
  
  set.seed(1)
  
  library(plyr)
  library(divergence)
  
  source("../src/util.R")
  
  source("../vars.R") # load DATA_DIR

  HallMarks <- data.matrix(read.csv(sprintf("%s/GENESETS/HALLMARKS.txt", DATA_DIR), 
                      head=TRUE, sep=" ",na.strings="NA",stringsAsFactors = FALSE))

  cat("Loading breast data\n")
  dataMat = data.matrix(readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_TPM.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))

  cat("Making groups\n")
  sel_n = which(dataPheno$sample_type == "Solid Tissue Normal")
  baseMat = dataMat[, sel_n]

  sel_t = list(
    which(dataPheno$PAM50_mRNA_nature2012 == "Luminal A"),
    which(dataPheno$PAM50_mRNA_nature2012 == "Luminal B"),
    which(dataPheno$PAM50_mRNA_nature2012 == "HER2-enriched"),
    which(dataPheno$PAM50_mRNA_nature2012 == "Basal-like"),
  
    which(dataPheno$ER_Status_nature2012 == "Positive"),
    which(dataPheno$ER_Status_nature2012 == "Negative")
  )

  sel_t_ids = unlist(sel_t)

  new_labels = c('LuminalA', 'LuminalB', 'HER2-E', 'Basal', 'Positive', 'Negative')

  groups_v = unlist(lapply(1:length(sel_t), function(j) rep(new_labels[j], length(sel_t[[j]]))))

  Mat = dataMat[, sel_t_ids]

  Groups = factor(groups_v, levels=new_labels)

  Pheno = data.frame(sample=colnames(Mat), group=Groups)

  cat("Computing quantiles\n")
  MatNormal = computeQuantileMatrix(baseMat)
  MatCancer = computeQuantileMatrix(Mat)
  
  
  cat("Running experiment..\n")

  #### use 1/2 normal samples to build support, and other 1/2 for test
  normtrain <- sample(colnames(MatNormal), size=floor(1/2*ncol(MatNormal)),replace=FALSE)
  normtest <- setdiff(colnames(MatNormal), normtrain)
  MatNormalTrain <- MatNormal[, normtrain]
  MatNormalTest  <- MatNormal[,normtest]


  #### find the optimal gamma and corresponding support

  FGAS <- findMultivariateGammaAndSupport(Mat=MatNormalTrain, genesets = HallMarks,
                                    gamma = seq(0.4, 0.99,0.02),beta =0.95, alpha= 0.01,method="euclidean")

  ### transform test data and cancer data into binary expression

  centermatrix <- FGAS$MultiSetSupport$Centermatrix_list
  radius <- FGAS$MultiSetSupport$Radius_list


  Binary_normTest <- computeMultivariateBinaryMatrix(Mat=MatNormalTest, genesets = HallMarks,
                                           Centermatrix_list = centermatrix , Radius_list = radius, method="euclidean" )
  Binary_Cancer <- computeMultivariateBinaryMatrix(Mat= MatCancer, genesets = HallMarks, 
                                             Centermatrix_list = centermatrix , Radius_list = radius, method="euclidean")


  #### output optimal gamma

  opt_gamma <- FGAS$optimal_gamma



  ##### boxplot for optimal gamma 

  PhenoNormTest <- data.frame(cbind(colnames(MatNormalTest), 'Normal'))
  colnames(PhenoNormTest)<- c('sample','group')


  MATALL2<- cbind(Binary_normTest, Binary_Cancer)
  PHENO2<- rbind(PhenoNormTest,Pheno)


  ### evaluate divergent pathway numbers for each sample
  Div_Test <- apply(MATALL2,2, sum)

  df = data.frame(sample=PHENO2$sample, group=PHENO2$group, div=Div_Test)
  df$group = factor(df$group,
                    levels=c("Normal", "LuminalA", "LuminalB", "HER2-E", "Basal", 
                             "Positive", "Negative"),
                    labels=c("Normal", "Luminal A", "Luminal B", "HER2-enriched", "Basal-like",
                             "ER+", "ER-"))
  df$set = rep("NORMAL", nrow(df))
  df$set[which(df$group %in% c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like"))] = "PAM50\nSubtype"
  df$set[which(df$group %in% c("ER+", "ER-"))] = "ER\nStatus"
  df$set = factor(df$set, levels=c("NORMAL", "PAM50\nSubtype", "ER\nStatus"))
  
  df$group_n = paste(df$group, make_n_factor(df$group), sep="\n")
  
  # ====== save ======
  save(df, file="obj/6_1.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()   
  
}, error = function(e){ print(e) })

sink()






