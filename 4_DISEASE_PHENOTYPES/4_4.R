
sink("log/4_4.log.txt", split=TRUE)

tryCatch({
  
  set.seed(1)
  
  library(plyr)
  library(rutils)
  library(divergence)
  
  source("../src/util.R")
  
  source("../vars.R") # load DATA_DIR
  
  # ====================================================
  # BRESAT, RNA-SEQ
  # ====================================================
  
  cat("Loading data..\n")
  
  dataMat = data.matrix(readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_TPM.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
  
  baseMat1 = dataMat[, which(dataPheno$sample_type %in% c("Solid Tissue Normal"))]
  Mat1 = dataMat[, which(dataPheno$sample_type %in% c("Primary Tumor"))]
  
  rm(dataMat, dataPheno)
  
  dataMat = data.matrix(readTable(sprintf("%s/BREAST/METHYLATION_450k/TCGA/TCGA_Beta.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/BREAST/METHYLATION_450k/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
  
  if(! all(is.finite(dataMat))){
    cat("Removing rows with missing values..\n")
    dataMat = dataMat[-which(apply(dataMat, 1, function(x) ! all(is.finite(x)))), ] #364130
  }
  
  baseMat2 = dataMat[, which(dataPheno$sample_type %in% c("Solid Tissue Normal"))]
  Mat2 = dataMat[, which(dataPheno$sample_type %in% c("Primary Tumor"))]
  
  rm(dataMat, dataPheno)
  
  common = intersect(colnames(Mat1), colnames(Mat2))  
  
  Mat1 = Mat1[, common]
  Mat2 = Mat2[, common]  
  
  # ================ compute divergence ================ 
  div1 = computeUnivariateDigitization(Mat=Mat1, baseMat=baseMat1,
                                       computeQuantiles = TRUE)
  
  div2 = computeUnivariateDigitization(Mat=Mat2, baseMat=baseMat2,
                                       computeQuantiles = FALSE)
  
  df = data.frame(
    N_RNASEQ=div1$div$count.div,
    N_METHYL=div2$div$count.div,
    #normalized
    N_RNASEQ_n=div1$div$count.div/nrow(Mat1),
    N_METHYL_n=div2$div$count.div/nrow(Mat2),
    SAMPLE=colnames(Mat1)
  )
  
  cat(sprintf("\nSpearman Correlation: %.3f\n", cor(df$N_RNASEQ, df$N_METHYL, method="spearman")))
  cat(sprintf("\nPearson Correlation: %.3f\n", cor(df$N_RNASEQ, df$N_METHYL, method="pearson")))
  
  # ====== save ======
  save(df, file="obj/4_4.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()  

  
}, error = function(e){ print(e) })

sink()






