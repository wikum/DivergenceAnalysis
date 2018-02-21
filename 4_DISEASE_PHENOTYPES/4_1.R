
sink("log/4_1.log.txt", split=TRUE)

tryCatch({

  set.seed(1)
  
  library(plyr)
  library(rutils)
  library(divergence)
  
  source("../src/util.R")
  
  source("../vars.R") # load DATA_DIR
  
  # ====================================================
  # PROSTATE, RNA-SEQ
  # ====================================================
  
  cat("Loading data..\n")
  
  dataMat = data.matrix(readTable(sprintf("%s/PROSTATE/RNASeq/TCGA/TCGA_TPM.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/PROSTATE/RNASeq/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
  
  baseMat1 = dataMat[, which(dataPheno$sample_type %in% c("Solid Tissue Normal"))]
  Mat1 = dataMat[, which(dataPheno$sample_type %in% c("Primary Tumor"))]

  rm(dataMat, dataPheno)
  
  dataMat = data.matrix(readTable(sprintf("%s/PROSTATE/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Beta.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/PROSTATE/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Pheno.csv.gz", DATA_DIR))
  
  baseMat2 = dataMat[, which(dataPheno$sample_type %in% c("Solid Tissue Normal"))]
  Mat2 = dataMat[, which(dataPheno$sample_type %in% c("Primary Tumor"))]
  
  rm(dataMat, dataPheno)
  
  common = intersect(colnames(Mat1), colnames(Mat2))  
  
  Mat1 = Mat1[, common]
  Mat2 = Mat2[, common]
  
  # 497 samples
  
  if(! all(is.finite(Mat2))){
    cat("Removing rows with missing values..\n")
    Mat2 = Mat2[-which(apply(Mat2, 1, function(x) ! all(is.finite(x)))), ]
  }
  if(! all(is.finite(baseMat2))){
    cat("Removing rows with missing values from baseline data..\n")
    baseMat2 = baseMat2[-which(apply(baseMat2, 1, function(x) ! all(is.finite(x)))), ]
  }
  
  common_cpg = intersect(rownames(Mat2), rownames(baseMat2))
  Mat2 = Mat2[common_cpg, ]
  baseMat2 = baseMat2[common_cpg, ]
  
  # > dim(Mat1)
  # [1] 20531   497
  # > dim(Mat2)
  # [1] 381192    497
  
  #> nrow(refMatP1)
  #[1] 20531
  #> nrow(refMat2)
  #[1] 381182
  
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
  save(df, file="obj/4_1.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()
  
}, error = function(e){ print(e) })

sink()






