
sink("log/4_4.log.txt", split=TRUE)

tryCatch({

  
  set.seed(1)
  
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(cowplot)
  library(RColorBrewer)
  library(ggrepel)
  library(plyr)
  
  currentPath = getwd()
  
  setwd("../../../../../PROBLEMS/code/")
  
  source("util.R")
  source("plotutil.R")
  source("divergence.R")
  
  setwd(currentPath)
  
  source("../../../code/temp_util.R")

  # ====================================================
  # BRESAT, RNA-SEQ
  # ====================================================
  
  cat("\nBREAST, RNA-SEQ\n")
  
  Mat = data.matrix(readTable("../../../../../PROBLEMS/BREAST/PROCESSED/A/MatData.csv.gz"))
  refMat = data.matrix(readTable("../../../../../PROBLEMS/BREAST/PROCESSED/A/MatRef.csv.gz"))
  Pheno = readTable("../../../../../DATA/BREAST/RNASeq/TCGA/TCGA_Pheno.csv.gz")
  rownames(Pheno) = Pheno$sample
  
  Pheno = Pheno[colnames(Mat), ]
  
  sel = which(Pheno$sample_type %in% c("Primary Tumor"))
  
  Mat = Mat[, sel]
  Pheno = Pheno[sel, ]
  
  MatP1 = getPercentileMat(Mat)
  refMatP1 = getPercentileMat(refMat)
  
  rm(Mat)
  rm(refMat)
  
  gc()
  
  cat("\nBREAST, Methylation\n")
  
  Mat2 = data.matrix(readTable("../../../../../PROBLEMS/BREAST/PROCESSED/B/MatData.csv"))
  refMat2 = data.matrix(readTable("../../../../../PROBLEMS/BREAST/PROCESSED/B/MatRef.csv"))

  sel = intersect(colnames(MatP1), colnames(Mat2))
  
  Mat2 = Mat2[, sel]
  MatP1 = MatP1[, sel]
  
  Pheno = Pheno[sel, ]

  # ================ compute divergence ================ 
  div1 = computeDivergences(Mat=MatP1, refMat=refMatP1, 
                           upper="perc.high", lower="perc.low")
  
  div2 = computeDivergences(Mat=Mat2, refMat=refMat2, 
                            upper="perc.high", lower="perc.low")
  
  save(div1, div2, file="4_4_data.rda")
  
  df = data.frame(
    N_RNASEQ=div1$N,
    N_METHYL=div2$N,
    #normalized
    N_RNASEQ_n=div1$N/nrow(refMatP1),
    N_METHYL_n=div2$N/nrow(refMat2),
    SAMPLE=colnames(MatP1)
  )
  
  cat("\nSpearman Correlation: %.3f\n", cor(df$N_RNASEQ, df$N_METHYL, method="spearman"))
  cat("\nPearson Correlation: %.3f\n", cor(df$N_RNASEQ, df$N_METHYL, method="pearson"))
  
  # ====================================================
  # ====================================================
  
  save(df, file="4_4.rda")  
  
}, error = function(e){ print(e) })

sink()






