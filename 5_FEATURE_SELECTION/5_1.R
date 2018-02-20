
sink("log/5_1.log.txt", split=TRUE)

tryCatch({

  set.seed(1)
  
  library(plyr)
  library(divergence)
  
  source("../src/util.R")
  
  source("vars.R") # load DATA_DIR

  # ====================================================
  # BREAST, RNA-SEQ, Luminal A vs Luminal B
  # ====================================================
  
  cat("Loading data..\n")
  
  dataMat = data.matrix(readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_TPM.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
  
  sel_n = which(dataPheno$sample_type == "Solid Tissue Normal")
  baseMat = dataMat[, sel_n]
    
  sel_t = which(dataPheno$PAM50_mRNA_nature2012 %in% c("Luminal A", "Luminal B"))
  
  Mat = dataMat[, sel_t]
  Pheno = dataPheno[sel_t, ]
  
  Groups = factor(Pheno$PAM50_mRNA_nature2012,
                  levels=c("Luminal A", "Luminal B"),
                  labels=c("LuminalA", "LuminalB"))
  
  # ================ compute divergence ================ 
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat,
                                      Groups=Groups, classes=c("LuminalA", "LuminalB"))
  
  save(div, Groups, file="obj/5_1_data.rda")
  
  cat("Computing Chi-squared tests\n")
  
  C = computeChiSquaredTest(Mat=div$Mat.div, Groups=Groups, classes=c("LuminalA", "LuminalB"))
  
  C$bon = p.adjust(C$pval, method="bon")
  
  df = data.frame(
    gene=rownames(div$Mat.div),
    PA=rowMeans(abs(div$Mat.div)[, which(Groups=="LuminalA")]),
    PB=rowMeans(abs(div$Mat.div)[, which(Groups=="LuminalB")]),
    C[rownames(div$Mat.div), ]
  )
  
  df$sig = factor(df$bon <= 0.05)
  
  df = df[is.finite(df$statistic), ]
  
  r = df$statistic * sign(df$PB - df$PA)
  
  r1 = rank(r)
  r2 = rank(-r)
  
  s = mapply(min, r1, r2)
  
  df$ssig = s <=10 & df$sig==TRUE
  
  # ====== save ======
  save(df, file="obj/5_1.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()    
  
}, error = function(e){ print(e) })

sink()






