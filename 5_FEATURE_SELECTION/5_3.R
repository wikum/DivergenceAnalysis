
sink("log/5_3.log.txt", split=TRUE)

tryCatch({

  set.seed(1)
  
  library(plyr)
  library(divergence)
  
  source("../src/util.R")
  
  source("../vars.R") # load DATA_DIR

  # ====================================================
  # GTEx
  # ====================================================
  
  cat("Loading data..\n")
  #Mat = data.matrix(readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Counts_normalized.csv", DATA_DIR)))
  load("GTEX_MatP.rda")
  
  Pheno = readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Pheno.csv", DATA_DIR))
  rownames(Pheno) = Pheno$sample
  
  tissues = unique(as.character(Pheno$tissue_type))
  # omit the "" tissue type (i.e. two samples with blank tissue types)
  tissues = tissues[nchar(tissues) > 0]
  
  if(FALSE){
    if(! all(is.finite(Mat))){
      cat("Removing rows with missing values..\n")
      Mat = Mat[-which(apply(Mat, 1, function(x) ! all(is.finite(x)))), ]
    }
    
    MatP = computeQuantileMatrix(Mat)
  }
  
  sel = which(Pheno$tissue_type %in% tissues)
  
  MatP = MatP[, sel]
  Pheno = Pheno[sel, ]
  
  sel_tissues = setdiff(tissues, "Lung")
  
  # ====================================================
  # train on Lung
  # ====================================================
  
  train_ids =   sample(which(Pheno$tissue_type == "Lung"), ceiling(length(which(Pheno$tissue_type == "Lung")) * 0.5 ))
  
  Groups = factor(ifelse(Pheno$tissue_type[-train_ids] == "Lung", "LUNG", "OTHER"))

  # ================ compute divergence ================ 
  div = computeUnivariateDigitization(Mat=MatP[, -train_ids], baseMat=MatP[, train_ids],
                                      computeQuantiles = FALSE,  gamma=1:9/10,
                                      Groups=Groups)
  
  cat("Computing Chi-squared tests\n")
  
  C = computeChiSquaredTest(Mat=div$Mat.div, Groups=Groups, classes=c("LUNG", "OTHER"))  
  
  C$bon = p.adjust(C$pval, method="bon")
  
  df = data.frame(
    gene=rownames(div$Mat.div),
    PA=rowMeans(abs(div$Mat.div)[, which(Groups=="LUNG")]),
    PB=rowMeans(abs(div$Mat.div)[, which(Groups=="OTHER")]),
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
  save(df, file="obj/5_3.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()    
  
}, error = function(e){ print(e) })

sink()






