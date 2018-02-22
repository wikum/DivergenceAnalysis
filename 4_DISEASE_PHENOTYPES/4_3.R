
sink("log/4_3.log.txt", split=TRUE)

tryCatch({

  source("../src/util.R")
  source("util_4.R")
  
  source("../vars.R") # load DATA_DIR
  
  # ====================================================
  # ADENOMA CARCINOMA
  # ====================================================

  cat("Loading data..\n")
  
  dataMat = readTable(sprintf("%s/COLON/MICROARRAY_GPL570/GSE20916/matGSE20916.csv.gz", DATA_DIR))
  Pheno = readTable(sprintf("%s/COLON/MICROARRAY_GPL570/GSE20916/phenoGSE20916.csv.gz", DATA_DIR))
  
  rownames(dataMat) = as.character(dataMat[, 1])
  dataMat = data.matrix(dataMat[, -1])
  
  if(! all(is.finite(dataMat))){
    cat("Removing rows with missing values..\n")
    dataMat = dataMat[-which(apply(dataMat, 1, function(x) ! all(is.finite(x)))), ] #8267
  }
  
  all(Pheno$SAMPLE == colnames(dataMat))
  
  baseMat = dataMat[, which(Pheno$TYPE == "NORMAL")]
  sel = which(Pheno$TYPE %in% c("ADENOMA", "CARCINOMA"))
  Mat = dataMat[, sel]
  Groups = factor(Pheno$TYPE[sel], levels=c("ADENOMA", "CARCINOMA"))
  
  # ================ compute divergence ================ 
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat,
                                      gamma = 1:9/10)  
  
  df = data.frame(N=div$div$count.div,
                  Groups=Groups, 
                  GroupsN=make_n_factor(Groups))
  
  dfMean = get_mean_df(df)
  
  # ====== save ======
  save(df, dfMean, file="obj/4_3.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()  
  
}, error = function(e){ print(e) })

sink()






