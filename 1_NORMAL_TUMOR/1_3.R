
sink("log/1_3.log.txt", split=TRUE)

tryCatch({
  
  library(knitr)

  source("../src/util.R")

  source("util_1.R")
  source("../vars.R")
  
  # ====================================================
  # GPL570
  # ====================================================
  
  cat("Loading gpl570 data..\n")
  
  Mat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL570/TUMOR_AGGREGATED/TUMOR_Expression.csv.gz", DATA_DIR)))
  baseMat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL570/NORMAL_AGGREGATED/NORMAL_Expression.csv.gz", DATA_DIR)))
  
  Pheno = readTable(sprintf("%s/BREAST/MICROARRAY_GPL570/TUMOR_AGGREGATED/TUMOR_Pheno.csv.gz", DATA_DIR))
  basePheno = readTable(sprintf("%s/BREAST/MICROARRAY_GPL570/NORMAL_AGGREGATED/NORMAL_Pheno.csv.gz", DATA_DIR))
  
  common = intersect(rownames(baseMat), rownames(Mat))
  
  Mat = Mat[common, ] #11806
  baseMat = baseMat[common, ]
  
  sel = sample(1:ncol(baseMat), ceiling(ncol(baseMat) * 0.5))
  
  Groups = factor(c(rep("TUMOR", ncol(Mat)), rep("NORMAL", ncol(baseMat[, -sel]))))
  Mat = cbind(Mat, baseMat[, -sel])
  baseMat = baseMat[, sel]
  
  # ================ compute divergence ================ 
  
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat, 
                                      gamma = 1:9/10,
                                      Groups=Groups, classes=c("NORMAL", "TUMOR"))
  
  df_GPL570 = data.frame(N=div$div$count.div, Groups=Groups)
  
  rm(Mat, baseMat, div, Pheno, Groups)
  
  # ====================================================
  # GPL96
  # ====================================================
  
  cat("Loading gpl96 data..\n")
  
  Mat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL96/TUMOR_AGGREGATED/TUMOR_Expression.csv.gz", DATA_DIR)))
  baseMat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL96/NORMAL_AGGREGATED/NORMAL_Expression.csv.gz", DATA_DIR)))
  
  Pheno = readTable(sprintf("%s/BREAST/MICROARRAY_GPL96/TUMOR_AGGREGATED/TUMOR_Pheno.csv.gz", DATA_DIR))
  basePheno = readTable(sprintf("%s/BREAST/MICROARRAY_GPL96/NORMAL_AGGREGATED/NORMAL_Pheno.csv.gz", DATA_DIR))
  
  common = intersect(rownames(baseMat), rownames(Mat))
  
  Mat = Mat[common, ] #11764
  baseMat = baseMat[common, ]
  
  sel = sample(1:ncol(baseMat), ceiling(ncol(baseMat) * 0.5))
  
  Groups = factor(c(rep("TUMOR", ncol(Mat)), rep("NORMAL", ncol(baseMat[, -sel]))))
  Mat = cbind(Mat, baseMat[, -sel])
  baseMat = baseMat[, sel]
  

  # ================ compute divergence ================ 
  
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat, 
                                      gamma = 1:9/10,
                                      Groups=Groups, classes=c("NORMAL", "TUMOR"))
  
  df_GPL96 = data.frame(N=div$div$count.div, Groups=Groups)
  
  rm(Mat, baseMat, div, Pheno, Groups)
  
  # ====================================================
  # GPL1708
  # ====================================================
  
  cat("Loading gpl1708 data..\n")
  
  Mat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL1708/TUMOR_AGGREGATED/TUMOR_Expression.csv.gz", DATA_DIR)))
  baseMat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL1708/NORMAL_AGGREGATED/NORMAL_Expression.csv.gz", DATA_DIR)))
  
  Pheno = readTable(sprintf("%s/BREAST/MICROARRAY_GPL1708/TUMOR_AGGREGATED/TUMOR_Pheno.csv.gz", DATA_DIR))
  basePheno = readTable(sprintf("%s/BREAST/MICROARRAY_GPL1708/NORMAL_AGGREGATED/NORMAL_Pheno.csv.gz", DATA_DIR))
  
  # remove sample with many missing values
  Mat = Mat[, -which(colnames(Mat) == "GSM508331")]
  
  Mat = Mat[-which(apply(Mat, 1, function(x) ! all(is.finite(x)))), ]
  
  common = intersect(rownames(baseMat), rownames(Mat))
  
  Mat = Mat[common, ] #18818
  baseMat = baseMat[common, ]
  
  sel = sample(1:ncol(baseMat), ceiling(ncol(baseMat) * 0.5))
  
  Groups = factor(c(rep("TUMOR", ncol(Mat)), rep("NORMAL", ncol(baseMat[, -sel]))))
  Mat = cbind(Mat, baseMat[, -sel])
  baseMat = baseMat[, sel]
  
  # ================ compute divergence ================ 
  
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat, 
                                      gamma = 1:9/10,
                                      Groups=Groups, classes=c("NORMAL", "TUMOR"))
  
  df_GPL1708 = data.frame(N=div$div$count.div, Groups=Groups)
  
  rm(Mat, baseMat, div, Pheno, Groups)
  
  # ====================================================
  # Microarray
  # ====================================================
  
  dfList = list(
    df_GPL570,
    df_GPL96,
    df_GPL1708
  )
  names(dfList) = c("GPL570", "GPL96", "GPL1708")
  
  # ====== save ======
  df1_3 = make_facet_df(dfList, prefix="Breast")
  
  df1_3$TissuePrefixed = factor(df1_3$TissuePrefixed,
                                levels=c( "Breast\nGPL570", "Breast\nGPL96", "Breast\nGPL1078"),
                                labels=c( "Breast\nGPL570", "Breast\nGPL96", "Breast\nGPL1708"))
  
  df1_3$Tissue = factor(df1_3$Tissue,
                        levels=c("GPL570", "GPL96", "GPL1078"),
                        labels=c("GPL570", "GPL96", "GPL1708"))
  
  df1_3$TissueSuffixed = factor(df1_3$TissueSuffixed,
                                levels=c("GPL570\n", "GPL96\n", "GPL1078\n"),
                                labels=c("GPL570\n", "GPL96\n", "GPL1708\n"))
  
  save(df1_3, file="obj/1_3.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()  
  
}, error = function(e){ print(e) })

sink()










