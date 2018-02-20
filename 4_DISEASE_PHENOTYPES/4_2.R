
sink("log/4_2.log.txt", split=TRUE)

tryCatch({
  
  set.seed(1)
  
  library(plyr)
  library(rutils)
  library(divergence)
  
  source("../src/util.R")
  source("util_4.R")
  
  source("vars.R") # load DATA_DIR

  dfList = list()
  dfMeanList = list()
  
  # ====================================================
  # PROSTATE, RNA-SEQ, PRIMARY GLEASON
  # ====================================================
  
  cat("Loading prostate data..\n")
  
  dataMat = data.matrix(readTable(sprintf("%s/PROSTATE/RNASeq/TCGA/TCGA_TPM.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/PROSTATE/RNASeq/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
  
  baseMat = dataMat[, which(dataPheno$sample_type %in% c("Solid Tissue Normal"))]
  
  sel = which(dataPheno$sample_type %in% c("Primary Tumor") & 
                dataPheno$primary_pattern %in% c(3, 4, 5))
  
  Mat = dataMat[, sel]
  Groups = factor(dataPheno$primary_pattern[sel], levels=3:5)
  
  rm(dataMat, dataPheno)
  
  # ================ compute divergence ================
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat,
                                      gamma = 1:9/10)  
  
  dfList[["GLEASON"]] = data.frame(N=div$div$count.div,
                                   Groups=Groups,
                                   GroupsN=utils.make_n_factor(Groups))
  
  
  dfMeanList[["GLEASON"]] = get_mean_df(dfList$GLEASON)
  
  cat("Wilcoxon P:\n")
  print(compute_p_mat(dfList$GLEASON, "Groups", "N"))

  rm(Mat, baseMat, div, Groups)
  gc()
  
  # ====================================================
  # LUNG, METHYLATION, SMOKING
  # ====================================================

  cat("Loading lung data..\n")
  
  dataMat = data.matrix(readTable(sprintf("%s/LUNG/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Beta.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/LUNG/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Pheno.csv.gz", DATA_DIR))
  
  if(! all(is.finite(dataMat))){
    cat("Removing rows with missing values..\n")
    dataMat = dataMat[-which(apply(dataMat, 1, function(x) ! all(is.finite(x)))), ]
  }
  
  baseMat = dataMat[, which(dataPheno$sample_type %in% c("Solid Tissue Normal"))]
  
  sel_levels = c(
    "Lifelong Non-smoker",
    "Current reformed smoker for > 15 years",
    "Current reformed smoker for < or = 15 years",
    "Current smoker"
  )
  new_levels = c(
    "Non-smoker",
    "Reformed-smoker(>15 years)",
    "Reformed-smoker(<15 years)",
    "Current smoker"
  )
  
  sel = which(dataPheno$sample_type %in% c("Primary Tumor") & 
                dataPheno$tobacco_smoking_history %in% sel_levels)
  
  Mat = dataMat[, sel]
  Groups = factor(dataPheno$tobacco_smoking_history[sel], 
                  levels=sel_levels, labels=new_levels)
  
  rm(dataMat, dataPheno)
  
  # ================ compute divergence ================ 
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat,
                                      computeQuantiles = FALSE,
                                      gamma = 1:9/10)    
  
  dfList[["SMOKING"]] = data.frame(N=div$div$count.div,
                                   Groups=Groups,
                                   GroupsN=utils.make_n_factor(Groups))
  
  dfMeanList[["SMOKING"]] = get_mean_df(dfList$SMOKING)
  
  cat("Wilcoxon P:\n")
  print(compute_p_mat(dfList$SMOKING, "Groups", "N"))
  
  rm(Mat, baseMat, div, Groups)
  gc()
  
  # ====================================================
  # BREAST, GPL96, HISTOLOGICAL GRADE
  # ====================================================

  cat("\nLoad breast data..\n")  
  
  Mat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL96/TUMOR_AGGREGATED/TUMOR_Expression.csv.gz", DATA_DIR)))
  baseMat = data.matrix(readTable(sprintf("%s/BREAST/MICROARRAY_GPL96/NORMAL_AGGREGATED/NORMAL_Expression.csv.gz", DATA_DIR)))
  
  Pheno = readTable(sprintf("%s/BREAST/MICROARRAY_GPL96/TUMOR_AGGREGATED/TUMOR_Pheno.csv.gz", DATA_DIR))
  rownames(Pheno) = Pheno$sample
  Pheno = Pheno[colnames(Mat), ]
  
  sel = which(Pheno$HIST_GRADE %in% 1:3)
  
  Mat = Mat[, sel]
  Pheno = Pheno[sel, ]

  Groups = factor(Pheno$HIST_GRADE, levels=1:3)
    

  # ================ compute divergence ================ 
  div = computeUnivariateDigitization(Mat=Mat, baseMat=baseMat,
                                      gamma = 1:9/10)    
  
  dfList[["HIST"]] = data.frame(N=div$N,
                                   Groups=Groups,
                                   GroupsN=utils.make_n_factor(Groups))
  
  
  dfMeanList[["HIST"]] = get_mean_df(dfList$HIST)
  
  cat("Wilcoxon P:\n")
  print(compute_p_mat(dfList$HIST, "Groups", "N"))
  
  rm(Mat, baseMat, div, Groups)
  gc()
  
  # ====== save ======
  save(dfList, dfMeanList, file="obj/4_2.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()  
  

}, error = function(e){ print(e) })

sink()






