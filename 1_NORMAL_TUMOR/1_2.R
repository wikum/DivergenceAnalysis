
sink("log/1_2.log.txt", split=TRUE)

tryCatch({

  library(plyr)
  library(knitr)
  library(divergence)
  
  source("../src/util.R")
  
  source("util_1.R")
  source("../vars.R")
  
  # ====================================================
  # TCGA Methylation 450k - LUAD
  # ====================================================
  
  Mat = data.matrix(readTable(sprintf("%s/LUNG/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Beta.csv.gz", DATA_DIR)))
  Pheno = readTable(sprintf("%s/LUNG/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Pheno.csv.gz", DATA_DIR))
    
  Mat = Mat[-which(apply(Mat, 1, function(x) ! all(is.finite(x)))), ] #370864
    
  keep = which(Pheno$sample_type %in% c("Solid Tissue Normal", "Primary Tumor"))
    
  Mat = Mat[, keep]
  Pheno = Pheno[keep, ]
    
  Groups = factor(factor(as.character(Pheno$sample_type)),
                    levels=c("Solid Tissue Normal", "Primary Tumor"),
                    labels=c("NORMAL", "TUMOR"))
    
  sel_train = sample(which(Groups == "NORMAL"), ceiling(sum(Groups == "NORMAL") * 0.5))
    
  baseMat = Mat[, sel_train]
  Mat = Mat[, -sel_train]
  Pheno = Pheno[-sel_train, ]
  Groups = Groups[-sel_train]
    
  #save(Mat, baseMat, file="TCGA_Methyl_LUAD_data.rda")
  #save(Pheno, Groups, file="TCGA_Methyl_LUAD_Pheno.rda")
    
  # ================ compute divergence ================ 

  div = computeUnivariateDigitization(Mat = Mat, baseMat = baseMat,
                                      computeQuantiles = FALSE, 
                                      Groups = Groups, classes = c("NORMAL", "TUMOR"))
  
  save(div, Groups, file="obj/TCGA_Methyl_LUAD_div.rda")
  
  df_LUAD = list(N=div$div$count.div, Groups=Groups)
  
  rm(Pheno, div)
  
  # ====================================================
  # TCGA Methylation 450k - HNSC
  # ====================================================
  
  Mat = data.matrix(readTable(sprintf("%s/HEADNECK/METHYLATION_450k/TCGA/TCGA_Beta.csv.gz", DATA_DIR)))
  Pheno = readTable(sprintf("%s/HEADNECK/METHYLATION_450k/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
    
  Mat = Mat[-which(apply(Mat, 1, function(x) ! all(is.finite(x)))), ] #375044
    
  keep = which(Pheno$sample_type %in% c("Solid Tissue Normal", "Primary Tumor"))
    
  Mat = Mat[, keep]
  Pheno = Pheno[keep, ]
    
  #all(colnames(Mat) == Pheno$sample)
    
  Groups = factor(factor(as.character(Pheno$sample_type)),
                    levels=c("Solid Tissue Normal", "Primary Tumor"),
                    labels=c("NORMAL", "TUMOR"))
    
  # table(Groups, Pheno$sample_type)
    
  sel_train = sample(which(Groups == "NORMAL"), ceiling(sum(Groups == "NORMAL") * 0.5))
    
  baseMat = Mat[, sel_train]
  Mat = Mat[, -sel_train]
  Pheno = Pheno[-sel_train, ]
  Groups = Groups[-sel_train]
    
  #save(Mat, baseMat, file="TCGA_Methyl_HNSC_data.rda")
  #save(Pheno, Groups, file="TCGA_Methyl_HNSC_Pheno.rda")
    
  # ================ compute divergence ================ 

  div = computeUnivariateDigitization(Mat = Mat, baseMat = baseMat,
                                      computeQuantiles = FALSE, 
                                      Groups = Groups, classes = c("NORMAL", "TUMOR"))
  
  
  save(div, Groups, file="obj/TCGA_Methyl_HNSC_div.rda")
  
  df_HNSC = list(N=div$div$count.div, Groups=Groups)
  
  rm(Pheno, div)
  
  # ====================================================
  # TCGA Methylation 450k - PRAD
  # ====================================================
  
  Mat = data.matrix(readTable(sprintf("%s/PROSTATE/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Beta.csv.gz", DATA_DIR)))
  Pheno = readTable(sprintf("%s/PROSTATE/METHYLATION_450k/TCGA_Adenocarcinoma/TCGA_Pheno.csv.gz", DATA_DIR))
    
  Mat = Mat[-which(apply(Mat, 1, function(x) ! all(is.finite(x)))), ] #381182
    
  keep = which(Pheno$sample_type %in% c("Solid Tissue Normal", "Primary Tumor"))
    
  Mat = Mat[, keep]
  Pheno = Pheno[keep, ]
    
  # all(colnames(Mat) == Pheno$sample)
    
  Groups = factor(factor(as.character(Pheno$sample_type)),
                    levels=c("Solid Tissue Normal", "Primary Tumor"),
                    labels=c("NORMAL", "TUMOR"))
    
  # table(Groups, Pheno$sample_type)
    
  sel_train = sample(which(Groups == "NORMAL"), ceiling(sum(Groups == "NORMAL") * 0.5))
    
  baseMat = Mat[, sel_train]
  Mat = Mat[, -sel_train]
  Pheno = Pheno[-sel_train, ]
  Groups = Groups[-sel_train]
    
  #save(Mat, baseMat, file="TCGA_Methyl_PRAD_data.rda")
  #save(Pheno, Groups, file="TCGA_Methyl_PRAD_Pheno.rda")
    
  # ================ compute divergence ================ 

  div = computeUnivariateDigitization(Mat = Mat, baseMat = baseMat,
                                      computeQuantiles = FALSE, 
                                      Groups = Groups, classes = c("NORMAL", "TUMOR"))
  
  
  save(div, Groups, file="obj/TCGA_Methyl_PRAD_div.rda")
  
  df_PRAD = list(N=div$div$count.div, Groups=Groups)
  
  rm(Pheno, div)
  
  # ====================================================
  # TCGA Methylation 450k
  # ====================================================
  
  dfList = list(LUAD=df_LUAD,
                HNSC=df_HNSC,
                PRAD=df_PRAD)
  names(dfList) = c("Lung", "Head and Neck", "Prostate")
  
  # ====== save ======
  df1_2 = make_facet_df(dfList, suffix="Methylation 450k")
  
  save(df1_2, file="obj/1_2.rda")
  # ==================

  sessionInfo()
  
  rm(list=ls())
  gc()  

}, error = function(e){ print(e) })

sink()










