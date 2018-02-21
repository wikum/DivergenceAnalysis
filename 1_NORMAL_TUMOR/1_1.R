
sink("log/1_1.log.txt", split=TRUE)

tryCatch({

  library(plyr)
  library(divergence)
  
  source("../src/util.R")

  source("util_1.R")
  source("../vars.R")
  
  # ====================================================
  # TCGA RNA-Seq  
  # ====================================================
  
  Mat = data.matrix(readTable(sprintf("%s/MULTI-TISSUE/RNASeq/TCGA/TCGA_Counts.csv.gz", DATA_DIR)))
  Pheno = readTable(sprintf("%s/MULTI-TISSUE/RNASeq/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
  
  sel_tissues = c("Breast", "Colon", "Kidney", "Liver", "Lung", "Prostate", "Thyroid Gland")
  
  sel_ids = which(Pheno$anatomical_origin %in% sel_tissues)
  Mat = Mat[, sel_ids]
  Pheno = Pheno[sel_ids, ]
  
  tissues = levels(Pheno$anatomical_origin)
  
  MatQ = computeQuantileMatrix(Mat)

  #sel_types = c("Solid Tissue Normal", "Primary Tumor")
  #all_ids = which(Pheno$anatomical_origin %in% sel_tissues & 
  #                  (Pheno$sample_type == "Solid Tissue Normal" | 
  #                     Pheno$sample_type == "Primary Tumor"))
  
  idList = list()
  resultList = list()
    
  for(t in sel_tissues){
      
      cat("(", t, ")\n")
      
      normal_ids = which(Pheno$anatomical_origin == t & Pheno$sample_type == "Solid Tissue Normal")
      tumor_ids = which(Pheno$anatomical_origin == t & Pheno$sample_type == "Primary Tumor")
      
      sel_train = sample(normal_ids, ceiling(length(normal_ids) * 0.5))
      sel_test = c(setdiff(normal_ids, sel_train), tumor_ids)
      
      Groups = factor(c(
        rep("NORMAL", length(setdiff(normal_ids, sel_train))), 
        rep("TUMOR", length(tumor_ids))), 
        levels=c("NORMAL", "TUMOR"))
      
      # ================ compute divergence ================ 
  
      #rr = findUnivariateGammaWithSupport(Mat=MatQ[, sel_train])

      div = computeUnivariateDigitization(Mat = MatQ[, sel_test], 
                                          baseMat = MatQ[, sel_train], 
                                          computeQuantiles = FALSE,
                                          Groups = Groups,
                                          classes = c("NORMAL", "TUMOR")
      )

      idList[[t]] = list(normal=normal_ids, tumor=tumor_ids, train=sel_train, test=sel_test)
      resultList[[t]] = list(div=div$div, features=div$features.div, Groups=Groups)
      
    }
    
  save(idList, resultList, file="obj/TCGA_RNASEQ_results.rda")
    
  dfList = lapply(1:length(resultList), function(j)
    data.frame(
      N=resultList[[j]]$div[, "count.div"], 
      Groups=resultList[[j]]$Groups
    )
  )
  names(dfList) = names(resultList)
  
  # ====== save ======
  df1_1 = make_facet_df(dfList, suffix="RNASeq")
  save(df1_1, file="obj/1_1.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()
  
}, error = function(e){ print(e) })

sink()










