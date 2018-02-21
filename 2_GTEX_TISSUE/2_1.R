
sink("log/2_1.log.txt", split=TRUE)

tryCatch({

  set.seed(1)
  
  library(plyr)
  library(divergence)
  library(rutils)
  
  source("../src/util.R")
  
  source("../vars.R")
  
  # ====================================================
  # GTEX
  # ====================================================
  
  if(! dir.exists("2_1"))
    dir.create("2_1")
  
  cat("Loading data..\n")
  
    Mat = data.matrix(readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Counts_normalized.csv", DATA_DIR)))
  
    Pheno = readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Pheno.csv", DATA_DIR))
    rownames(Pheno) = Pheno$sample
    # all(rownames(Pheno) == colnames(Mat)) #TRUE
  
    sel_types = c("Breast - Mammary Tissue",
                "Adipose - Subcutaneous",
                "Adipose - Visceral (Omentum)",
                "Skin - Sun Exposed (Lower leg)",
                "Skin - Not Sun Exposed (Suprapubic)",
                "Brain - Cortex",
                "Brain - Cerebellum")
  
    new_types = c("Breast-Mammary Tissue",
      "Adipose-Subcutaneous",
      "Adipose-Visceral",
      "Skin-Sun Exposed",
      "Skin-Not Sun Exposed",
      "Brain-Cortex",
      "Brain-Cerebellum")
  
    sel = which(Pheno$tissue_subtype %in% sel_types)
  
    Mat = Mat[, sel]
    Pheno = Pheno[sel, ]
    
    Groups = factor(factor(as.character(Pheno$tissue_subtype)),
                    levels=sel_types,
                    labels=new_types)
    
  train_types = c(
    "Breast-Mammary Tissue",
    "Skin-Sun Exposed"
  )
  
  dfList = list()
    
  for(t in train_types){
    
    cat("\n(", t, ")\n")
    
    t_dir = sprintf("2_1/%s", t)
    
    if(! dir.exists(t_dir))
      dir.create(t_dir)
    
    sel_train = sample(which(Groups == t), ceiling(sum(Groups == t) * 0.5))
    
    tissueGroups = Groups[-sel_train]
    
    # ================ compute divergence ================ 
    cat("Computing divergence..\n")
    div = computeUnivariateDigitization(Mat=Mat[, -sel_train], baseMat=Mat[, sel_train],
                                        gamma = 1:9/10)
    
    save(div, tissueGroups, file=sprintf("%s/div.rda", t_dir))
    
    dfList[[t]] = data.frame(N=div$div$count.div,
                             Groups=tissueGroups,
                             GroupsN=utils.make_n_factor(tissueGroups),
                             trainGroup=factor(tissueGroups == t))
    
  }
  
  # ====== save ======
  save(dfList, file="obj/2_1.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()  
  
}, error = function(e){ print(e) })

sink()






