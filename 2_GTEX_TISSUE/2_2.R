
sink("log/2_2.log.txt", split=TRUE)

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
  
  if(! dir.exists("2_2"))
    dir.create("2_2")
  
  cat("Loading data..\n")
  
    Mat = data.matrix(readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Counts_normalized.csv", DATA_DIR)))
  
    Pheno = readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Pheno.csv", DATA_DIR))
    rownames(Pheno) = Pheno$sample
    # all(rownames(Pheno) == colnames(Mat)) #TRUE
  
    Groups = Pheno$tissue_subtype

    train_type = "Brain - Cortex"
    
    sel_train = sample(which(Groups == train_type), 
                       ceiling(sum(Groups == train_type) * 0.5))
    
    tissueGroups = Groups[-sel_train]
    
    # ================ compute divergence ================ 
    cat("Computing divergence..\n")
    div = computeUnivariateDigitization(Mat=Mat[, -sel_train], baseMat=Mat[, sel_train],
                                        gamma = 1:9/10)
    
    t_dir = sprintf("2_2/%s", train_type)
    
    if(! dir.exists(t_dir))
      dir.create(t_dir)
    
    save(div, tissueGroups, file=sprintf("%s/div.rda", t_dir))
    
    df = data.frame(N=div$div$count.div,
                    Groups=tissueGroups,
                    GroupsN=utils.make_n_factor(tissueGroups),
                    GroupsSuffixed=paste(tissueGroups, " (", utils.make_n_factor(tissueGroups), ")", sep=""),
                    trainGroup=factor(tissueGroups == train_type))
    
    # ====== save ======
    save(df, file="obj/2_2.rda")
    # ==================
    
    sessionInfo()
    
    rm(list=ls())
    gc()  
    
  
}, error = function(e){ print(e) })

sink()






