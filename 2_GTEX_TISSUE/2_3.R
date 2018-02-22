
sink("log/2_3.log.txt", split=TRUE)

tryCatch({

  source("../src/util.R")
  
  source("../vars.R")
  
  lapply_c = function(x, ...){
    
    result = lapply(x, ...)
    names(result) = x
    
    result
    
  }
  
  # ====================================================
  # GTEX
  # ====================================================
  
  if(! dir.exists("2_3"))
    dir.create("2_3")
  
  cat("Loading data..\n")
  Mat = data.matrix(readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Counts_normalized.csv", DATA_DIR)))
  
  Pheno = readTable(sprintf("%s/MULTI-TISSUE/RNASeq/GTEX/GTEX_Pheno.csv", DATA_DIR))
  rownames(Pheno) = Pheno$sample
  # all(rownames(Pheno) == colnames(Mat)) #TRUE
  
  n = 50
  # select n samples from each subtype for training
  
  u = table(Pheno$tissue_subtype) # 53 subtypes
  v = u[u >= n] # 48 subtypes
  sel_types = names(v)
  
  sel = which(Pheno$tissue_subtype %in% sel_types)
  
  Mat = Mat[, sel]
  Pheno = Pheno[sel, ]
  
  print(table(is.finite(Mat)))
  
  if(! all(is.finite(Mat))){
    cat("Removing rows with missing values..\n")
    Mat = Mat[-which(apply(Mat, 1, function(x) ! all(is.finite(x)))), ]
  }
  
  MatP = computeQuantileMatrix(Mat)

  rm(Mat)
  gc()
  
  Groups = factor(as.character(Pheno$tissue_subtype))
  
  trainIds = lapply_c(levels(Groups), function(t)
                    sample(which(Groups == t), n))

  sel_train = unlist(trainIds)
  
  tissueGroups = Groups[-sel_train]
  
  n_map = make_n_factor_map(tissueGroups)
  rownames(n_map) = n_map$levels
  n_map = n_map
  
  save(Groups, trainIds, tissueGroups, n_map, file="2_3/tissueGroups.rda")
  
  divMat = c()
  
  for(t in levels(Groups)){
    
    cat("\n(", t, ")\n")
    
    t_dir = sprintf("2_3/%s", t)
    
    if(! dir.exists(t_dir))
      dir.create(t_dir)
    
    # ================ compute divergence ================ 

    div = computeUnivariateDigitization(Mat=MatP[, -sel_train], baseMat=MatP[, trainIds[[t]] ],
                                        gamma = 1:9/10)
    
    divMat = cbind(divMat, div$div$count.div)
    
  }
  
  colnames(divMat) = levels(Groups)
  
  predictions = levels(Groups)[apply(divMat, 1, function(x) which.min(x)[1])]
  
  accMat = matrix(0, length(levels(Groups)), length(levels(Groups)))
  rownames(accMat) = levels(Groups)
  colnames(accMat) = levels(Groups)
  
  accMat = t(sapply(levels(Groups), function(t){
    x = which(tissueGroups == t)
    sapply(levels(Groups), function(p) sum(predictions[x] == p))/length(predictions[x])
  }))
  
  rownames(accMat) = as.character(n_map[rownames(accMat), "levelsSuffixed"])

  save(divMat, accMat,  file="2_3/results.rda")
  
  meltedAccMat = meltMat(accMat)
  
  # ====== save ======
  save(meltedAccMat, file="obj/2_3.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()    
  
}, error = function(e){ print(e) })

sink()






