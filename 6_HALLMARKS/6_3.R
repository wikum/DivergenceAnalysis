

sink("log/6_3.log.txt", split=TRUE)

tryCatch({
  
  set.seed(1)
  
  library(plyr)
  library(divergence)
  library(rpart)
  
  source("../src/util.R")
  
  source("../vars.R") # load DATA_DIR
  
  DATA_DIR = "/Volumes/WORK/MechPred/DATA"
  
  HallMarks <- data.matrix(read.csv(sprintf("%s/GENESETS/HALLMARKS.txt", DATA_DIR), 
                                    head=TRUE, sep=" ",na.strings="NA",stringsAsFactors = FALSE))
  
  cat("Loading breast data\n")
  dataMat = data.matrix(readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_TPM.csv.gz", DATA_DIR)))
  dataPheno = readTable(sprintf("%s/BREAST/RNASeq/TCGA/TCGA_Pheno.csv.gz", DATA_DIR))
  
  cat("Making groups\n")
  sel_n = which(dataPheno$sample_type == "Solid Tissue Normal")
  baseMat = dataMat[, sel_n]
  
  sel_t = list(
    which(dataPheno$PAM50_mRNA_nature2012 == "Luminal A"),
    which(dataPheno$PAM50_mRNA_nature2012 == "Luminal B"),
    which(dataPheno$PAM50_mRNA_nature2012 == "HER2-enriched"),
    which(dataPheno$PAM50_mRNA_nature2012 == "Basal-like"),
    
    which(dataPheno$ER_Status_nature2012 == "Positive"),
    which(dataPheno$ER_Status_nature2012 == "Negative")
  )
  
  sel_t_ids = unlist(sel_t)
  
  new_labels = c('LuminalA', 'LuminalB', 'HER2-E', 'Basal', 'Positive', 'Negative')
  
  groups_v = unlist(lapply(1:length(sel_t), function(j) rep(new_labels[j], length(sel_t[[j]]))))
  
  Mat = dataMat[, sel_t_ids]
  
  Groups = factor(groups_v, levels=new_labels)
  
  Pheno = data.frame(sample=colnames(Mat), group=Groups)
  
  cat("Computing quantiles\n")
  MatNormal = computeQuantileMatrix(baseMat)
  MatCancer = computeQuantileMatrix(Mat)
  
  cat("Running experiment..\n")
  
  # select only luminal A and luminal B samples
  sel_lum = which(Pheno$group %in% c('LuminalA', 'LuminalB'))
  MatLuminal = MatCancer[, sel_lum]
  PhenoLum = Pheno[sel_lum, ]
  rownames(PhenoLum) <- PhenoLum$sample
  PhenoLum$group = factor(as.character(PhenoLum$group), levels=c('LuminalA', 'LuminalB'))
  
  table(PhenoLum$group)
  # LuminalA LuminalB 
  # 231      127 
  
  cat("Running experiment..\n")
  
  # compute decision tree
  
  
  ############## compute decision tree  
  
  
  
  # find optimal gamma and corresponding support for 50 HallMarks Pathways with given beta and alpha
  MatNormal <- MatNormal[,c(1:110,112)]
  
  FGAS <- findMultivariateGammaAndSupport(Mat=MatNormal, genesets = HallMarks,
                                          gamma = seq(0.02, 0.99,0.02),beta =0.95, alpha= 0.001,method="euclidean")
  
  
  # select optimal gamma
  
  opt_gamma <- FGAS$optimal_gamma
  
  
  # transform all cancer data into binary
  
  centermatrix <- FGAS$MultiSetSupport$Centermatrix_list
  radius <- FGAS$MultiSetSupport$Radius_list
  
  Binary_Cancer <- computeMultivariateBinaryMatrix(Mat= MatCancer, genesets = HallMarks, 
                                               Centermatrix_list = centermatrix , Radius_list = radius, method="euclidean")
  
  # select LuminalA and  LuminalB Binary Expression 
  
  BinaryLuminal <- Binary_Cancer[, sel_lum]

  
  # compute decision tree with 100 random split
  
  BinaryLuminal <- t(BinaryLuminal)
  BinaryLuminal[BinaryLuminal==1] <- "divergent"
  BinaryLuminal[BinaryLuminal==0] <- "nondivergent"
  classes <- c('LuminalB', 'LuminalA')
  
  
  
  BinaryLuminal<- as.data.frame(unclass(BinaryLuminal))
  rownames(BinaryLuminal)<- PhenoLum$sample
  
  TEST_ACCURACY<- c()
  
  
  for (j in 1:100){
    
    
    A<- as.matrix(PhenoLum$sample[PhenoLum$group==classes[1]])
    B<- as.matrix(PhenoLum$sample[PhenoLum$group==classes[2]])
    
    
    # sampling the train data 
    train_A<- sample(A, size=floor(1/2*length(A)))
    train_B<- sample(B, size = floor(1/2*length(B)))
    # balance sample size
    train_AC <-rep(train_A, 2)
    train_BC <- rep(train_B,1)
    
    
    train_sample<-c(train_AC, train_BC)
    
    train_Pheno<- PhenoLum[train_sample,]
    
    traindata <- BinaryLuminal[train_sample,]
    
    test_A<- setdiff(A, train_A)
    test_B<- setdiff(B,train_B)
    test_AC <-rep(test_A, 2)
    test_BC <- rep(test_B,1)
    
    test_sample<- c(test_AC,test_BC)
    test_Pheno<- PhenoLum[test_sample,]
    
    testdata <- BinaryLuminal[test_sample,]
    
    
    
    
    #### train the decision tree on training set
    
    trainDT <- rpart(train_Pheno$group ~ ., as.data.frame(traindata))
    
    
    train_pred_label <- predict(trainDT, as.data.frame(traindata), type='class')
    
    #### test the tree on test set 
    test_pred_label<-predict(trainDT,as.data.frame(testdata), type='class')
    
    train_accuracy<-getPredictionStats(train_pred_label,train_Pheno$group,classes )
    test_accuracy<-getPredictionStats(test_pred_label,test_Pheno$group,classes)
    
    
    TEST_ACCURACY<- rbind(TEST_ACCURACY, c(train_accuracy,test_accuracy))
    
    
  }
  
  # print out average of 100 iteration
  
  e<-apply(TEST_ACCURACY,2, mean)
  
  # ====== save ======
  save(trainDT, file="obj/6_3.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()     
  
}, error = function(e){ print(e) })

sink()








