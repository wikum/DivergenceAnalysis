

sink("log/6_2.log.txt", split=TRUE)

tryCatch({
  
  set.seed(1)
  
  library(plyr)
  library(divergence)
  
  source("../src/util.R")
  
  source("../vars.R") # load DATA_DIR
  
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

  FGAS <- findMultivariateGammaAndSupport(Mat=MatNormal[, -111], genesets = HallMarks,
                                    gamma = seq(0.02, 0.99,0.02),beta =0.95, alpha= 0.001,method="euclidean")

  centermatrix <- FGAS$MultiSetSupport$Centermatrix_list
  radius <- FGAS$MultiSetSupport$Radius_list
  
  Binary_norm <- computeMultivariateBinaryMatrix(Mat=MatNormal, genesets = HallMarks,
                                             Centermatrix_list = centermatrix , Radius_list = radius, method="euclidean" )
  
  Binary_Cancer <- computeMultivariateBinaryMatrix(Mat= MatCancer, genesets = HallMarks, 
                                               Centermatrix_list = centermatrix , Radius_list = radius, method="euclidean")
  
  
  #### select optimal gamma
  
  opt_gamma <- FGAS$optimal_gamma
  
  
  #############heatmap####
  
  classes <- matrix(c('LuminalA', 'LuminalB', 'Basal', 'HER2-E', 'Positive', 'Negative'),nrow=1)
  Mat_g <- Binary_Cancer
  
  ### calculate divergent probability for each single pathway
  
  P <- apply(classes,2, function(x){
    Mat_c <- Mat_g[,Pheno$group ==x]
    pa <- apply(Mat_c,1, function(y){
      sum(y==1)/length(y)
    })
    
  })
  
  
  
  PathwayName <- colnames(HallMarks)
  PathwayName <- gsub('*HALLMARK_','', PathwayName)
  
  rownames(P)<- PathwayName
  colnames(P)<- classes
  
  #### reorder pathway according to the variable of Divergent probability
  
  VARhallmark <- apply(P,1, function(x){
    var(x)
  })
  
  order1 <- sort(VARhallmark,decreasing = TRUE,index.return=TRUE)
  
  
  
  P_reorder <- as.matrix(P[order1$ix,])
  
  if(FALSE){
  library(gplots)
  pdf("heatmap.pdf")
  heatmap.2(P_reorder,col=redblue(100),density.info="none", trace="none",
            Rowv = FALSE, Colv=FALSE, cexRow = 0.5, cexCol = 0.8,
            main='Divergent Probability Heatmap(50%)')
  dev.off()
  }
  
  
  # ====== save ======
  save(P_reorder, file="obj/6_2.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()   
  
}, error = function(e){ print(e) })

sink()






