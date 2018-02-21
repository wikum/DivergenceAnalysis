
sink("log/4_3.log.txt", split=TRUE)

tryCatch({

  set.seed(1)
  
  library(plyr)
  library(rutils)
  library(divergence)
  
  source("../src/util.R")
  source("util_4.R")
  
  source("vars.R") # load DATA_DIR


  # ====================================================
  # ADENOMA CARCINOMA
  # ====================================================

  cat("\nADENOMA CARCINOMA\n")
  
  if(FALSE){
  
  Mat1 = readTable("/Volumes/WORK/Dropbox (MechPred)/Additional Data/MAT/matGSE4183.csv.gz")
  Pheno1 = readTable("/Volumes/WORK/Dropbox (MechPred)/Additional Data/PHENO/phenoGSE4183.csv.gz")
  rownames(Mat1) = as.character(Mat1[, 1])
  Mat1 = log2(data.matrix(Mat1[, -1])+1)
  
  all(Pheno1$SAMPLE == colnames(Mat1))
  
  Mat2 = readTable("/Volumes/WORK/Dropbox (MechPred)/Additional Data/MAT/matGSE20916.csv.gz")
  Pheno2 = readTable("/Volumes/WORK/Dropbox (MechPred)/Additional Data/PHENO/phenoGSE20916.csv.gz")
  rownames(Mat2) = as.character(Mat2[, 1])
  Mat2 = data.matrix(Mat2[, -1])
  
  all(Pheno2$SAMPLE == colnames(Mat2))
  
  common = intersect(rownames(Mat1), rownames(Mat2))
  
  Mat1 = Mat1[common, ]
  Mat2 = Mat2[common, ]
  
  Mat1P = getPercentileMat(Mat1)
  Mat2P = getPercentileMat(Mat2)
  
  rm(Mat1)
  rm(Mat2)
  
  refMat1P = Mat1P[, which(Pheno1$TYPE == "NORMAL")]
  sel1 = which(Pheno1$TYPE %in% c("ADENOMA", "CRC"))
  Mat1P = Mat1P[, sel1]
  Pheno1 = Pheno1[sel1, ]
  
  Groups1 = factor(Pheno1$TYPE, levels=c("ADENOMA", "CRC"),
                   labels=c("ADENOMA", "CARCINOMA"))
  
  
  refMat2P = Mat2P[, which(Pheno2$TYPE == "NORMAL")]
  sel2 = which(Pheno2$TYPE %in% c("ADENOMA", "CARCINOMA"))
  Mat2P = Mat2P[, sel2]
  Pheno2 = Pheno2[sel2, ]
  
  Groups2 = factor(Pheno2$TYPE, levels=c("ADENOMA", "CARCINOMA"))
  
  
  M = cbind(Mat1P, refMat1P, Mat2P, refMat2P)
  G = factor(c(
    paste("1_", Groups1, sep=""),
    rep("1_NORMAL", ncol(refMat1P)),
    paste("2_", Groups2, sep=""),
    rep("2_NORMAL", ncol(refMat2P))
  ))
  
  
  save(M, G, file="4_3_data.rda")
  
  }else{
    load("4_3_data.rda")
  }
  
  # ================ compute divergence ================ 
  div0 = computeDivergences(Mat=M, 
                           refMat=M[, which(G %in% c("1_NORMAL", "2_NORMAL"))], 
                           upper="perc.high", lower="perc.low",
                           classes=G)
  
  div1 = computeDivergences(Mat=M, 
                            refMat=M[, which(G %in% c("1_NORMAL"))], 
                            upper="perc.high", lower="perc.low",
                            classes=G)
  
  div2 = computeDivergences(Mat=M, 
                            refMat=M[, which(G %in% c("2_NORMAL"))], 
                            upper="perc.high", lower="perc.low",
                            classes=G)
  
  save(div0, div1, div2, file="4_3_div.rda")
  
  source("../../../code/temp_util.R")
  
  
  sel = which(G %in% c("2_ADENOMA", "2_CARCINOMA"))
  
  df = data.frame(N=div2$N[sel],
                  Groups=factor(G[sel], levels=c("2_ADENOMA", "2_CARCINOMA"),
                                labels=c("Adenoma", "Carcinoma")),
                  GroupsN=make_n_factor(G[sel]))
  
  
  dfMean = get_mean_df(df)
  
  save(df, dfMean, file="4_3.rda")
  
  cat("Wilcoxon P:\n")
  print(compute_p_mat(df, "Groups", "N"))
  
  #boxplot(div$N ~ G)
  
  #rm(Mat1P, Mat2P, refMat1P, refMat2P, Pheno1, Pheno2, div)
  #gc()
  
  if(FALSE){
    
    load("4_3_div.rda")
    load("4_3_data.rda")
    
    #par(mfrow=c(1, 3))
    library(ggplot2)
    
    GG = paste(G, make_n_factor(G), sep="\n")
    
    pdf("adenoma_carcinoma.pdf")
    
    ggplot(data.frame(x=div0$N, y=GG), aes(x=y, y=x, fill=y))+geom_boxplot()+
      xlab("group")+ylab("size of divergent set")+guides(fill=FALSE)+ggtitle("trained on 1_NORMAL and 2_NORMAL")
    
    ggplot(data.frame(x=div1$N, y=GG), aes(x=y, y=x, fill=y))+geom_boxplot()+
      xlab("group")+ylab("size of divergent set")+guides(fill=FALSE)+ggtitle("trained on 1_NORMAL")

    ggplot(data.frame(x=div2$N, y=GG), aes(x=y, y=x, fill=y))+geom_boxplot()+
      xlab("group")+ylab("size of divergent set")+guides(fill=FALSE)+ggtitle("trained on 2_NORMAL")
    
    #boxplot(div0$N ~ G, main="trained on both", las=2)
    #boxplot(div1$N ~ G, main="trained on 1", las=2)
    #boxplot(div2$N ~ G, main="trained on 2", las=2)
    dev.off()
    
  }
  

}, error = function(e){ print(e) })

sink()






