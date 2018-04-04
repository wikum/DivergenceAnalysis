
sink("log/3_1.log.txt", split=TRUE)

tryCatch({

  set.seed(1)

  library(plyr)
  library(divergence)
  library(ClusterR)
  library(kernlab)
  
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
    which(dataPheno$PAM50_mRNA_nature2012 == "Basal-like")
  )

  sel_t_ids = unlist(sel_t)

  new_labels = c('LuminalA', 'LuminalB', 'HER2-E', 'Basal')

  groups_v = unlist(lapply(1:length(sel_t), function(j) rep(new_labels[j], length(sel_t[[j]]))))

  Mat = dataMat[, sel_t_ids]

  Groups = factor(groups_v, levels=new_labels)

  Pheno = data.frame(sample=colnames(Mat), group=Groups)
  Pheno$basal = factor(ifelse(Pheno$group == "Basal", "Basal-like", "Non-basal-like"))

  cat("Computing quantiles\n")
  MatNormal = computeQuantileMatrix(baseMat)
  TCGA.Quantile = computeQuantileMatrix(Mat)

  cat("Running experiment..\n")

  #rownames(dataPheno) = dataPheno$sample
  rownames(Pheno) = Pheno$sample
  TCGA.Pheno = Pheno[colnames(TCGA.Quantile), ]

  normal_intervals = computeUnivariateSupport(MatNormal)
  NormalInterval.min = normal_intervals$Ranges$baseline.low
  NormalInterval.max = normal_intervals$Ranges$baseline.high

  cat("Digitization..\n")
  TCGA.Ternary = computeUnivariateTernaryMatrix(Mat=TCGA.Quantile, Ranges=normal_intervals$Ranges)

  ###### PAM50 Samples
  # Determine number of clusters
  I.PAM50=c(which(TCGA.Pheno$group=="Basal"),which(TCGA.Pheno$group=="HER2-E"),
          which(TCGA.Pheno$group=="LuminalA"),which(TCGA.Pheno$group=="LuminalB"))
  TCGA.Quantile.PAM50=TCGA.Quantile[,I.PAM50]
  TCGA.Ternary.PAM50=TCGA.Ternary[,I.PAM50]
  wss.Q.PAM50 <- replicate(10,NA)
  wss.T.PAM50 <- replicate(10,NA)
  bss_tss.Q.PAM50 <- replicate(10,NA)
  bss_tss.T.PAM50 <- replicate(10,NA)
  for (i in 1:10){
    wss.Q.PAM50[i] <- sum(kmeans(t(TCGA.Quantile.PAM50), centers=i)$withinss)
    wss.T.PAM50[i] <- sum(kmeans(t(TCGA.Ternary.PAM50), centers=i)$withinss)
    bss_tss.Q.PAM50[i] <- kmeans(t(TCGA.Quantile.PAM50), centers=i)$betweenss/kmeans(t(TCGA.Quantile.PAM50), centers=i)$totss
    bss_tss.T.PAM50[i] <- kmeans(t(TCGA.Ternary.PAM50), centers=i)$betweenss/kmeans(t(TCGA.Ternary.PAM50), centers=i)$totss
  }

  cat("Running ks test..\n")
  I.train.pam50 = unlist(lapply(levels(TCGA.Pheno$group), function(x) 
    sample(which(TCGA.Pheno$group == x), ceiling(sum(TCGA.Pheno$group == x) * 0.5)) ))
  I.test.pam50=setdiff(c(1:length(TCGA.Pheno$group)),I.train.pam50)

  kw = t(sapply(rownames(TCGA.Quantile.PAM50), function(x){
    u = kruskal.test(x=TCGA.Quantile.PAM50[x, I.train.pam50], g=TCGA.Pheno$group[I.train.pam50])
    c(stat=u$statistic, pval=u$p.value)
  }))
  kw = data.frame(kw, pbon=p.adjust(kw[, "pval"], method="bon"))
  kw = kw[order(kw$pbon, decreasing = FALSE, method = "radix"), ]

  I.KWtest.PAM50.train = rownames(kw)[1:100]

  cat("PCA..\n")
  dat = t(TCGA.Quantile.PAM50[I.KWtest.PAM50.train,I.train.pam50])
  pca_dat = stats::princomp(dat)$scores[, 1:2]
  km = KMeans_rcpp(pca_dat, clusters = 2, num_init = 5, max_iters = 100)
  spec=specc(pca_dat,centers=2)
  
  df_quantile = data.frame(pca_dat, spec, TCGA.Pheno[rownames(pca_dat), ])
  colnames(df_quantile)[1:3] = c("PC1", "PC2", "cluster")

  #c = computeChiSquaredTest(Mat=TCGA.Ternary.PAM50[, I.train.pam50],
  #                    Groups=factor(TCGA.Pheno$group[I.train.pam50] == "Basal"),
  #                    classes=c("FALSE", "TRUE"))

  chisq = function(x, y, z=c(-1, 0, 1)){
    c = sapply(unique(y), function(a) sapply(z, function(b) sum(x[y == a] == b)))
    colnames(c) = unique(y)
    rownames(c) = z
    cx = chisq.test(c)
    c(statistic=cx$statistic, pval=cx$p.value)
  }
  C = t(apply(TCGA.Ternary.PAM50[, I.train.pam50], 1, 
            function(x) chisq(x, TCGA.Pheno$group[I.train.pam50], z=c(-1, 0, 1)))
  )
  C = C[order(C[, "pval"]), ]
  
  chisq.pvalue.PAM50.train = rownames(C)[1:100]
  #chisq.pvalue.PAM50.train = sample(rownames(C), 100)

  dat = t(TCGA.Ternary.PAM50[chisq.pvalue.PAM50.train, I.train.pam50])
  pca_dat = stats::princomp(dat)$scores[, 1:2]
  km = KMeans_rcpp(pca_dat, clusters = 2, num_init = 5, max_iters = 100)
  spec=specc(pca_dat,centers=2)
  #table(spec)

  df_ternary = data.frame(pca_dat, spec, Pheno[rownames(pca_dat), ])
  colnames(df_ternary)[1:3] = c("PC1", "PC2", "cluster")

  # ====== save ======
  save(df_quantile, df_ternary, file="obj/3_1.rda")
  # ==================
  
  sessionInfo()
  
  rm(list=ls())
  gc()   
  
}, error = function(e){ print(e) })

sink()



