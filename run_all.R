
if(! dir.exists("log"))
  dir.create("log")

sink("log/run_all.log.txt", split=TRUE)

# install necessary packages
source("install.R")

# run tumor vs normal univariate experiments
cat("Running tumor vs. normal experiments..\n")
cat("--------------------------------------\n")
  setwd("1_NORMAL_TUMOR/")
  source("run.R")
  setwd("../")

cat("Running GTEx tissue divergence experiments..\n")
cat("--------------------------------------\n")
  setwd("2_GTEX_TISSUE/")
  source("run.R")
  setwd("../")

cat("Running clustering experiments..\n")
cat("--------------------------------------\n")
  setwd("3_CLUSTERING/")
  source("run.R")
  setwd("../")
  
cat("Running disease phenotype comparison experiments..\n")
cat("--------------------------------------\n")
  setwd("4_DISEASE_PHENOTYPES/")
  source("run.R")
  setwd("../")

cat("Running feature selection experiments..\n")
cat("--------------------------------------\n")
  setwd("5_FEATURE_SELECTION/")
  source("run.R")
  setwd("../")

cat("Running hallmark geneset experiments..\n")
cat("--------------------------------------\n")
  setwd("6_HALLMARKS/")
  source("run.R")
  setwd("../")
  
cat("Making figures..\n")  
cat("----------------\n")
  source("figures.R")

sink()
