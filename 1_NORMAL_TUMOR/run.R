

if(FALSE){
  library(devtools)
  current_path = getwd()
  setwd("../../../devel_divergence/divergenceDevel")
  install("divergence")
  reload("divergence")
  setwd(current_path)
}

if(! dir.exists("obj"))
  dir.create("obj")

if(! dir.exists("log"))
  dir.create("log")

tryCatch({
  cat("\n RNA-Seq NORMAL vs. TUMOR \n ===================== \n")
  source("1_1.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n Methylation NORMAL vs. TUMOR \n ===================== \n")
  source("1_2.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n Microarray NORMAL vs. TUMOR \n ===================== \n")
  source("1_3.R")
}, error = function(e){print(e)})


