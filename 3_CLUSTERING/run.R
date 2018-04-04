

if(! dir.exists("obj"))
  dir.create("obj")

if(! dir.exists("log"))
  dir.create("log")

tryCatch({
  cat("\n 1 \n ===================== \n")
  source("3_1.R")
}, error = function(e){print(e)})

