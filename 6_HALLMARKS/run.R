

if(! dir.exists("obj"))
  dir.create("obj")

if(! dir.exists("log"))
  dir.create("log")

tryCatch({
  cat("\n 1 \n ===================== \n")
  source("6_1.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n 2 \n ===================== \n")
  source("6_2.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n 3 \n ===================== \n")
  source("6_3.R")
}, error = function(e){print(e)})

