
library(tidyr)
library(plyr)
library(divergence)

### ========================================================
### Utility functions
### ========================================================

writeTable = function(obj, file, ...){
  write.table(obj, file, sep=",", row.names=TRUE, col.names=TRUE, ...)
}

readTable = function(file){
  if(length(grep("*.gz$", file)) > 0 && file.exists(file))
    read.table(gzfile(file), sep=",", header=TRUE, check.names=FALSE)
  else if(file.exists(file))
    read.table(file, sep=",", header=TRUE, check.names=FALSE)
  else if(file.exists(sprintf("%s.gz", file)))
    read.table(gzfile(sprintf("%s.gz", file)), sep=",", header=TRUE, check.names=FALSE)
  else if(file.exists(sprintf("%s.tar.gz", file)))
    read.table(gzfile(sprintf("%s.tar.gz", file)), sep=",", header=TRUE, check.names=FALSE)
  else
    stop(sprintf("ERROR: File not found either in regular or compressed form:%s\n", file))
}

meltMat = function(Mat){
  Mat = as.data.frame(Mat)
  meltedMat = gather(Mat)
  colnames(meltedMat) = c("colKey", "value")
  meltedMat$rowKey = rep(rownames(Mat), ncol(Mat))
  meltedMat
}


list_to_df = function(l){
  data.frame(x=unlist(l), y=rep(names(l), sapply(l, length)))
}


# e.g. x = c("A", "A", "A", "B", "B", "B", "C", "C", "D", "D", "D", "E", "E", "F")
make_n_factor = function(x){
  
  if(! is.factor(x))
    x = as.factor(x)
  
  y = levels=levels(x)
  z = sapply(y, function(u) sprintf("n=%d", sum(x == u, na.rm=TRUE)))
  
  # add extra spaces at the beginning for groups with same number of samples
  while(sum(duplicated(z)) > 0){
    z[duplicated(z)] = paste(" ", z[duplicated(z)], sep="")
  }
  
  factor(x, levels=y, labels=z)
  
}

make_n_factor_map = function(x){
  
  if(! is.factor(x))
    x = as.factor(x)
  
  y = levels=levels(x)
  z = sapply(y, function(u) sprintf("n=%d", sum(x == u, na.rm=TRUE)))
  
  data.frame(levels=y, count=z, levelsSuffixed=paste(y, " (", z, ")", sep=""))
}

lapply_c = function(x, ...){
  
  result = lapply(x, ...)
  names(result) = x
  
  result
  
}

