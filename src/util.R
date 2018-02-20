

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


