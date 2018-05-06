
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

getPredictionStats = function(predictions, truth, classes=NULL, decision_values=NULL){
  
  if(length(predictions) != length(truth))
    stop("Predictions and true label vectors are not of the same length")
  
  # accuracy
  accuracy = NULL
  # balanced accuracy = (sensitivity + specificity) / 2
  balanced_accuracy = NULL
  # sensitivity = accuracy of class 1 (positive class)
  sensitivity = NULL
  # specificity = accuracy of class 0
  specificity = NULL
  
  # If classes are not supplied, take the factor levels of Groups argument; 
  # assume the alphanumeric order is class 0 followed by class 1, so reverse the factor levels
  # since we expect the classes argument to be (class 1, class 0)
  
  if(is.null(classes))
    classes = rev(levels(factor(truth)))
  
  # calculate the stats
  accuracy = sum(predictions == truth)/length(predictions)
  sensitivity = sum(predictions[ truth == classes[1] ] == truth[ truth == classes[1] ])/sum(truth == classes[1])
  specificity = sum(predictions[ truth == classes[2] ] == truth[ truth == classes[2] ])/sum(truth == classes[2])
  balanced_accuracy = (sensitivity + specificity) / 2
  
  if(is.null(decision_values)){
    
    c(accuracy=accuracy, sensitivity=sensitivity, specificity=specificity, balanced_accuracy=balanced_accuracy)
    
  }else{
    
    auc = NULL
    
    tryCatch({
      auc = auc(roc(truth, decision_values, levels=classes, direction=">"))
    }, error = function(e){
      #print(e)
    })
    
    c(accuracy=accuracy, sensitivity=sensitivity, specificity=specificity, balanced_accuracy=balanced_accuracy, auc=auc)
  }
  
  
}
