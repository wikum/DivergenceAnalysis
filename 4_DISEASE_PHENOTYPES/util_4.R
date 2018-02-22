
source("../src/util.R")

compute_p_mat = function(df, key, value, levels=NULL){
  
  key = df[, key]
  value = df[, value]
  
  if(! is.factor(key))
    key = factor(key)
  
  if(is.null(levels))
    levels = levels(key)
  
  m = length(levels)
  
  R = matrix(NA, nrow=m-1, ncol=m-1)
  rownames(R) = levels[1:(m-1)]
  colnames(R) = levels[2:m]
  
  for(si in 1:nrow(R)){
    for(ti in (si):ncol(R)){
      s = rownames(R)[si]
      t = colnames(R)[ti]
      R[s, t] = wilcox.test(value[which(key == s)], value[which(key == t)])$p.value
    }
  }
  
  R
}

get_mean_df = function(df, u="N", v="Groups", w="GroupsN", alpha=0.05){
  
  GroupsN = sapply(levels(df[, v]), function(x) as.character(df[which(df[, v] == x), w])[1])
  
  df = data.frame(
    mean=sapply(levels(df[, v]), function(x) mean(df[which(df[, v] == x), u])),
    Groups=factor(levels(df[, v]), ordered=TRUE),
    GroupsN=factor(GroupsN,
                   levels=GroupsN,
                   ordered=TRUE),
    sd=sapply(levels(df[, v]), function(x) sd(df[which(df[, v] == x), u])),
    n=sapply(levels(df[, v]), function(x) length(which(df[, v] == x)))
  )
  
  df$sem = df$sd/sqrt(df$n)
  df$me = qt(1-alpha/2, df$n)*df$sem
  
  df
}
