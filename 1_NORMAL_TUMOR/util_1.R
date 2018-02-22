
source("../src/util.R")

make_facet_df = function(dfList, suffix="", prefix=""){
  
  temp = lapply(1:length(dfList), function(i){ 
    
    df = data.frame(N=dfList[[i]]$N, 
                    Groups=dfList[[i]]$Groups,
                    GroupsN=make_n_factor(dfList[[i]]$Groups),
                    Tissue=factor(rep(names(dfList)[i], length(dfList[[i]]$Groups))))

    df$TissueSuffixed =  factor(paste(df$Tissue, suffix, sep="\n"))
    df$TissuePrefixed =  factor(paste(prefix, df$Tissue, sep="\n"))
    
    df
  })
  
  names(temp) = names(dfList)
  
  Reduce(rbind, temp)
  
}
