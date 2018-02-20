

tryCatch({
  
}, error = function(e){ print(e) })

# make figures

library(knitr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(RColorBrewer)
library(ggrepel)

source("../src/util.R")
source("../src/plotutil.R")

source("util_1.R")
source("vars.R")

# =========================================================
#  Figure 2A
# =========================================================

tryCatch({
  
  loaded = load("obj/1_1.rda")
  df = df1_1[df1_1$Tissue %in% c("Colon"), ] 
  ggplot(df, aes(x=GroupsN, y=N, fill=Groups))+
    geom_boxplot()+ggtitle("A")+
    ylab("size of divergent set")+xlab("")+
    gtheme+THEME_X_HORIZ+
    scale_fill_manual(values=NORMAL_TUMOR_COLS)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center",
          legend.title=element_blank())+
    facet_wrap(~TissueSuffixed, nrow=1, scales="free_x")
  rm(df)
  rm(list=loaded)
  
  loaded = load("obj/")
  
  pdf("figure1.pdf", width=5, height=10)
  try(
  )
  
  dev.off()
  
  
  
  
}, error = function(e){ print(e) })

# =========================================================
#  Supplementary Figure 1A
# =========================================================

tryCatch({
  
  loaded = load("obj/1_1.rda")
  
  df = df1_1[! df1_1$Tissue %in% c("Colon", "Head and Neck region"), ] #, "Head and Neck region"
  
  pdf("figureS1.pdf", width=14, height=10)
  try(
    ggplot(df, aes(x=GroupsN, y=N, fill=Groups))+
    geom_boxplot()+ggtitle("A")+
    ylab("size of divergent set")+xlab("")+
    gtheme+THEME_X_HORIZ+
    scale_fill_manual(values=NORMAL_TUMOR_COLS)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center",
          legend.title=element_blank())+
    facet_wrap(~TissueSuffixed, nrow=1, scales="free_x")
  )
  dev.off()
  
  rm(df)
  rm(list=loaded)
  
}, error = function(e){ print(e) })


# =========================================================
# 
# =========================================================

tryCatch({
  
}, error = function(e){ print(e) })

# =========================================================
# 
# =========================================================

tryCatch({
  
}, error = function(e){ print(e) })

# =========================================================
# 
# =========================================================

tryCatch({
  
}, error = function(e){ print(e) })


