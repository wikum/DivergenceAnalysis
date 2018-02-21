

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
library(rutils)

# =========================================================
#  Figure 2
# =========================================================

tryCatch({
  
  NORMAL_TUMOR_COLS = c("dodgerblue3", "salmon")
  
  
  loaded = load("1_NORMAL_TUMOR/obj/1_1.rda")
  df1 = df1_1[df1_1$Tissue %in% c("Colon"), ] 
  rm(list=loaded)
  loaded = load("1_NORMAL_TUMOR/obj/1_2.rda")
  df2 = df1_2[df1_2$Tissue %in% c("Head and Neck"), ] 
  rm(list=loaded)
  loaded = load("1_NORMAL_TUMOR/obj/1_3.rda")
  df3 = df1_3[df1_3$Tissue %in% c("GPL96"), ]
  rm(list=loaded)
  cols = c("N", "Groups", "GroupsN", "Tissue", "TissueSuffixed", "TissuePrefixed")
  df = rbind(df1[, cols], df2[, cols], df3[, cols])
  df$set = factor(c(as.character(df1$TissueSuffixed), 
                    as.character(df2$TissueSuffixed), 
                    as.character(df3$TissuePrefixed)),
                  levels=c("Colon\nRNASeq", "Head and Neck\nMethylation 450k", "Breast\nGPL96")
  )
  plot2A = ggplot(df, aes(x=GroupsN, y=N, fill=Groups))+
    geom_boxplot()+ggtitle("A")+
    ylab("size of divergent set")+xlab("")+
    gtheme.GENERIC()+
    scale_fill_manual(values=NORMAL_TUMOR_COLS)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center",
          legend.title=element_blank())+
    facet_wrap(~set, nrow=1, ncol=6, scales="free")  

  
  loaded = load("6_HALLMARKS/obj/6_1.rda")
  NORMAL_palette = c("dodgerblue3")
  PAM50_palette = c("maroon", "yellow4", "turquoise4", "cyan", "magenta")
  ER_palette = c("pink", "maroon")
  cols = c(NORMAL_palette, PAM50_palette, ER_palette)
  plot2B = ggplot(df, aes(x=group_n, y=div, fill=group))+
    geom_boxplot()+
    gtheme.GENERIC()+ggtitle("B")+
    xlab("")+ylab("size of divergent set")+
    scale_fill_manual(name="", values=cols)+
    facet_grid(~set, scales="free_x", space="free_x")+
    guides(fill=FALSE)+
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=15),
          plot.title=element_text(hjust=0, size=30))
  rm(df)
  rm(list=loaded)

  pdf("figures/figure1.pdf", width=16, height=6)
  try(
    grid.arrange(plot2A, plot2B, layout_matrix=matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2), nrow=1))
  )
  dev.off()
  
}, error = function(e){ print(e) })

# =========================================================
#  Supplementary Figure 1
# =========================================================

tryCatch({
  
  loaded = load("1_NORMAL_TUMOR/obj/1_1.rda")
  df = df1_1[! df1_1$Tissue %in% c("Colon"), ] #, "Head and Neck region"
  plotS1A = ggplot(df, aes(x=GroupsN, y=N, fill=Groups))+
    geom_boxplot()+ggtitle("A")+
    ylab("size of divergent set")+xlab("")+
    gtheme.GENERIC()+
    scale_fill_manual(values=NORMAL_TUMOR_COLS)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center",
          legend.title=element_blank())+
    facet_wrap(~TissueSuffixed, nrow=1, scales="free_x")
  rm(df)
  rm(list=loaded)
  
  loaded = load("1_NORMAL_TUMOR/obj/1_2.rda")
  df = df1_2[! df1_2$Tissue %in% c("Head and Neck"), ]
  plotS1B = ggplot(df, aes(x=GroupsN, y=N, fill=Groups))+
    geom_boxplot()+ggtitle("B")+
    ylab("size of divergent set")+xlab("")+
    gtheme.GENERIC()+
    scale_fill_manual(values=NORMAL_TUMOR_COLS)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center",
          legend.title=element_blank())+
    facet_wrap(~TissueSuffixed, nrow=1, ncol=6, scales="free_x")
  rm(df)
  rm(list=loaded)
  
  loaded = load("1_NORMAL_TUMOR/obj/1_3.rda")
  df = df1_3[! df1_3$Tissue %in% c("GPL96"), ]
  plotS1C = ggplot(df, aes(x=GroupsN, y=N, fill=Groups))+
    geom_boxplot()+ggtitle("C")+
    ylab("size of divergent set")+xlab("")+
    gtheme.GENERIC()+
    scale_fill_manual(values=NORMAL_TUMOR_COLS)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center",
          legend.title=element_blank())+
    facet_wrap(~TissuePrefixed, nrow=1, ncol=6, scales="free_x")
  rm(df)
  rm(list=loaded)

  pdf("figures/figureS1.pdf", width=18, height=18)
  try(
    grid.arrange(plotS1A, plotS1B, plotS1C,
               layout_matrix=rbind(c(1, 1, 1, 1, 1, 1), c(2, 2, 3, 3, NA, NA)))
  )
  dev.off()
  
}, error = function(e){ print(e) })


# =========================================================
# Figure 6
# =========================================================

tryCatch({
  
  THEME_TITLE_LEFT = theme(plot.title = element_text(colour="grey20",size=30,angle=0,face="plain", hjust=0))
  
  loaded = load("2_GTEX_TISSUE/obj/2_1.rda")
  cols=c("turquoise4", "yellow2")
  
  df = dfList$`Breast-Mammary Tissue`
  
  plot6A = ggplot(df, aes(x=reorder(Groups, N, median), y=N, fill=trainGroup))+
    geom_boxplot(outlier.size=1, lwd=1)+gtheme.GENERIC()+
    theme(axis.text.x=element_text(angle=90, hjust=1))+
    scale_fill_manual(values=cols)+
    ggtitle("A")+xlab("tissue")+ylab("size of divergent set")+guides(fill=FALSE)+
    THEME_TITLE_LEFT
  
  df = dfList$`Skin-Sun Exposed`
  
  plot6B = ggplot(df, aes(x=reorder(Groups, N, median), y=N, fill=trainGroup))+
    geom_boxplot(outlier.size=1, lwd=1)+gtheme.GENERIC()+
    theme(axis.text.x=element_text(angle=90, hjust=1))+
    scale_fill_manual(values=cols)+
    ggtitle("B")+xlab("tissue")+ylab("size of divergent set")+guides(fill=FALSE)+
    THEME_TITLE_LEFT

  rm(list=loaded)
  
  pdf("figures/figure6.pdf", width=10, height=8)
  try(grid.arrange(plot6A, plot6B, ncol=2))
  dev.off()
  
}, error = function(e){ print(e) })

# =========================================================
#  Supplementary figure 4
# =========================================================

tryCatch({
  
  loaded = load("2_GTEX_TISSUE/obj/2_2.rda")
  
  plotS4 = ggplot(df, aes(x=reorder(GroupsSuffixed, N, median), y=N, fill=trainGroup))+
    geom_boxplot(outlier.size=1, lwd=1)+
    scale_fill_manual(values=c("turquoise4", "yellow2"))+
    ggtitle("")+xlab("tissue subtype")+
    ylab("size of divergent set")+guides(fill=FALSE)+
    gtheme.GENERIC()+
    theme(
      axis.text.x=element_text(angle=90, hjust=1),
      plot.title=element_text(hjust=0, vjust=1, size=30), 
      legend.position="bottom", 
      legend.direction="horizontal",
      legend.justification="center",
      legend.title=element_blank())
  
  rm(list=loaded)
  
  pdf("figures/figureS4.pdf", width=18, height=12)
  try(plotS4)
  dev.off()
  
}, error = function(e){ print(e) })

# =========================================================
#  Supplementary figure 5
# =========================================================

tryCatch({
  
  loaded = load("2_GTEX_TISSUE/obj/2_3.rda")
  remove_leading_zero = function(x) sub('^(-)?0[.]', '\\1.', x)
  plotS5 = ggplot(meltedAccMat, aes(x=colKey, y=rowKey, fill=value))+geom_tile(col="gray50")+
    scale_fill_gradient(low="white", high="steelblue", name="proportion of samples\npredicted")+
    geom_text(data=subset(meltedAccMat, value > 0.01), 
              aes(x=colKey, y=rowKey, label=remove_leading_zero(round(value, 2))), size=5)+
    xlab("predicted subtype")+ylab("true subtype")+
    gtheme.GENERIC()+ggtitle("")+
    theme(plot.title = element_text(colour="grey20",size=30,angle=0,face="plain", hjust=-.35),
          axis.title.x=element_text(colour="grey20",size=25,angle=0,face="bold"),
          axis.title.y=element_text(colour="grey20",size=25,angle=90,face="bold"),
          axis.text.x=element_text(colour="grey20",size=20,angle=90,face="plain", hjust=1),
          legend.position="bottom", legend.justification="center", legend.key.width=unit(45, "pt"))
  rm(list=loaded)
  
  pdf("figures/figureS5.pdf", width=24, height=24)
  try(plotS5)
  dev.off()
  
}, error = function(e){ print(e) })


