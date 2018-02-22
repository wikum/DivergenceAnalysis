

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

source("src/util.R")
source("src/plotutil.R")

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
    themeGENERIC()+
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
    themeGENERIC()+ggtitle("B")+
    xlab("")+ylab("size of divergent set")+
    scale_fill_manual(name="", values=cols)+
    facet_grid(~set, scales="free_x", space="free_x")+
    guides(fill=FALSE)+
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1, size=15),
          plot.title=element_text(hjust=0, size=30))
  rm(df)
  rm(list=loaded)

  pdf("figures/figure2.pdf", width=16, height=6)
  try(
    grid.arrange(plot2A, plot2B, layout_matrix=matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2), nrow=1))
  )
  dev.off()
  
}, error = function(e){ print(e) })

# =========================================================
#  Figure 3
# =========================================================

tryCatch({
  
  loaded = load("5_FEATURE_SELECTION/obj/5_1.rda")
  plot5A = ggplot(df, aes(x=PA, y=PB, shape=sig, alpha=sig, col=sig))+
    geom_point(size=2)+
    xlab("Luminal A")+ylab("Luminal B")+
    ggtitle("A")+
    geom_text_repel(data=df[df$ssig==TRUE,], aes(x=PA, y=PB, label=gene), 
                    col="maroon", size=6)+
    scale_shape_manual(values=16:17)+
    scale_alpha_manual(values=c(.2, .9))+
    scale_color_manual(#name=expression(chi^2 ~ test ~ P <= 0.05 ~ (Bonferroni)), 
      values=c("lightblue", "steelblue"))+
    geom_abline(slope=1, size=.3)+
    themeGENERIC()+guides(alpha=FALSE, shape=FALSE, col=FALSE)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center")
  rm(list=loaded)
  
  loaded = load("5_FEATURE_SELECTION/obj/5_2.rda")
  plot5B = ggplot(df, aes(x=PA, y=PB, shape=sig, alpha=sig, col=sig))+
    geom_point(size=2)+
    ylab("Prostate")+xlab("Other")+
    ggtitle("B")+
    #xlim(0, .2)+
    geom_text_repel(data=df[df$ssig==TRUE,], aes(x=PA, y=PB, label=gene),
                    col="maroon", size=6)+
    scale_shape_manual(#name=expression(chi^2 ~ test ~ P <= 0.05 ~ (Bonferroni)), 
      values=16:17)+
    scale_alpha_manual(values=c(.2, .9))+
    scale_color_manual(values=c("lightblue", "steelblue"))+
    geom_abline(slope=1, size=.3)+
    themeGENERIC()+guides(color=FALSE, alpha=FALSE, shape=FALSE)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center")
  rm(list=loaded)
  
  pdf("figures/figure3.pdf", width=16, height=8)
  try(grid.arrange(plot5A, plot5B, ncol=2))
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
    geom_boxplot(outlier.size=1, lwd=1)+themeGENERIC()+
    theme(axis.text.x=element_text(angle=90, hjust=1))+
    scale_fill_manual(values=cols)+
    ggtitle("A")+xlab("tissue")+ylab("size of divergent set")+guides(fill=FALSE)+
    THEME_TITLE_LEFT
  
  df = dfList$`Skin-Sun Exposed`
  
  plot6B = ggplot(df, aes(x=reorder(Groups, N, median), y=N, fill=trainGroup))+
    geom_boxplot(outlier.size=1, lwd=1)+themeGENERIC()+
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
#  Supplementary Figure 1
# =========================================================

tryCatch({
  
  loaded = load("1_NORMAL_TUMOR/obj/1_1.rda")
  df = df1_1[! df1_1$Tissue %in% c("Colon"), ] #, "Head and Neck region"
  plotS1A = ggplot(df, aes(x=GroupsN, y=N, fill=Groups))+
    geom_boxplot()+ggtitle("A")+
    ylab("size of divergent set")+xlab("")+
    themeGENERIC()+
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
    themeGENERIC()+
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
    themeGENERIC()+
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
#  Supplementary Figure 2
# =========================================================

tryCatch({

  loaded = load("5_FEATURE_SELECTION/obj/5_3.rda")
  plotS2A = ggplot(df, aes(x=PA, y=PB, shape=sig, alpha=sig, col=sig))+
    geom_point(size=2)+
    xlab("Lung")+ylab("Other")+
    ggtitle("A")+
    geom_text_repel(data=df[df$ssig==TRUE,], aes(x=PA, y=PB, label=gene), 
                    col="maroon", size=6)+
    scale_shape_manual(values=16:17)+
    scale_alpha_manual(values=c(.2, .9))+
    scale_color_manual(#name=expression(chi^2 ~ test ~ P <= 0.05 ~ (Bonferroni)), 
      values=c("lightblue", "steelblue"))+
    geom_abline(slope=1, size=.3)+
    themeGENERIC()+guides(alpha=FALSE, shape=FALSE, col=FALSE)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center")  
  rm(list=loaded)
  
  loaded = load("5_FEATURE_SELECTION/obj/5_4.rda")
  plotS2B = ggplot(df, aes(x=PA, y=PB, shape=sig, alpha=sig, col=sig))+
    geom_point(size=2)+
    xlab("Breast")+ylab("Other")+
    ggtitle("B")+
    geom_text_repel(data=df[df$ssig==TRUE,], aes(x=PA, y=PB, label=gene), 
                    col="maroon", size=6)+
    scale_shape_manual(values=16:17)+
    scale_alpha_manual(values=c(.2, .9))+
    scale_color_manual(#name=expression(chi^2 ~ test ~ P <= 0.05 ~ (Bonferroni)), 
      values=c("lightblue", "steelblue"))+
    geom_abline(slope=1, size=.3)+
    themeGENERIC()+guides(alpha=FALSE, shape=FALSE, col=FALSE)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center")  
  rm(list=loaded)
  
  pdf("figures/figureS2.pdf", width=16, height=8)
  try(grid.arrange(plotS2A, plotS2B, ncol=2, nrow=1))
  dev.off()
  
}, error = function(e){ print(e) })




# =========================================================
#  Supplementary Figure 3
# =========================================================

tryCatch({
  
  loaded = load("5_FEATURE_SELECTION/obj/5_5.rda")
  plotS3A = ggplot(df, aes(x=PA, y=PB, shape=sig, alpha=sig, col=sig))+
    geom_point(size=2)+
    xlab("Other")+ylab("Lung")+
    ggtitle("A")+
    geom_text_repel(data=df[df$ssig==TRUE,], aes(x=PA, y=PB, label=gene), 
                    col="maroon", size=6)+
    scale_shape_manual(values=16:17)+
    scale_alpha_manual(values=c(.2, .9))+
    scale_color_manual(#name=expression(chi^2 ~ test ~ P <= 0.05 ~ (Bonferroni)), 
      values=c("lightblue", "steelblue"))+
    geom_abline(slope=1, size=.3)+
    themeGENERIC()+guides(alpha=FALSE, shape=FALSE, col=FALSE)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center")
  rm(list=loaded)
  
  loaded = load("5_FEATURE_SELECTION/obj/5_6.rda")
  plotS3B = ggplot(df, aes(x=PA, y=PB, shape=sig, alpha=sig, col=sig))+
    geom_point(size=2)+
    xlab("Other")+ylab("Breast")+
    ggtitle("B")+
    geom_text_repel(data=df[df$ssig==TRUE,], aes(x=PA, y=PB, label=gene), 
                    col="maroon", size=6)+
    scale_shape_manual(values=16:17)+
    scale_alpha_manual(values=c(.2, .9))+
    scale_color_manual(#name=expression(chi^2 ~ test ~ P <= 0.05 ~ (Bonferroni)), 
      values=c("lightblue", "steelblue"))+
    geom_abline(slope=1, size=.3)+
    themeGENERIC()+guides(alpha=FALSE, shape=FALSE, col=FALSE)+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="horizontal",
          legend.justification="center")
  rm(list=loaded)
  
  pdf("figures/figureS3.pdf", width=16, height=8)
  try(grid.arrange(plotS3A, plotS3B, ncol=2, nrow=1))
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
    themeGENERIC()+
    theme(
      axis.text.x=element_text(angle=90, hjust=1),
      plot.title=element_text(hjust=0, vjust=1, size=30), 
      legend.position="bottom", 
      legend.direction="horizontal",
      legend.justification="center",
      legend.title=element_blank())
  
  rm(list=loaded)
  
  pdf("figures/figureS4.pdf", width=18, height=12)
  try(
    grid.arrange(plotS4)
  )
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
    themeGENERIC()+ggtitle("")+
    theme(plot.title = element_text(colour="grey20",size=30,angle=0,face="plain", hjust=-.35),
          axis.title.x=element_text(colour="grey20",size=25,angle=0,face="bold"),
          axis.title.y=element_text(colour="grey20",size=25,angle=90,face="bold"),
          axis.text.x=element_text(colour="grey20",size=20,angle=90,face="plain", hjust=1),
          legend.position="bottom", legend.justification="center", legend.key.width=unit(45, "pt"))
  rm(list=loaded)
  
  pdf("figures/figureS5.pdf", width=24, height=24)
  try(
    grid.arrange(plotS5)
  )
  dev.off()
  
}, error = function(e){ print(e) })

# =========================================================
#  Supplementary figure 6
# =========================================================

tryCatch({
  
  loaded = load("4_DISEASE_PHENOTYPES/obj/4_2.rda")
  plotS6A = ggplot(dfMeanList$GLEASON, 
                aes(x=GroupsN, y=mean, fill=Groups))+
    geom_col()+ggtitle("A")+
    geom_errorbar(aes(ymin=mean-me, ymax=mean+me), width=.2)+
    ylab("average size of divergent set")+xlab("")+
    themeGENERIC()+
    scale_fill_manual(name="Gleason (Primary)", values=brewer.pal(3, "Set2"))+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="vertical",
          legend.justification="center",
          legend.key.height=unit(2,"line"))
  plotS6B = ggplot(dfMeanList$SMOKING, aes(x=GroupsN, y=mean, fill=Groups))+
    geom_col()+ggtitle("B")+
    geom_errorbar(aes(ymin=mean-me, ymax=mean+me), width=.2)+
    ylab("average size of divergent set")+xlab("")+
    themeGENERIC()+
    scale_fill_manual(name="Smoking", values=brewer.pal(4, "Set3"))+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="vertical",
          legend.justification="center",
          legend.key.height=unit(2,"line"))+
    scale_y_continuous(labels=scales::comma)  
  plotS6C = ggplot(dfMeanList$HIST, aes(x=GroupsN, y=mean, fill=Groups))+
    geom_col()+ggtitle("C")+
    geom_errorbar(aes(ymin=mean-me, ymax=mean+me), width=.2)+
    ylab("average size of divergent set")+xlab("")+
    themeGENERIC()+
    scale_fill_manual(name="Histological grade", values=c("yellow2", "turquoise4", "steelblue"))+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="vertical",
          legend.justification="center",
          legend.key.height=unit(2, "line"))  
  rm(list=loaded)
  loaded = load("4_DISEASE_PHENOTYPES/obj/4_3.rda")
  plotS6D = ggplot(dfMean, aes(x=GroupsN, y=mean, fill=Groups))+
    geom_col()+ggtitle("D")+
    geom_errorbar(aes(ymin=mean-me, ymax=mean+me), width=.2)+
    ylab("average size of divergent set")+xlab("")+
    themeGENERIC()+
    scale_fill_manual(name="Type", values=c("lightpink1", "maroon"))+
    theme(plot.title=element_text(hjust=0, size=30), 
          legend.position="bottom", 
          legend.direction="vertical",
          legend.justification="center",
          legend.key.height=unit(2, "line"))
  rm(list=loaded)
  
  pdf("figures/figureS6.pdf", width=24, height=10)
  try(
    print(plot_grid(plotS6A, plotS6B, plotS6C, plotS6D, align="h", nrow=1, ncol=4))
  )
  dev.off()

}, error = function(e){ print(e) })


# =========================================================
#  Supplementary figure 7
# =========================================================

tryCatch({

  loaded = load("4_DISEASE_PHENOTYPES/obj/4_1.rda")
  plotS7A = ggplot(df,
                aes(x=N_RNASEQ, y=N_METHYL))+
    geom_point(col="violetred2")+
    ggtitle("A")+
    xlab("size of divergent set\n(expression profile)")+
    ylab("size of divergent set\n(methylation profile)")+
    themeGENERIC()+
    theme(legend.position="bottom", legend.direction="vertical",
          plot.title=element_text(hjust=0, size=30),
          legend.title=element_blank(), 
          legend.box.spacing=unit(10,"line"),
          legend.justification="center")+
    scale_y_continuous(labels=scales::comma)  
  rm(list=loaded)
  
  loaded = load("4_DISEASE_PHENOTYPES/obj/4_4.rda")
  plotS7B = ggplot(df,
                aes(x=N_RNASEQ, y=N_METHYL))+
    geom_point(col="olivedrab3")+
    ggtitle("B")+
    xlab("size of divergent set\n(expression profile)")+
    ylab("size of divergent set\n(methylation profile)")+
    themeGENERIC()+
    theme(legend.position="bottom", legend.direction="vertical",
          plot.title=element_text(hjust=0, size=30),
          legend.title=element_blank(), 
          legend.box.spacing=unit(10,"line"),
          legend.justification="center")+
    scale_y_continuous(labels=scales::comma)  
  rm(list=loaded)
  
  pdf("figures/figureS7.pdf", width=12, height=6)
  try(
    print(plot_grid(plotS7A, plotS7B, align="h", nrow=1, ncol=2))
  )
  dev.off()
  
  
}, error = function(e){ print(e) })







