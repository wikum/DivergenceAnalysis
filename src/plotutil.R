
require(ggplot2)

themeGENERIC = function(){
  theme(axis.text.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
               axis.text.y = element_text(colour="grey20",size=20,angle=0,hjust=1,vjust=0,face="plain"),  
               axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
               axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
               plot.title = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
               legend.text = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
               legend.title = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
               legend.key.height=unit(1.5, "line"),
               panel.grid.major = element_line(colour="darkgray"))
}

themeNO_X_LABS = function(){
  theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
}

themeX_HORIZ = function(){
  theme(axis.text.x = element_text(angle=0))
}

themeX_VERT = function(){
  theme(axis.text.x = element_text(angle=90))
}

