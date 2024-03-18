rm(list=ls())
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cowplot)
library(gridExtra)
require(grid)
for (l in 5){
  data <- read.delim2(file=paste0("SignatureEstimation_denovo_",l,"_signature.xls"),sep="\t",stringsAsFactors=F,header = T)
  # colnames(data) <- unlist(lapply(colnames(data),function(x) paste0(unlist(strsplit(x,'[.]'))[1],"[",unlist(strsplit(x,'[.]'))[2],"]",unlist(strsplit(x,'[.]'))[3],"-",unlist(strsplit(x,'[.]'))[4])))
  # data <- as.data.frame(t(data))
  colnames(data) <- paste0("RS",c(1:ncol(data)))
  # data <- data[,c("S1","S2","S5","S4","S3")]
  num2 <- ncol(data)
  num1 <- num2-1
  print(class(num1))
  print(class(num2))
  
  ###plot the bar of six type mutational directions
  a<-1:32
  b<-rep(1.5,32)
  bar_color=c("#880000","#98F898","#4682B4","#696969","#880000","#98F898","#4682B4","#696969","white")
  c<-c(rep("a",5),rep("b",5),rep("c",5),rep("d",1),rep("e",5),rep("f",5),rep("g",5),rep("h",1))
  bar<-as.data.frame(cbind(a,b,c),stringsAsFactors =F)
  bar$a<-as.numeric(bar$a)
  bar$b<-as.numeric(bar$b)
  bar$c<-factor(bar$c,levels = unique(bar$c))
  bar_label<-c("\nDel","\nDup","\nInv","\nT","\nDel","\nDup","\nInv","\nT")
  bar_x<-c(3,8,13,15.7,19,24,28.5,32)
  bar_plot<-ggplot(bar,aes(x=a,y=b,fill=c))+geom_bar(stat="identity",width = 1)+theme_bw()+theme(panel.grid=element_blank(),                                                                                                                               panel.border=element_blank(),plot.margin=unit(c(5,2.5,0,10), unit="mm"),axis.title = element_blank(),axis.text = element_blank(),                                                                                                                             axis.line = element_blank(),axis.ticks = element_blank(),legend.position = "none")+
    scale_fill_manual(values=bar_color)
  for(i in 1:8){
    bar_plot<-bar_plot+annotate("text",x=bar_x[i],y=1.7,label=bar_label[i],color="black",size=3)
  }
  
  data <- as.data.frame(data)
  data$type1<-as.character(rownames(data))
  data<-separate(data = data, col = type1, into = c("A5", "type1"), sep = "\\[")
  data<-separate(data = data, col = type1, into = c("type1","A3"), sep = "\\]")
  #data$type2<-paste(data$A5,data$A3,sep='.')
  data$type3 <- data$A3
  data<-subset(data,select = -c(A5,A3))
  data$type2 <- c(rep("clustered",16),rep("non-clustered",16))
  data1<-data[,c("type1","type2")]
  data2<-data[,1:(ncol(data)-1)]
  data<-as.data.frame(cbind(data1,data2))
  #prepare the x-axis text
  sig <- data
  #sig$type3 <- sig$type2
  sig$order <- 1:nrow(sig)
  sig$type1 <- factor(sig$type1,levels=unique(sig$type1))
  xlable<-paste(sig$type2,sig$type3)
  
  ###translate the formate of data to plot
  sig.t <- melt(sig,id.vars = c("type1","type2","type3","order"))
  sig.t$value<-as.numeric(sig.t$value)
  #fmt_dcimals <- function(decimals=0){
  # return a function responpsible for formatting the
  # axis labels with a given number of decimals
  #function(x) as.character(round(x,decimals))
  #}
  
  ylims <- list(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
  ylable <- c(1,1,1,1,1,1,1)
  bar<-as.matrix(c(1,1,1,1,1,1))
  
  ###plot signature bar for n(signature)-1, which needs to adjust aitificialy for each Sig 
  p1 <- lapply(1:num1,function(x){
    data <- sig.t[which(sig.t$variable==unique(sig.t$variable)[x]),]
    p<-ggplot(data,aes(x=order,y=value,fill=type1))+
      geom_bar(width=0.6,stat="identity")+
      scale_x_continuous(expand = c(0,0))+
      labs(y="",x="")+
      lims(y=ylims[[x]])+
      scale_fill_manual(values=c("#880000","#98F898","#4682B4","#696969","#880000","#98F898","#4682B4","#696969"))+
      annotate("text",x=2,y=ylable[x]-0.05,label=unique(sig.t$variable)[x])+
      theme_bw()+
      theme(plot.margin=unit(c(0,9,0,2), unit="mm"),
            legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            # axis.ticks.x=element_blank(),
            axis.line.y=element_line(),
            text=element_text(size=12),
            axis.text.y=element_text(size=12,color="black"),
            plot.title=element_text(hjust=0.5),
            strip.text = element_text(colour = 'white',face='bold',size=12))
    # if(x==1 |x==3){
    #   p<-p+scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8),labels=c("0.00","0.20","0.40","0.60","0.80"))+
    #     annotate("text",x=2,y=0.72,label=unique(sig.t$variable)[x])
    # }else{
    #   p <- p+annotate("text",x=2,y=ylable[x]-0.01,label=unique(sig.t$variable)[x])
    # }
    g <- ggplot_gtable(ggplot_build(p))
    strip_both <- which(grepl('strip-', g$layout$name))
    fills <- c("#1EBFF0","#05070B","#E62725","#CBCACB","#A1CF64","#EDC8C5")
    k <- 1
    for (i in strip_both) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    return(g)
  })
  
  ###plot the last signature
  p1_2 <- lapply(num2,function(x){
    data <- sig.t[which(sig.t$variable==unique(sig.t$variable)[x]),]
    print(nrow(data))
    print(length(xlable))
    p<-ggplot(data,aes(x=order,y=value,fill=type1))+
      geom_bar(width=0.6,stat="identity")+
      #facet_grid(~type1)+
      scale_x_continuous(breaks=1:nrow(data),labels=xlable,expand = c(0,0))+
      labs(y="",x="")+
      lims(y=ylims[[x]])+
      scale_fill_manual(values=c("#880000","#98F898","#4682B4","#696969","#880000","#98F898","#4682B4","#696969"))+
      
      annotate("text",x=2,y=ylable[x]-0.02,label=unique(sig.t$variable)[x])+
      theme_bw()+
      theme(plot.margin=unit(c(0,9,0,2), unit="mm"),legend.position="none",
            panel.grid=element_blank(),
            panel.border=element_blank(),
            # axis.text.x=element_text(angle=90,size=8,vjust=0.5),
            axis.text.x=element_text(angle=90,size=5,vjust=0.5,hjust = 1,color="black"),
            axis.line=element_line(),
            text=element_text(size=12),
            axis.text.y=element_text(size=12,color="black"),
            plot.title=element_text(hjust=0.5),
            strip.text = element_text(colour = 'white',face='bold',size=12),
            strip.background =element_rect(fill="white",color="white"))
    #p<-p+scale_y_continuous(breaks=c(0,0.02,0.04,0.06,0.08,0.1),labels=c(0,0.02,0.04,0.06,0.08,0.1))
  })
  
  ###plot the final figure
  pdf(file = paste0("denovoSignatures_sv_",l,"_plot.pdf"),width=6,height=num2+1.2)
  p <- plot_grid(plotlist=c(p1,p1_2),ncol=1,rel_heights=c(rep(0.85,(num2-1)),(0.75+1)))
  pp<-grid.arrange(bar_plot,p,ncol=1,heights=c((num2+3)*0.06,(num2+1)*0.925))
  grid.draw(pp)
  dev.off()
}
