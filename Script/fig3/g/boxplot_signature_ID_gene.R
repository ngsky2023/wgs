rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggsci)
library("ggpubr")
library(ggsignif)

signatures<-read.table("ID1-2_ID83B_weight_vs_NRAS_BRAF.txt",header = T,quote = "",sep = "\t")
signatures$ID1_2<-rowSums(signatures[,c("ID1","ID2")])

genes=c("NRAS","BRAF")

signatures$NRAS_BRAF <-ifelse(signatures$Gene %in% genes,"BRAF/NRAS-driven tumor","Non-BRAF/NRAS-driven tumor") 

signatures$NRAS_BRAF <-factor(signatures$NRAS_BRAF,levels=c("Non-BRAF/NRAS-driven tumor","BRAF/NRAS-driven tumor"))

table(signatures$NRAS_BRAF)#Type	hypertension	anamnesis	MutNum
variables <- c("NRAS_BRAF")
burden <- c('ID1_2','ID83B')

for(i in burden){
  p_list <- lapply(1:length(variables),function(x){
    clin <- variables[x]
    ymax <- max(signatures[,i])
    data <- signatures[!is.na(signatures[,clin]),]
    
    pv <- wilcox.test(data[,i]~data[,clin])$p.value
    
    p <- ggplot()+
      geom_boxplot(data=data,aes(x=factor(data[,clin]),y=data[,i]),color="black",fill="white",outlier.shape=NA)+
      geom_jitter(data=data,aes(x=as.character(data[,clin]),y=data[,i],color=data[,clin]),width=0.15,size=0.6)+
      geom_segment(aes(x=1,y=ymax,xend=1,yend=ymax+ymax/20),color='black')+
      geom_segment(aes(x=2,y=ymax,xend=2,yend=ymax+ymax/20),color='black')+
      geom_segment(aes(x=1,y=ymax+ymax/20,xend=2,yend=ymax+ymax/20),color='black')+
      annotate("text",x=1.5,y=ymax+ymax/10,label=paste0("P = ",signif(pv,3)),size=2.7)+
      theme_bw()+
      theme(panel.border=element_blank(),
            panel.grid=element_blank(),
            axis.line=element_line(),
            legend.position='none',
            legend.key.size = unit(0.2, "lines"),
            legend.margin=margin(0,0,0,0),
            legend.background=element_rect(fill="transparent",color="transparent"),
            legend.text=element_text(size=8,color='black'),
            text=element_text(size=8,color='black'),
            axis.text.x=element_text(size=8,color='black',angle=-10,vjust = 0.5,hjust = 0.5),
            axis.text.y=element_text(size=8,color='black'),
            axis.title=element_text(size=8,color='black'),
            plot.margin=margin(5.5,0,5.5,0,unit="pt"))
    if ( i == "ID1_2"){
    p <- p+labs(x="",y=paste0("Proportion of ID1/2"),color="")+
    scale_color_manual(values=c("#8DAFDC","#506D9C"))
    }else if ( i == "ID83B"){
    p <- p+labs(x="",y=paste0("Proportion of ID83B"),color="")+
    scale_color_manual(values=c("#8DAFDC","#506D9C"))
    }
    
  })
  
   p_leg=as_ggplot(get_legend(p_list[[1]]))
   p_list=lapply(p_list,function(x) x=x+theme(legend.position="none") )
  
   plot_grid(plotlist=p_list,ncol=length(variables),rel_widths=c(1.1,rep(1,length(variables)-1)))
   ggsave(paste0("boxplot_",i,"_NRAS_BRAF.pdf"), width = length(variables)*2.8, height = 3)# dev.off()
}
