rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggsci)
library("ggpubr")
library(ggsignif)

tmb<-read.table("PUCH_AMGP.xls",header = T,quote = "",sep = "\t")
colnames(tmb)=c("sample","TMB","sv_number","cnv_length","Cohort")
tmb$sv_number=log(tmb$sv_number)

tmb$Cohort <- factor(tmb$Cohort,levels = c("PUCH","AMGP"))

table(tmb$Cohort)#Type	hypertension	anamnesis	MutNum
variables <- c("Cohort")
burden <- c('TMB','cnv_length','sv_number')

for(i in burden){
  p_list <- lapply(1:length(variables),function(x){
    clin <- variables[x]
    ymax <- max(tmb[,i])
    data <- tmb[!is.na(tmb[,clin]),]
    
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
            axis.text=element_text(size=8,color='black'),
            axis.title=element_text(size=8,color='black'),
            plot.margin=margin(5.5,0,5.5,0,unit="pt"))
  if((x %% 2) != 1){
    if ( i == "TMB"){
    p <- p+labs(x="",y=paste0("Mutation per megabase"),color="")+
    scale_color_manual(values=c("#8DB8FF","#5371A3"))
    }else if ( i == "sv_number"){
    p <- p+labs(x="",y=paste0("Number of SVs\n(log scale)"),color="")+
    scale_color_manual(values=c("#8DB8FF","#5371A3"))
    }else if ( i == "cnv_length"){
    p <- p+labs(x="",y=paste0("CNV(except CN2 CN3-5) %genome"),color="")+
    scale_color_manual(values=c("#8DB8FF","#5371A3"))
    }
  }else{
    if ( i == "TMB"){
    p <- p+labs(x="",y=paste0("Mutation per megabase"),color="")+
    #scale_color_manual(values=c("#8A98B8","#F09386"))
    scale_color_manual(values=c("#003C9D","#FF3333"))
    }else if ( i == "sv_number"){
    p <- p+labs(x="",y=paste0("Number of SVs\n(log scale)"),color="")+
    #scale_color_manual(values=c("#8A98B8","#F09386"))
    scale_color_manual(values=c("#003C9D","#FF3333"))
    }else if ( i == "cnv_length"){
    p <- p+labs(x="",y=paste0("CNV(except CN2 CN3-5) %genome"),color="")+
    #scale_color_manual(values=c("#8A98B8","#F09386"))
    scale_color_manual(values=c("#003C9D","#FF3333"))
    }
  }
    
  })
  
   p_leg=as_ggplot(get_legend(p_list[[1]]))
   p_list=lapply(p_list,function(x) x=x+theme(legend.position="none") )
  
   plot_grid(plotlist=p_list,ncol=length(variables),rel_widths=c(1.1,rep(1,length(variables)-1)))
   ggsave(paste0("boxplot_",i,"_cohort.pdf"), width = length(variables)*2, height = 3)# dev.off()
}
