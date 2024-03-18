rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggsci)
library("ggpubr")
library(ggsignif)

tmb<-read.table("mutation_burden.xls",header = T,quote = "",sep = "\t")
tmb <- tmb[,c(1,5)]
names(tmb)[1] <- "sample"
wgd <- read.table("ALL.purity.txt",header = T,quote = "",sep = "\t")
wgd <- wgd[,c(1,18)]
names(wgd) <- c("sample","WGD")
wgd$WGD[wgd$WGD=="TRUE"] <- "WGD"
wgd$WGD[wgd$WGD=="FALSE"] <- "non-WGD"

df <- read.table("clinical_55.txt",header = T,quote = "",sep = "\t")

names(df)[1] <- "sample"
df$stage[df$Stage %in% c("I","II")]<-"I+II"
df$stage[df$Stage %in% c("III","IV")]<-"III+IV"
df$age[df$Age >=58 ]<-">=58"
df$age[df$Age <58]<-"<58"

#################
sv <- read.table("sv_filter_BND_number.xls",sep="\t", quote="", header=T,stringsAsFactors = F)
names(sv) <- c("sample","TRA","DEL","DUP","INV")
sv$sv_number <- rowSums(sv[,c(2,3,4,5)])
# sv <- sv[,c(1,6)]
tmb <- merge(tmb,sv,by="sample")
tmb <- merge(df,tmb,by="sample")
cnv <- read.table("cnv_type_length.xls",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(cnv) <- c("sample","Deletion","Loss_CN1","Copy_Neutral_LOH","CN2","CN3-5","CN6")
cnv$cnv_length <- rowSums(cnv[,c(2:7)])
cnv[,c(2:7)] <- cnv[,c(2:7)]/rowSums(cnv[,c(2:7)])
cnv[,"CNV(except CN2)"] <- cnv$`CN3-5`+cnv$CN6 + cnv$Deletion + cnv$Loss_CN1 +cnv$Copy_Neutral_LOH
cnv[,"CNV(except CN2 CN3-5)"] <-cnv$CN6 + cnv$Deletion + cnv$Loss_CN1+cnv$Copy_Neutral_LOH

tmb <- merge(tmb,cnv,by = "sample")
tmb <- merge(tmb,wgd,by="sample")
##规范数据


tmb$Specimen.Type <- factor(tmb$Specimen.Type,levels = c("primary","metastasis"))
tmb$Gender <- factor(tmb$Gender,levels = c("F","M"))
tmb$Primary.Site[which(tmb$Primary.Site=="other")]<- "other(foot)"
tmb$Primary.Site <- factor(tmb$Primary.Site,levels = c("hand","sole","heel","other(foot)"))
tmb$WGD <- factor(tmb$WGD,levels = c("non-WGD","WGD"))
tmb$sv_number <- log(tmb$sv_number)
tmb$INV <- log(tmb$INV)
tmb$DEL <- log(tmb$DEL)
tmb$DUP <- log(tmb$DUP)
tmb$TRA <- log(tmb$TRA)

table(tmb$Specimen.Type)#Type	hypertension	anamnesis	MutNum
variables <- c("age","Gender","Specimen.Type","stage","WGD")
# burden <- c("perMb","sv_number","cnv_length",colnames(cnv[2:7]))
burden <- c('CNV(except CN2)','CNV(except CN2 CN3-5)')
#
# p_list=lapply(seq(1,length(loc),2),function(x){
for(n in 1:length(burden)){
  i <- burden[n]
  p_list <- lapply(1:length(variables),function(x){
    clin <- variables[x]
    ymax <- max(tmb[,i])
    data <- tmb[!is.na(tmb[,clin]),]
    
    
    pv <- wilcox.test(data[,i]~data[,clin])$p.value
    # pv <- t.test(data$perMb~data[,clin])$p.value
    
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
  if((x %% 2)== 1){
    p <- p+labs(x="",y=paste0("Percentage of ",i),color="")+
    scale_color_manual(values=c("#8DB8FF","#5371A3"))
  }else{
    p <- p+labs(x="",y=paste0("Percentage of ",i),color="")+
      scale_color_manual(values=c("#8A98B8","#F09386"))
  }
  
    
  })
  
  # pdf("AM_boxplot_TMB_clinical.pdf",width=length(variables)*2-1,height=3)
  # p_list[[6]] <- p2
  p_leg=as_ggplot(get_legend(p_list[[1]]))
  p_list=lapply(p_list,function(x) x=x+theme(legend.position="none") )
  # pdf(paste0("AM_boxplot_",i,"_clinical.pdf"),width=8,height =8 )
  plot_grid(plotlist=p_list,ncol=3)
  ggsave(paste0("AM_boxplot_",i,"_clinical.pdf"), width = floor(length(variables)/3)*5, height = ceiling(length(variables)/9)*5)# dev.off()
}
#########################四组比较
library(reshape2)
data2 <- tmb[,c('sample','Primary.Site',"Deletion","Loss_CN1","Copy_Neutral_LOH","CN2","CN3-5","CN6",'CNV(except CN2)','CNV(except CN2 CN3-5)')]
burden <- c("Deletion","Loss_CN1","Copy_Neutral_LOH","CN2","CN3-5","CN6",'CNV(except CN2)','CNV(except CN2 CN3-5)')
for(n in 1:length(burden)){
  i <- burden[n]
  clin <- "Primary.Site"
  ymax <- max(data2[,i])
  data <- data2[!is.na(data2[,clin]),]
  # data$cnv_count <- log(data$cnv_count)
  my_comparisons=list(c("hand","sole"),c("hand","heel"),c("hand","other(foot)"),c("sole","heel"),c("sole","other(foot)"),c("heel","other(foot)")) 
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                      symbols = c("****", "***", "**", "*", "ns"))
  
  g <- ggboxplot(data,"Primary.Site",i,palette = 
                   c("#CC9933","#D28EFF","#00A087FF","#6633FF"),
                 fill = "Primary.Site") +
    # geom_jitter(data=data,aes(x=as.character(data[,clin]),y=data[,i],color=data[,clin]),width=0.15,size=0.6)+
    stat_compare_means(comparisons = my_comparisons,method="wilcox.test",)+
    stat_compare_means(label.y=0.1*ceiling(ymax)+1.5*ymax)+
    # stat_compare_means(comparisons=my_comparisons,method="wilcox.test")+ # #,aes(label=..p.signif..)
    labs(x=i,y=paste0("Percentage of ",i),color="")+
    # scale_fill_lancet()+
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(),
          axis.text.x=element_text(size=14,color="black"),
          axis.text.y=element_text(size=14,color="black"),
          axis.title.y=element_text(size=14),
          axis.title.x=element_blank(),
          legend.position = "right")
  
  ggsave(paste0("AM_boxplot_Site_",i,"_solo_clinical.pdf"),width = 6,height = 5)
}




########################四组比较
library(reshape2)
data2 <- tmb[,c('sample','Primary.Site','CNV(except CN2)','CNV(except CN2 CN3-5)')]
data3 <- melt(data2,id.vars = c("sample","Primary.Site"),variable.name = "CNV",value.name = "number")

burden <- ("CNV")
for(n in 1:length(burden)){
  i <- burden[n]
  clin <- "Primary.Site"
  ymax <- max(data3[,'number'])
  data <- data3[!is.na(data3[,clin]),]
  # data$cnv_count <- log(data$cnv_count)
  my_comparisons=list(c("hand","sole"),c("hand","heel"),c("hand","other(foot)"),c("sole","heel"),c("sole","other(foot)"),c("heel","other(foot)")) 
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                      symbols = c("****", "***", "**", "*", "ns"))
  
  g <- ggboxplot(data,"Primary.Site","number",palette = 
                   c("#CC9933","#D28EFF","#00A087FF","#6633FF"),
                 fill = "Primary.Site" ,facet.by = "CNV",) +
    # geom_jitter(data=data,aes(x=as.character(data[,clin]),y=data[,i],color=data[,clin]),width=0.15,size=0.6)+
    stat_compare_means(comparisons = my_comparisons,method="wilcox.test",)+
    stat_compare_means(label.y=1.8)+
    # stat_compare_means(comparisons=my_comparisons,method="wilcox.test")+ # #,aes(label=..p.signif..)
    labs(x=i,y="Percentage of CNV",color="")+
    # scale_fill_lancet()+
    theme(panel.background=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(),
          axis.text.x=element_text(size=14,color="black"),
          axis.text.y=element_text(size=14,color="black"),
          axis.title.y=element_text(size=14),
          axis.title.x=element_blank(),
          legend.position = "right")
  
  ggsave(paste0("AM_boxplot_Site_",i,"_clinical.pdf"),width = 8,height = 6)
}
