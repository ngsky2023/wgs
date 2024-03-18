###add top genes on by frequency
rm(list=ls())
library(reshape2)
library(RColorBrewer)
library(plotrix)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(grid)
library(ggsignif)
library(ggsci)
library(dplyr)
library(pheatmap)
library(gridExtra)
library(scales)
##################################func####################################
topgene<-16
show_col(c("#F7B671","#EE7E32","#EF8E8F","#D8322F","#A1CF89","#38A245","#9BC6D4","#2B73A9","#B9A5C5","#5B448B","#848683"))
###set block color
show_col(c("#F7B671","#EE7E32","#EF8E8F","#D8322F","#A1CF89","#38A245","#9BC6D4","#2B73A9","#B9A5C5","#5B448B","#848683"))
###set block color
mut<-c("Translocation","Inversion","Deletion","Duplication","amp","CN>=6","loss","del")
col<-c("#9BC6D4","#F7B671","#5B448B","#003C9D","#FFCCCC","#D8322F","#66FFCC","#339900")
coll <- col
mut_col<-as.data.frame(cbind(mut,col))
##################################func####################################
Plot_func<-function(gene=gene,sample=sample,df=df,topgene=topgene,mut_col=mut_col,cex1=2,cex2=1,cex3=0.7,cex4=0.4,xdist=0,ydist=0.2){
  y<-topgene
  plot(x=nrow(sample)+1,y=topgene,xlim = c(0,nrow(sample)-1.5),ylim=c(0,topgene),cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
  for(i in as.character(gene$gene)){
    x<-0
    for(j in as.character(sample$sample)){
      tmp<-df[which(df$gene==i & df$sample==j),]
      tmp<-tmp[order(tmp[,"class"],decreasing=TRUE),]
      if(nrow(tmp)==0){
        points(x=x,y=y,pch=15,col="#EDEDED",cex=cex1)
      }else if(nrow(tmp)>0 & nrow(tmp)<=2){
        for(k in 1:nrow(tmp)){
          if(k==1){
            points(x=x,y=y,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[1,"type"]),"col"]),cex=cex1)
          }else if(k==2){
            points(x=x-xdist,y=y,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[2,"type"]),"col"]),cex=cex2)
          }
        }
      }else if(nrow(tmp)==3){
        for(k in 1:nrow(tmp)){
          if(k==1){
            points(x=x,y=y,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[1,"type"]),"col"]),cex=cex1)
          }else if(k==2){
            points(x=x-xdist,y=y+0.2,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[2,"type"]),"col"]),cex=cex3)
          }else if(k==3){
            points(x=x-xdist,y=y-0.2,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[3,"type"]),"col"]),cex=cex3)
          }
        }
      }else if(nrow(tmp)==4){
        for(k in 1:nrow(tmp)){
          if(k==1){
            points(x=x,y=y,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[1,"type"]),"col"]),cex=cex1)
          }else if(k==2){
            points(x=x-xdist,y=y+0.5,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[2,"type"]),"col"]),cex=cex4)
          }else if(k==3){
            points(x=x-xdist,y=y+0.2,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[3,"type"]),"col"]),cex=cex4)
          }else if(k==4){
            points(x=x-xdist,y=y-0.2,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[4,"type"]),"col"]),cex=cex4)
          }
        }
      }
      x<-x+1
    }
    #print(nrow(sample)+1)
    #print(y)
    points(x=nrow(sample),y=y,pch=15,col="white",cex=cex1)
    y<-y-1
  }
  for(i in 0:nrow(sample)){
    points(x=i,y=0,pch=15,col="white",cex=cex1)
  }
}
#
Pick_topgene<-function(df,topgene=topgene){
  df1<-df[,c("sample","gene")]
  df1<-df1[!duplicated(df1),]
  gene<-as.data.frame(table(df1$gene))
  gene<-gene[order(gene[,2],decreasing = TRUE),]
  print(head(gene,n=50L))
  gene<-gene[c(1:topgene),]
  colnames(gene)<-c("gene","Freq")
  return(gene)
}

Subdata<-function(df=df,gene=gene){
  df<-merge(df,gene,by="gene",sort=FALSE)
  names(df) <- c("gene","sample","type","Freq")
  df<-df[,c("gene","sample","type")]
  df$type<-as.character(df$type)
  df<-df[!duplicated(df),]
  ###indicate mutation types by numbers
  # df$class[df$type=="chromothripsis"]<-1
  df$type[df$type=="BND"]<-"Translocation"
  df$type[df$type=="INV"]<-'Inversion'
  df$type[df$type=="DEL"]<-'Deletion'
  df$type[df$type=="DUP"]<-'Duplication'
  df$class[df$type=="amp"] <- 1
  df$class[df$type=="CN>=6"]<-2
  df$class[df$type=="loss"]<-3
  df$class[df$type=="del"]<-4
  df$class[df$type=="Translocation"]<-8
  df$class[df$type=="Inversion"]<-7
  df$class[df$type=="Deletion"]<-6
  df$class[df$type=="Duplication"]<-5
  
  return(df)
}

###sample order in plot
Sample_order<-function(df,gene){
  ##for genes of each sample with more than one mutation type, retain the greatest "class"
  df<-df[order(df[,"gene"],df[,"sample"],df[,"class"]),]
  df$filter<-paste(df$gene,df$sample,sep="-")
  df_uniq_type<-df[!duplicated(df$filter,fromLast=TRUE),]
  df_uniq_type<-df_uniq_type[,c("gene","sample","class")]
  ##order each line by mutation type according to class number
  df_uniq_type_2<-dcast(df_uniq_type,gene~sample)
  df_uniq_type_2<-merge(gene,df_uniq_type_2,by="gene",sort = FALSE)
  rownames(df_uniq_type_2)<-df_uniq_type_2[,1]
  df_uniq_type_2<-df_uniq_type_2[,-c(1,2)]
  df_uniq_type_2[is.na(df_uniq_type_2)]<-0
  order <- "df_uniq_type_2 <- df_uniq_type_2[,order("
  for(i in 1:nrow(df_uniq_type_2)){
    order <- paste(order,"df_uniq_type_2[",i,",],",sep = "")
  }   
  order <- paste(order,"decreasing = TRUE)]")
  df_uniq_type_2 <- eval(parse(text=order))
  sample<-as.data.frame(colnames(df_uniq_type_2))
  names(sample)<-"sample"
  return(sample)
}

###frequency
Frequency<-function(df=df,gene=gene){
  freq <-NULL
  for(i in as.character(gene$gene)){
    tmp<-df[which(df$gene==i),]
    sample<-unique(tmp$sample)
    num<-length(sample)
    freq <- c(freq,num)
  }
  return(freq)
}


#read data
driver_gene <- read.delim2("snv_sv_cnv_driver.txt",header=F,stringsAsFactors = F)
colnames(driver_gene) <- "gene"
#
sv <- read.table("AC_snv_sv_cnv_driver_sv_result_20220221.xls",header=T,stringsAsFactors = F)
sv <- sv[which(sv$gene %in% driver_gene$gene),]
names(sv)[3] <- "mutation"
#
cnv <- read.table("AC_cnv_all_gene_20220511.xls",header=T,stringsAsFactors = F)
cnv <- cnv[,c(1,8,9)]
cnv <- cnv[which(cnv$gene %in% driver_gene$gene),]
names(cnv)[3] <- "mutation"
#
F1 <- sv
F1$TYPE <- "driver"
names(F1) <- c("sample","gene","type","TYPE")

gene_freq <- Pick_topgene(df=F1,topgene = topgene)
F1 <- merge(F1,gene_freq,by="gene")
gene_freq$gene <- as.character(gene_freq$gene)
F1 <- Subdata(F1,gene_freq)
#
F2 <- cnv
F2$TYPE <- "driver"
names(F2) <- c("sample","gene","type","TYPE")

gene_freq2 <- Pick_topgene(df=F2,topgene = length(unique(F2$gene)))
F2 <- merge(F2,gene_freq2,by="gene")
gene_freq2$gene <- as.character(gene_freq2$gene)
F2 <- Subdata(F2,gene_freq2)
#
#WGS
chromo <- read.delim2("ac_chromothripsis_class.xls", sep="\t", quote="", header=T,stringsAsFactors = F)
chromo$chrom[chromo$chrom=="X"] <- 23
# chromo <- chromo[-which(chromo$class %in% c("non","unknown")),]
chromothripsis <- unique(chromo[,c(2,1)])

chromothripsis_hight <- chromothripsis[grep("high",chromothripsis$class),]
chromothripsis_hight$class <- "high"
chromothripsis_low <- chromothripsis[grep("low",chromothripsis$class),]
chromothripsis_low$class <- "low"
chromothripsis <- unique(rbind(chromothripsis_hight,chromothripsis_low))


# chromothripsis$class2[chromothripsis$class %in% c("high1;high2","high1","high2","high3")] <- "high"
# chromothripsis <- chromothripsis[which(chromothripsis$class %in% c("high1")),]
b <- unique(chromothripsis)
b_data <- lapply(unique(b$sample),function(x) b[match(x,b$sample),] )
b_data <- unique(do.call(rbind,b_data))
names(b_data) <- c("sample","chromothripsis")
##

#WGS
wgd <- read.table("01.purity_WGD.txt",sep="\t", quote="", header=T,stringsAsFactors = F)
wgd <- wgd[,c(1,4)]
names(wgd) <- c("sample","WGD")
#
# #sbs_cluster
# sbs_cluster <- read.table("SBS_cluster.xls",sep = "\t",quote = "",header = TRUE,na.strings = "")
# names(sbs_cluster)[1] <- "sample"
#RS_cluster
RS_cluster <- read.table("RS_cluster.xls",sep = "\t",quote = "",header = TRUE,na.strings = "")
names(RS_cluster)[1] <- "sample"
#
kata_type <-read.table("ac_kata_sample_seq.xls", sep="\t", quote="", header=TRUE,na.strings = "")
names(kata_type) <- c("sample","kataegis")
#sample driver
sample_twt <- read.table("sample_TWT.txt",header = T,stringsAsFactors = F,sep = "\t")
names(sample_twt) <- c("sample","TWT")
#clinical


clinical<-read.table("clinical_55.txt", sep="\t", quote="", header=TRUE,na.strings = "")
clinical <- clinical[,c(1,2,5,4,10,7,9)]
names(clinical) <- c("sample","name","Gender","Age","Site","Stage","Type")
clinical <- merge(clinical,wgd,by="sample")
clinical <- merge(clinical,kata_type,by="sample")
clinical <- merge(clinical,b_data,by="sample",all=TRUE)
clinical <- merge(clinical,RS_cluster,by="sample")#clinical
clinical <- merge(clinical,sample_twt,by="sample",all=T)
clinical$TWT[is.na(clinical$TWT)] <- "TWT"

#
clinical$Type[which(clinical$Type=="metastasis")] <-20
clinical$Type[which(clinical$Type=="primary")] <- 19
clinical$Type[which(clinical$Type=="NA")] <- 25



#Age

age <- median(na.omit(clinical$Age) )
print(age)
clinical$age[which(clinical$Age >=58)] <- ">58"#"70-79y"
clinical$age[which(clinical$Age <58)] <- "<58"
clinical$age[which(clinical$Age =="NA")] <- "NA"
#
clinical$Age[which(clinical$age =="NA")] <- 25#"70-79y"
clinical$Age[which(clinical$age %in% c(">58"))] <- 7
clinical$Age[which(clinical$age %in% c("<58"))] <- 4
#Stage
#Stage
clinical$Stage[which(clinical$Stage=="Ib")] <- "I"
clinical$Stage[which(clinical$Stage=="IIIC")] <- "III"
clinical$Stage[which(clinical$Stage=="IIc")] <- "II"
clinical$Stage[which(clinical$Stage=="IIa")] <- "II"
clinical$Stage[which(clinical$Stage=="IIIB")] <- "III"
clinical$Stage[which(clinical$Stage=="IIb")] <- "II"
#
clinical$Stage[which(clinical$Stage=="NA")] <-25
clinical$Stage[which(clinical$Stage=="I")] <- 8#
clinical$Stage[which(clinical$Stage=="II")] <- 9#
clinical$Stage[which(clinical$Stage=="III")] <- 10#
clinical$Stage[which(clinical$Stage=="IV")] <- 11#

clinical$Gender[which(clinical$Gender=="F")] <- 13
clinical$Gender[which(clinical$Gender=="M")] <- 14
clinical$Gender[which(clinical$Gender=="NA")] <- 25
# clinical$Gender <- factor(clinical$Gender,levels = c("F","M"))
#Site
clinical$Site[which(clinical$Site=="NA")] <- 25

clinical$Site[which(clinical$Site=="hand")] <- 15
clinical$Site[which(clinical$Site=="heel")] <- 16
clinical$Site[which(clinical$Site=="sole")] <- 17
clinical$Site[which(clinical$Site=="other")] <- 18
#wgd

clinical$WGD[which(clinical$WGD==FALSE)] <- 23
clinical$WGD[which(clinical$WGD==TRUE)] <- 24
#
clinical$kata[which(clinical$kataegis > 20) ] <- 31
clinical$kata[which(clinical$kataegis>10 & clinical$kataegis<=20)] <- 30
clinical$kata[which(clinical$kataegis <=10)] <- 28
#
clinical$chromothripsis[is.na(clinical$chromothripsis)] <-25 #没有染色体碎裂

clinical$chromothripsis[which(clinical$chromothripsis=="low")] <- 32
clinical$chromothripsis[which(clinical$chromothripsis=="high")] <- 33


#APOBEC
# clinical$APOBEC[clinical$APOBEC=="Yes"] <- 34
# clinical$APOBEC[clinical$APOBEC=="No"] <- 35
#
clinical$TWT[which(clinical$TWT=="BRAF")] <- 34
clinical$TWT[which(clinical$TWT=="NRAS")] <- 35
clinical$TWT[which(clinical$TWT=="NF1")] <- 36
clinical$TWT[which(clinical$TWT=="TWT")] <- 37


#sbs_cluster
# clinical$sbs_cluster[which(clinical$sbs_cluster=="Polymerase")] <- 38
# clinical$sbs_cluster[which(clinical$sbs_cluster=="Unknown")] <- 39
# clinical$sbs_cluster[which(clinical$sbs_cluster=="BER_deficiency")] <- 40
# clinical$sbs_cluster[which(clinical$sbs_cluster=="Age")] <- 41
# clinical$sbs_cluster[which(clinical$sbs_cluster=="UV_light_exposure_indirect_effect")] <- 46
# clinical$sbs_cluster[which(clinical$sbs_cluster=="UV_light_exposure")] <- 42
# clinical$sbs_cluster[which(clinical$sbs_cluster=="APOBEK_activity")] <- 43
# clinical$sbs_cluster[which(clinical$sbs_cluster=="Age+UV_light_exposure")] <- 44
# clinical$sbs_cluster[which(clinical$sbs_cluster=="drug_mutagenesis")] <- 45
#
#RS_cluter
clinical$RS_cluster[which(clinical$RS_cluster=="A")] <- 47
clinical$RS_cluster[which(clinical$RS_cluster=="B")] <- 48
clinical$RS_cluster[which(clinical$RS_cluster=="C")] <- 49
clinical$RS_cluster[which(clinical$RS_cluster=="D")] <- 50
clinical$RS_cluster[which(clinical$RS_cluster=="E")] <- 51
clinical$RS_cluster[which(clinical$RS_cluster=="F")] <- 52



#
############


###merge
# data <- merge(df,clinical,by="Sample")
data <- clinical
data <- data[order(data$Age,data$Gender,data$Stage,decreasing=F),]
data <- data[order(data$Site,decreasing=F),]
# data <- data[order(data$cds_numtaion,decreasing=T),]
############################RS cluster##############################################

RS_sig <- read.delim("sv_5_samples.xls",header = T,stringsAsFactors = F)
colnames(RS_sig)<-paste0("RS",c(1:(ncol(RS_sig))))
names(RS_sig)[1] <- "sample"

sig4 <-RS_sig
# sig4 <- sig4[match(sample_order1$sample,sig4$sample),]
rownames(sig4) <- sig4[,1]
sig4 <- sig4[,-1]
sig4 <- as.matrix(t(sig4))
result8 <- pheatmap(sig4,
                    scale="column",cluster_rows = TRUE,
                    cutree_cols=1,cutree_rows = 1,show_colnames = TRUE,
                    color = colorRampPalette(colors = c("#fad2e1","white","#EE7E32"))(20))#"#fff0f3","#c9184a"
sample_order7 <- as.data.frame(colnames(sig4)[result8$tree_col$order])
cutree(result8$tree_col,k=5)
colnames(sample_order7) <- "sample"
sample_order <- sample_order7
ggsave("RS_noname.pdf",result8,height = 5,width = 10)
#
p <- as.data.frame(cutree(result8$tree_col,k=6))
colnames(p) <- c("cluster")
p$sample <- rownames(p)
p <- p[match(sample_order$sample,p$sample),]
cluster <- p
print(cluster$cluster)
cluster$cluster[which(cluster$cluster %in% c(2))] <- "A"
cluster$cluster[which(cluster$cluster %in% c(4))] <- "B"
cluster$cluster[which(cluster$cluster %in% c(6))] <- "C"
cluster$cluster[which(cluster$cluster %in% c(1))] <- "D"
cluster$cluster[which(cluster$cluster %in% c(5))] <- "E"
cluster$cluster[which(cluster$cluster %in% c(3))] <- "F"
RS_cluster <- as.data.frame(cluster[,c("sample","cluster")])
names(RS_cluster) <- c("sample","RS_cluster")
write.table(RS_cluster,file = "RS_cluster.xls",row.names = F,quote = F,sep = "\t")
print(RS_cluster)

annotation_col = as.data.frame(RS_cluster[,c("RS_cluster")])  #增加Time，CellType分组信息
names(annotation_col) <- c("RS_cluster")
ann_colors = list(RS_cluster=c(A="#E64B35FF",B="#F39B7FFF",C="#4DBBD5FF",D="#00A087FF",E="#03827f",F="#f94144"))
rownames(annotation_col) <- RS_cluster$sample
resultRs <- pheatmap(sig4,annotation_col = annotation_col,annotation_colors = ann_colors,
                     scale="column",cluster_rows = TRUE,
                     cutree_cols=1,cutree_rows = 1,show_colnames = TRUE)
ggsave("RS_noname.pdf",resultRs,height = 5,width = 10)
#####################

data <- data[match(sample_order$sample,data$sample),]

############
data <- data[match(sample_order$sample,data$sample),]
clin <- data[,c("Stage","Type","Age","Gender","WGD","kata","chromothripsis","Site","TWT","RS_cluster")]
names(clin) <- c("Stage","Specimen type","Age","Gender","WGD","kataegeis","chromothripsis","Primary site","Driver","RS_cluster")
############
for(i in 1:ncol(clin)){
  clin[,i]<-as.numeric(clin[,i])
}

class<-t(clin)
col2<-matrix(rep("#000000",length(class)),nrow=nrow(class))
col2[class==1]<-"#5371A3"
#Age"#fff0f3","#ffccd5","#ff8fa3","#ff4d6d","#c9184a","#590d22"
col2[class==2]<-"#fff0f3"
col2[class==3]<-"#ffccd5"
col2[class==4]<-"#ff8fa3"
col2[class==5]<-"#ff4d6d"
col2[class==6]<-"#c9184a"
col2[class==7]<-"#590d22"
#Stage"#EDE0D4","#E6CCB2","#DDB892","#9C6644"
col2[class==8]<-"#EDE0D4"
col2[class==9]<-"#E6CCB2"
col2[class==10]<-"#DDB892"
col2[class==11]<-"#9C6644"
#
col2[class==12]<-"lightgrey"
#sex"#fbb1bd","#bbd0ff"
col2[class==13]<-"#fbb1bd"
col2[class==14]<-"#bbd0ff"
#site
shows=c("#CC9933","#D28EFF","#00A087FF","#6633FF")
show_col(shows)
col2[class==15]<-"#CC9933" #hand
col2[class==16]<-"#D28EFF"#heel
col2[class==17]<-"#00A087FF"#sole
col2[class==18]<-"#6633FF"#other
#type"#03827f","#f94144"
col2[class==19]<-"#03827f"
col2[class==20]<-"#f94144"
#Region "#48D1CC","#996666"
col2[class==21]<-"#48D1CC"
col2[class==22]<-"#996666"
#WGD "#555555","#99CC99"
col2[class==23]<-"#555555"
col2[class==24]<-"#99CC99"
#
col2[class==25]<-"#d3d3d3"
cols <-  c("#ed1299", "#ddd53e", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92")
#kata
#c("<10","10-20","20-30","30-40","40-100",">100")
col2[class==26]<-"#fff0f3"
col2[class==27]<-"#ffccd5"
col2[class==28]<-"#ff8fa3"
col2[class==29]<-"#ff4d6d"
col2[class==30]<-"#c9184a"
col2[class==31]<-"#590d22"


#chromothri
#"high","low","Non"
#"#ff8fa3","#c9184a","#d3d3d3"
col2[class==32]<-"#ff8fa3" #low
col2[class==33]<-"#c9184a" #hig

#twt
#"#FF33FF","#990066","#99FF99","#CC9999"
col2[class==34]<-"#FF33FF"
col2[class==35]<-"#990066"
col2[class==36]<-"#99FF99"
col2[class==37]<-"#CC9999"
# #sbs_cluster
# col2[class==38]<-"#00A087FF" #Polymerase
# col2[class==39]<-"#A9A9A9" #Unknown
# col2[class==40]<-"#D28EFF" #BER_deficiency
# col2[class==41]<-"#590d22" #Age
# col2[class==46]<-"#7B68EE"#UV_light_exposure_indirect_effect
# col2[class==42]<-"#4B0082"#UV_light_exposure
# col2[class==43]<-"#FF0000" #APOBEK_activity
# col2[class==44]<-"#2F4F4F" #Age+UV_light_exposure
# col2[class==45]<-"#7FFFD4" #drug_mutagenesis

#RS_cluster
# "A","B","C","D","E","F"
#"#E64B35FF","#F39B7FFF","#4DBBD5FF","#00A087FF","#03827f","#f94144"
col2[class==47]<-"#E64B35FF"
col2[class==48]<-"#F39B7FFF"
col2[class==49]<-"#4DBBD5FF"
col2[class==50]<-"#00A087FF"
col2[class==51]<-"#03827f"
col2[class==52]<-"#f94144"
############RScluster end
pdf("AM_china_sv_complex_clin.pdf",width=15, height=17)

layout(matrix(c(1:8), 4, 2, byrow = TRUE),  heights=c(0.9,3,0.4,0.5),widths = c(3,1))

###sv
par(mar=c(0, 6.5, 2, 1), xpd=TRUE)
Plot_func(gene=gene_freq,sample=sample_order,df=F1,topgene=nrow(gene_freq), mut_col=mut_col, cex1=4.5, cex2=1.8, cex3=1.6, xdist=0.2,ydist=0.2)
axis(2, at=1:nrow(gene_freq)+0.2, labels=rev(gene_freq$gene), line=-1, tick=FALSE, las=1, cex.axis=1,font=1)
# axis(1, at=c(1:length(data$Sample))-1,labels=data$Sample, line=-2, tick=FALSE, las=2, cex.axis=1,font=1,mgp=c(0, 1,1))
mtext(side = 2, text = "sv", line = 5, cex=1,font = 1)
axis(4, at=1:nrow(gene_freq), labels=paste(rev(round(gene_freq$Freq/nrow(sample_order)*100,2)),"%",sep=""), line=-1, tick=FALSE, las=1, cex.axis=1,font=1,mgp=c(3, 0, 0))
for(i in 0:nrow(sample_order)){
  lines(x=c(i-0.65,i-0.65),y=c(0.7,nrow(gene_freq)+0.7),col="black",lwd=0.5)
}
for(i in 0:nrow(gene_freq)){
  lines(x=c(-0.65,nrow(sample_order)-0.65),y=c(i+0.7,i+0.7),col="black",lwd=0.5)
}

#
par(mar=c(0.1, 1, 0, 0.5), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,8, pch=15,  col=mut_col$col[1:4],legend=mut_col$mut[1:4], title="", title.adj = 0, border=FALSE, cex=1.2,horiz=F,bty="n")
##############cnv
par(mar=c(0, 6.5, 0, 1), xpd=TRUE)
Plot_func(gene=gene_freq2 , sample=sample_order,df=F2,topgene=nrow(gene_freq2), mut_col=mut_col, cex1=4.5, cex2=1.8, cex3=1.6, xdist=0.2,ydist=0.2)
axis(2, at=1:nrow(gene_freq2), labels=rev(gene_freq2$gene), line=-1, tick=FALSE, las=1, cex.axis=1,font=1)
mtext(side = 2, text = "cnv", line = 5, cex=1,font = 1)
axis(4, at=1:nrow(gene_freq2), labels=paste(rev(round(gene_freq2$Freq/nrow(sample_order)*100,2)),"%",sep=""), line=0, tick=FALSE, las=1, cex.axis=1,font=1,mgp=c(3, 0, 0))
for(i in 0:nrow(sample_order)){
  lines(x=c(i-0.65,i-0.65),y=c(0.7,nrow(gene_freq2)+0.7),col="black",lwd=0.5)
}
for(i in 0:nrow(gene_freq2)){
  lines(x=c(-0.65,nrow(sample_order)-0.65),y=c(i+0.7,i+0.7),col="black",lwd=0.5)
}
#leng
par(mar=c(0.1, 1, 0, 0.5), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,8, pch=15,  col=mut_col$col[5:8],legend=mut_col$mut[5:8], title="", title.adj = 0, border=FALSE, cex=1.2,horiz=F,bty="n")

#clin
par(mar=c(0.1,8.5, 0, 2.8), xpd=TRUE)
color2D.matplot(class, cellcolors=col2, border = NA, main="", xlab="", ylab="", axes=FALSE)
axis(1,at=0.5:(nrow(data)-0.5),labels =data$name,las=2,cex.axis=1,font=1,line = F,tick = F)

axis(2,at=0.5:(nrow(class)-0.5),labels = rev(rownames(class)),las=2,cex.axis=1,font=1,line = F,tick = F)
for(i in (0:(ncol(class)))){
  segments(i,0,i,nrow(class),col = "black",lwd=0)
}
for(i in (0:(nrow(class)))){
  segments(0,i,ncol(class),i,col = "black",lwd=1)
}
par(mar=c(0, 0, 0, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,12, pch=15,  col=c("#EDE0D4","#E6CCB2","#DDB892","#9C6644","#d3d3d3"),legend=c("I","II","III","IV","NA"), title="", title.adj = 0, border=FALSE, cex=0.7,horiz=T,bty="n")
legend(3,11.5, pch=15,  col=c("#03827f","#f94144"),legend=c("Primary","Metastasis"),title="", title.adj = 0, border=FALSE,  cex=0.7,horiz=T,bty="n")
legend(3,11, pch=15,  col=c("#ff8fa3","#590d22","#d3d3d3"),legend=c("<58",">=58","NA"), title="",title.adj = 0, border=FALSE, cex=0.7,horiz=T,bty="n")
legend(3,10.5, pch=15,  col=c("#fbb1bd","#bbd0ff"),legend=c("F","M"),title="", title.adj = 0, border=FALSE,  cex=0.7,horiz=T,bty="n")
legend(3,10, pch=15, col=c("#555555","#99CC99"),legend=c("non-WGD","WGD"), title="",title.adj = 0, border=FALSE, cex=0.7,horiz=F,ncol = 2,bty="n")
legend(3,9, pch=15, col=c("#ff8fa3","#c9184a","#590d22"),legend=c("<10","10-20",">20"), title="kataegeis",title.adj = 0, border=FALSE, cex=0.7,horiz=F,ncol = 3,bty="n")
legend(3,8, pch=15, col=c("#ff8fa3","#c9184a","#d3d3d3"),legend=c("low","high","Non"), title="chromothripsis",title.adj = 0, border=FALSE, cex=0.7,horiz=T,bty="n")
legend(3,7.5, pch=15, col=c("#CC9933","#D28EFF","#00A087FF","#6633FF"),legend=c("hand","heel","sole","other(foot)"), title="",title.adj = 0, border=FALSE, cex=0.7,horiz=F,ncol = 2,bty="n")
legend(3,6.5, pch=15, col=c("#FF33FF","#990066","#99FF99","#CC9999"),legend=c("BRAF","NRAS","NF1","Triple WT"), title="",title.adj = 0,text.width = 0.5, border=FALSE, cex=0.7,horiz=T,bty="n")
legend(3,5.8, pch=15,  col=c("#E64B35FF","#F39B7FFF","#4DBBD5FF","#00A087FF","#03827f","#f94144"),
       legend=c("A","B","C","D","E","F"),title="", title.adj = 0, border=FALSE,  cex=0.7,horiz = T,bty="n")

dev.off()
#################################

