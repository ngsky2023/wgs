###add top genes on by frequency
rm(list=ls())
library(reshape2)
library(RColorBrewer)
library(plotrix)
library(scales)

####
mut<-c("Missense","In-frame indel","Frameshift indel","stopgain","Splice-site","RAS G12/G13/Q61","V600E")
col<-c("#2B406D","#7030A0",border=brewer.pal(9,"Set1")[1],"#5F7133","#008080","#99CC99","#ffafcc")
coll <- col
mut_col<-as.data.frame(cbind(mut,col))
##################################func####################################
Plot_func<-function(gene=gene,sample=sample,df=df,topgene=topgene,mut_col=mut_col,cex1=2,cex2=1,cex3=0.7,xdist=0,ydist=0.2){
  y<-topgene
  plot(x=nrow(sample)+1,y=topgene,xlim = c(0,nrow(sample)-1.5),ylim=c(0,topgene),cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
  for(i in as.character(gene$gene)){
    x<-0
    for(j in as.character(sample$Sample)){
      tmp<-df[which(df$gene==i & df$Sample==j),]
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
  df1<-df[,c("Sample","gene")]
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
  names(df) <- c("gene","Sample","TYPE","type","Freq")
  df<-df[,c("gene","Sample","TYPE","type")]
  df$type<-as.character(df$type)
  df$type[df$type=="frameshift insertion"]<-"Frameshift indel"
  df$type[df$type=="frameshift deletion"]<-"Frameshift indel"
  df$type[df$type=="frameshift substitution"]<-"Frameshift indel"
  df$type[df$type=="nonframeshift insertion"]<-"In-frame indel"
  df$type[df$type=="nonframeshift deletion"]<-"In-frame indel"
  df$type[df$type=="nonframeshift substitution"]<-"In-frame indel"
  df$type[df$type=="."]<-"Splice-site"
  df$type[df$type=="nonsynonymous SNV"]<-"Missense"
  df<-df[!duplicated(df),]
  # mut<-c("Missense","In-frame indel","Frameshift indel","stopgain","fusion","Splice-site")
  # col<-c("#2B406D","#7030A0",border=brewer.pal(9,"Set1")[1],"#5F7133",border=brewer.pal(9,"Set1")[8],"#009999")
  ###indicate mutation types by numbers
  df$class[df$type=="Missense"]<-1
  df$class[df$type=="Splice-site"]<-6
  df$class[df$type=="fusion"]<-5
  df$class[df$type=="stopgain"]<-4
  df$class[df$type=="In-frame indel"]<-2
  df$class[df$type=="Frameshift indel"]<-3
  df$class[df$type=="stoploss"]<-7
  df$class[df$type=="RAS G12/G13/Q61"]<-8
  df$class[df$type=="V600E"]<-9
  return(df)
}

###sample order in plot
Sample_order<-function(df,gene){
  ##for genes of each sample with more than one mutation type, retain the greatest "class"
  df<-df[order(df[,"gene"],df[,"Sample"],df[,"class"]),]
  df$filter<-paste(df$gene,df$Sample,sep="-")
  df_uniq_type<-df[!duplicated(df$filter,fromLast=TRUE),]
  df_uniq_type<-df_uniq_type[,c("gene","Sample","class","TYPE")]
  ##order each line by mutation type according to class number
  df_uniq_type_2<-dcast(df_uniq_type,gene~Sample)
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
  names(sample)<-"Sample"
  return(sample)
}

###frequency
Frequency<-function(df=df,gene=gene){
  freq <-NULL
  for(i in as.character(gene$gene)){
    tmp<-df[which(df$gene==i),]
    sample<-unique(tmp$Sample)
    num<-length(sample)
    freq <- c(freq,num)
  }
  return(freq)
}

#read data
df<-read.table("mutation_burden.xls", sep="\t", quote="", header=T,stringsAsFactors = F)
df<- df[,c(1,6,9,10)]
df <- df[order(df$cds_numtaion,decreasing = T),]
names(df)[1] <- c("Sample")
#WGS
wgd <- read.table("01.purity_WGD.txt",sep="\t", quote="", header=T,stringsAsFactors = F)
wgd <- wgd[,c(1,4)]
names(wgd) <- c("Sample","WGD")
#clinical

clinical<-read.table("clinical_55.txt", sep="\t", quote="", header=TRUE,na.strings = "")
clinical <- clinical[,c(1,2,5,4,10,7,9)]
names(clinical) <- c("Sample","name","Gender","Age","Site","Stage","Type")
clinical <- merge(clinical,wgd,by="Sample")

sbs_cluster <- read.table("SBS_cluster.xls",sep = "\t",quote = "",header = TRUE,na.strings = "")
names(sbs_cluster)[1] <- "Sample"
clinical <- merge(clinical,sbs_cluster,by="Sample")
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
# clinical$Stage[which(clinical$Stage=="Ib")] <- "I"
# clinical$Stage[which(clinical$Stage=="IIIC")] <- "III"
# clinical$Stage[which(clinical$Stage=="IIc")] <- "II"
# clinical$Stage[which(clinical$Stage=="IIa")] <- "II"
# clinical$Stage[which(clinical$Stage=="IIIB")] <- "III"
# clinical$Stage[which(clinical$Stage=="IIb")] <- "II"
#
clinical$Stage[which(clinical$Stage=="NA")] <-25
clinical$Stage[which(clinical$Stage=="I")] <- 8#
clinical$Stage[which(clinical$Stage=="II")] <- 9#
clinical$Stage[which(clinical$Stage=="III")] <- 10#
clinical$Stage[which(clinical$Stage=="IV")] <- 11#
#Gender
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
#sbs_cluster 38¿ªÊ¼
#'Polymerase','Unknown','BER_deficiency','Age','UV_light_exposure_indirect_effect','UV_light_exposure','APOBEK_activity','Age+UV_light_exposure','drug_mutagenesis'

clinical$sbs_cluster[which(clinical$sbs_cluster=="Polymerase")] <- 38
clinical$sbs_cluster[which(clinical$sbs_cluster=="Unknown")] <- 39
clinical$sbs_cluster[which(clinical$sbs_cluster=="BER_deficiency")] <- 40
clinical$sbs_cluster[which(clinical$sbs_cluster=="Age")] <- 41
clinical$sbs_cluster[which(clinical$sbs_cluster=="UV_light_exposure_indirect_effect")] <- 46
clinical$sbs_cluster[which(clinical$sbs_cluster=="UV_light_exposure")] <- 42
clinical$sbs_cluster[which(clinical$sbs_cluster=="APOBEK_activity")] <- 43
clinical$sbs_cluster[which(clinical$sbs_cluster=="Age+UV_light_exposure")] <- 44
clinical$sbs_cluster[which(clinical$sbs_cluster=="drug_mutagenesis")] <- 45

############


###merge
data <- merge(df,clinical,by="Sample")
data <- data[order(data$Age,data$Gender,data$Stage,decreasing=F),]
data <- data[order(data$Site,decreasing=F),]
# data <- data[order(data$Region,decreasing=F),]
data <- data[order(data$cds_numtaion,decreasing=T),]
sample_list <- as.data.frame(data$Sample)
names(sample_list) <- "Sample"
###mutation

#
gene <- read.delim2("gene_type.txt",header = T,stringsAsFactors = F)
snv <- read.delim2("acral_coding_snv_indel.xls",header = T,stringsAsFactors = F,sep = "\t")
snv <- snv[,c("SampleName","Gene.refGene","ExonicFunc.refGene")]
snv <- merge(snv,gene,by="Gene.refGene")
names(snv) <- c("gene","Sample","type","TYPE")
snv <- snv[,c("Sample","gene","TYPE","type")]
gene <- Pick_topgene(snv,topgene = length(unique(snv$gene)))
gene$gene <- as.character(gene$gene)
# gene$Freq[is.na(gene$Freq)] <- 0
snv <- Subdata(snv,gene)
for (i in c("gene","Sample","TYPE","type")) {
  snv[,i]<-as.character(snv[,i])
}
sample_order<-Sample_order(snv,gene)
sample_order$Sample<-as.character(sample_order$Sample)
sample_order<-merge(sample_order,sample_list,by="Sample",sort=FALSE,all=T)
F1 <- snv
F2 <- F1[which(F1$gene %in% c("NRAS","BRAF","NF1")),c(2,1)]
write.table(file = "sample_TWT.txt",F2,sep = "\t",quote = F,row.names = F)
#SMG

SMG <- F1[-which(F1$TYPE %in% c("ONCOGENE","TSG")),]
SMG <- subset(SMG,select = -TYPE)
# names(SMG) <- c("Sample","gene","type")
# SMG <- merge(sample_list,SMG,by="Sample",all=TRUE)
gene1<-Pick_topgene(SMG,topgene=length(unique(SMG$gene)))

# smg <- c("NRAS","BRAF","KIT","NF1","TYRP1","PTEN")
# smg <- as.data.frame(smg)
# names(smg) <- "gene"
# gene1<- merge(smg,gene1,by="gene",all = T)
# gene1 <- gene1[match(smg$gene,gene1$gene),]
# gene1$gene<-as.character(gene1$gene)
# gene1$Freq[is.na(gene1$Freq)] <- 0
#pick subset data from df by topgenes
# SMG<-Subdata(SMG,gene=gene1)
for (i in c("gene","Sample","type")) {
  SMG[,i]<-as.character(SMG[,i])
}
#
TSG <- F1[-which(F1$TYPE %in% c("ONCOGENE","SMG")),]
TSG <- subset(TSG,select = -TYPE)
gene2<-Pick_topgene(TSG,topgene=length(unique(TSG$gene)))
# tsg <- c("KMT2D","ATRX","NOTCH2","KMT2C","APC","ATM","TP53")
# tsg <- as.data.frame(tsg)
# names(tsg) <- "gene"
# gene2<- merge(tsg,gene2,by="gene",all = T)
# gene2[is.na(gene2)] <- 0
# gene2 <- gene2[match(tsg$gene,gene2$gene),]
# gene2$gene<-as.character(gene2$gene)
for (i in c("gene","Sample","type")) {
  TSG[,i]<-as.character(TSG[,i])
}

#ONCOGENE

ONCOGENE <- F1[-which(F1$TYPE %in% c("TSG","SMG")),]
ONCOGENE <- subset(ONCOGENE,select = -TYPE)
# names(ONCOGENE) <- c("Sample","gene","type")
# ONCOGENE <- merge(sample_list,ONCOGENE,by="Sample",all=TRUE)
gene3<-Pick_topgene(ONCOGENE,topgene=length(unique(ONCOGENE$gene)))
# oncogene <- c("ANK3","FGFR2","KRAS","MAP2K1","CTNNB1","HRAS","RAC1","RAF1")
# oncogene <- as.data.frame(oncogene)
# names(oncogene) <- "gene"
# gene3<- merge(oncogene,gene3,by="gene",all = T)
# gene3 <- gene3[match(oncogene$gene,gene3$gene),]
# gene3$gene<-as.character(gene3$gene)
# gene3[is.na(gene3)] <- 0
# gene3$gene<-as.character(gene3$gene)
#pick subset data from df by topgenes
# ONCOGENE<-Subdata(ONCOGENE,gene=gene1)
for (i in c("gene","Sample","type")) {
  ONCOGENE[,i]<-as.character(ONCOGENE[,i])
}
############
data <- data[match(sample_order$Sample,data$Sample),]
clin <- data[,c("Stage","Type","Age","Gender","sbs_cluster","WGD","Site")]
names(clin) <- c("Stage","Specimen type","Age","Gender","SBS_cluster","WGD","Primary site")
tmb <- data[,c(1,4,3)]
rownames(tmb) <- tmb[,1]
tmb <- tmb[,-1]
b <- t(tmb)
# colnames(tmb.t) <- tmb.t[1,]
# rownames(tmb.t) <- tmb.t[,1]
# tmb.t <- tmb.t[,-1]


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
library(scales)
#sbs_cluster
#"#00A087FF","#A9A9A9","#D28EFF","#590d22","#7B68EE","#4B0082","#FF0000","#2F4F4F","#7FFFD4"
#'Polymerase','Unknown','BER_deficiency','Age','UV_light_exposure_indirect_effect','UV_light_exposure','APOBEK_activity','Age+UV_light_exposure','drug_mutagenesis'

# cols <-  c("#ed1299", "#ddd53e", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92")
col2[class==38]<-"#00A087FF" #Polymerase
col2[class==39]<-"#A9A9A9" #Unknown
col2[class==40]<-"#D28EFF" #BER_deficiency
col2[class==41]<-"#590d22" #Age
col2[class==46]<-"#7B68EE"#UV_light_exposure_indirect_effect
col2[class==42]<-"#4B0082"#UV_light_exposure
col2[class==43]<-"#FF0000" #APOBEK_activity
col2[class==44]<-"#2F4F4F" #Age+UV_light_exposure
col2[class==45]<-"#7FFFD4" #drug_mutagenesis


##########################################################################
pdf("AM_SMG_clin.pdf",width=11.5, height=8.5)
layout(matrix(c(1:12), 6, 2, byrow = TRUE),  heights=c(2,0.6,1.4,1.4,1.8,1.5),widths = c(3,1))
#TMB
par(mar=c(1, 6.2, 2, 0.12), xpd=TRUE)
barplot(as.matrix(b), col = c("#00688B","#00EEEE"), border = "black",space=0.1, 
        names.arg=NULL,xaxt='n',yaxt='n',ylab="", xlab="", main="", cex.lab=1.4, cex.axis=1.7, font=2, las=1);
axis(2, at=c(0,100,200,300,400), labels=c(0,100,200,300,400), las=1, cex.axis=0.7,font = 1,line = -1);
mtext(side = 2, text = "Number of coding\n snv-indels", line = 1.9, cex=1,font = 1)
#
par(mar=c(0.3, 0.5, 0.1, 0.5), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,8, pch=15,  col=c("#00688B","#00EEEE"),legend=c(rownames(b)), title="", title.adj = 0, border=FALSE, cex=1.2,horiz=F,bty="n")

###SMG
par(mar=c(0.4, 7.2, 1, 0.5), xpd=TRUE)
Plot_func(gene=gene1,sample=sample_order,df=SMG,topgene=nrow(gene1), mut_col=mut_col, cex1=4.5, cex2=1.8, cex3=1.6, xdist=0.2,ydist=0.2)
axis(2, at=1:nrow(gene1), labels=rev(gene1$gene), line=-1, tick=FALSE, las=1, cex.axis=1.2,font=1)
mtext(side = 2, text = "SMG", line = 4.2, cex=1,font = 1)
axis(4, at=1:nrow(gene1), labels=paste(rev(round(gene1$Freq/nrow(sample_order)*100,2)),"%",sep=""), line=-1, tick=FALSE, las=1, cex.axis=1.2,font=1,mgp=c(3, 0, 0))
for(i in 0:nrow(sample_order)){
  lines(x=c(i-0.85,i-0.85),y=c(0.6,nrow(gene1)+0.6),col="black",lwd=0.5)
}
for(i in 0:nrow(gene3)){
  lines(x=c(-0.85,nrow(sample_order)-0.85),y=c(i+0.6,i+0.6),col="black",lwd=0.5)
}
#
par(mar=c(0.5, 0.5, 0.5, 0.5), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
#### axis(2, at=greater$at-0.2, labels=greater$greater, line=FALSE, tick=FALSE, las=1, cex.axis=1.7,font=1,col.axis="red",mgp=c(3, 2.5, 0))
#ONCOGENE
par(mar=c(0.4, 7.2, 0.9, 0.5), xpd=TRUE)
Plot_func(gene=gene3,sample=sample_order,df=ONCOGENE,topgene=nrow(gene3), mut_col=mut_col, cex1=4.5, cex2=1.8, cex3=1.6, xdist=0.2,ydist=0.2)
axis(2, at=1:nrow(gene3), labels=rev(gene3$gene), line=-1, tick=FALSE, las=1, cex.axis=1.2,font=1)
mtext(side = 2, text = "ONCOGENE", line = 4.2, cex=1,font = 1)
axis(4, at=1:nrow(gene3), labels=paste(rev(round(gene3$Freq/nrow(sample_order)*100,2)),"%",sep=""), line=-1, tick=FALSE, las=1, cex.axis=1.2,font=1,mgp=c(3, 0, 0))
for(i in 0:nrow(sample_order)){
  lines(x=c(i-0.85,i-0.85),y=c(0.6,nrow(gene3)+0.6),col="black",lwd=0.5)
}
for(i in 0:nrow(gene3)){
  lines(x=c(-0.85,nrow(sample_order)-0.85),y=c(i+0.6,i+0.6),col="black",lwd=0.5)
}
#
par(mar=c(0, 1, 1, 0.5), xpd=TRUE)
plot(x=5,y=7,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=coll,legend=mut, title="", title.adj = 0, border=FALSE, cex=1 ,horiz=F,ncol = 1,bty="n")

#### axis(2, at=greater$at-0.2, labels=greater$greater, line=FALSE, tick=FALSE, las=1, cex.axis=1.7,font=1,col.axis="red",mgp=c(3, 2.5, 0))

# for(i in 0:nrow(gene)){
#   lines(x=c(-1,nrow(sample_order)),y=c(i+0.5,i+0.5),col="white",lwd=0.5)
# }
#TSG
par(mar=c(0.4, 7.2, 0.9, 0.5), xpd=TRUE)
Plot_func(gene=gene2,sample=sample_order,df=TSG,topgene=nrow(gene2), mut_col=mut_col, cex1=4.5, cex2=1.8, cex3=1.6, xdist=0.2,ydist=0.2)
axis(2, at=1:nrow(gene2), labels=rev(gene2$gene), line=-1, tick=FALSE, las=1, cex.axis=1.2,font=1)
mtext(side = 2, text = "TSG", line = 4.2, cex=1,font = 1)
axis(4, at=1:nrow(gene2), labels=paste(rev(round(gene2$Freq/nrow(sample_order)*100,2)),"%",sep=""), line=-1, tick=FALSE, las=1, cex.axis=1.2,font=1,mgp=c(3, 0, 0))
for(i in 0:nrow(sample_order)){
  lines(x=c(i-0.85,i-0.85),y=c(0.6,nrow(gene2)+0.6),col="black",lwd=0.5)
}
for(i in 0:nrow(gene2)){
  lines(x=c(-0.85,nrow(sample_order)-0.85),y=c(i+0.6,i+0.6),col="black",lwd=0.5)
}
#
par(mar=c(0, 0.5, 0, 0.5), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
# legend(3,11, pch=15,  col=coll,legend=mut, title="", title.adj = 0, border=FALSE, cex=1,horiz=F,bty="n")

#clin
par(mar=c(0.4, 8.3, 0.1, 2), xpd=TRUE)
color2D.matplot(class, cellcolors=col2, border = NA, main="", xlab="", ylab="", axes=FALSE) 
axis(2,at=0.5:(nrow(class)-0.5),labels = rev(rownames(class)),las=2,cex.axis=1.2,font=1,line = F,tick = F)
for(i in (0:(ncol(class)))){
  segments(i,0,i,nrow(class),col = "black",lwd=0)
}
for(i in (0:(nrow(class)))){
  segments(0,i,ncol(class),i,col = "black",lwd=1)
}
# color2D.matplot(plot_df, cellcolors=col, border = "white", main="", xlab="", ylab="", axes=FALSE)
axis(1, at=c(1:length(data$name))-0.5,labels=data$name, line=FALSE, tick=FALSE, las=2, cex.axis=1.2,font=1,mgp=c(0, 1,1))

par(mar=c(0, 0, 0, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,12, pch=15,  col=c("#EDE0D4","#E6CCB2","#DDB892","#9C6644","#d3d3d3"),legend=c("I","II","III","IV","NA"), title="", title.adj = 0, border=FALSE, cex=1,horiz=T,bty="n")
legend(3,11.5, pch=15,  col=c("#03827f","#f94144"),legend=c("Primary","Metastasis"),title="", title.adj = 0, border=FALSE,  cex=1,horiz=T,bty="n")
legend(3,11, pch=15,  col=c("#ff8fa3","#590d22","#d3d3d3"),legend=c("<58",">=58","NA"), title="",title.adj = 0, border=FALSE, cex=1,horiz=T,bty="n")
legend(3,10.5, pch=15,  col=c("#fbb1bd","#bbd0ff"),legend=c("F","M"),title="", title.adj = 0, border=FALSE,  cex=1,horiz=T,bty="n")
legend(3,10, pch=15,  col=c("#00A087FF","#A9A9A9","#D28EFF","#590d22","#7B68EE","#4B0082","#FF0000","#2F4F4F","#7FFFD4"),
       legend=c('Polymerase','Unknown','BER_deficiency','Age','UV_light_exposure_indirect_effect','UV_light_exposure','APOBEK_activity','Age+UV_light_exposure','drug_mutagenesis'),
       title="", title.adj = 0, border=FALSE,  cex=1,ncol = 2,bty="n")
legend(3,7.3, pch=15, col=c("#555555","#99CC99"),legend=c("non-WGD","WGD"), title="",title.adj = 0, border=FALSE, cex=1,horiz=F,ncol = 2,bty="n")
legend(3,6.8, pch=15, col=c("#CC9933","#D28EFF","#00A087FF","#6633FF"),legend=c("hand","heel","sole","other(foot)"), title="",title.adj = 0, border=FALSE, cex=1,horiz=F,ncol = 2,bty="n")


dev.off()
#################################
#"#00A087FF","#A9A9A9","#D28EFF","#590d22","#7B68EE","#4B0082","#FF0000","#2F4F4F","#ed1299"
#'Polymerase','Unknown','BER_deficiency','Age','UV_light_exposure_indirect_effect','UV_light_exposure','APOBEK_activity','Age+UV_light_exposure','drug_mutagenesis'

