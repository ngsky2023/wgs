###add top genes on by frequency
rm(list=ls())
library(reshape2)
library(RColorBrewer)
library(plotrix)
library(ggsci)
library(scales)
##################################func####################################
df<-read.table("mutation_burden.xls", sep="\t", quote="", header=T,stringsAsFactors = F)
df<- df[,c(1,5,6,7)]
df <- df[order(df$perMb,decreasing = T),]

sbs_sig <- read.table("Decomposed_Solution_Activities_SBS96.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(sbs_sig)[1] <- "sample"

##UV_light_exposure(SBS7a,7b,7c,7d)
# sbs_sig$UV_light_exposure <- sbs_sig$SBS7a+sbs_sig$SBS7b+sbs_sig$SBS7c+sbs_sig$SBS7d
# #UV_light_exposure_indirect_effect(SBS38)
# sbs_sig$UV_light_exposure_indirect_effect <- sbs_sig$SBS38
# #
# sbs_sig$aging <- sbs_sig$SBS1+sbs_sig$SBS5+sbs_sig$SBS40


#APOBEC_activity(SBS2 and SBS13)
# sbs_sig$APOBEC_activity <- sbs_sig$SBS2+sbs_sig$SBS13
# #Polymerase eta somatic hypermutation
# sbs_sig$Polymerase <- sbs_sig$SBS9
# #BER deficiency
# sbs_sig$BER_deficiency <- sbs_sig$SBS36
#drug mutagenesis
sbs_sig$drug_mutagenesis <- sbs_sig$SBS17a+sbs_sig$SBS17b+sbs_sig$SBS31+sbs_sig$SBS32+sbs_sig$SBS35
#Unknown
# sbs_sig$Unknown <- sbs_sig$SBS33+sbs_sig$SBS34+sbs_sig$SBS19+sbs_sig$SBS54 #新的结果没有sbs19
sbs_sig$Unknown <- sbs_sig$SBS33+sbs_sig$SBS34+sbs_sig$SBS54
#
sbs_sig <- sbs_sig[,c("sample","SBS7a","SBS7b","SBS7c","SBS7d","SBS38","SBS2",
                      "SBS13","drug_mutagenesis","SBS9","Unknown","SBS36","SBS1","SBS5","SBS40")]

# sbs_sig <- sbs_sig[,c("sample","UV_light_exposure","UV_light_exposure_indirect_effect","APOBEC_activity","Polymerase","BER_deficiency","drug_mutagenesis","Unknown","aging")]
#
#DBS
dbs <- read.table("Decomposed_Solution_Activities_SBSDINUC.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(dbs)[1] <- "sample"
#ID
id <- read.table("Decomposed_Solution_Activities_SBSINDEL.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(id)[1] <- "sample"
#WGS
wgd <- read.table("01.purity_WGD.txt",sep="\t", quote="", header=T,stringsAsFactors = F)
wgd <- wgd[,c(1,4)]
names(wgd) <- c("sample","WGD")
#clinical
clinical<-read.table("clinical_55.txt", sep="\t", quote="", header=TRUE,na.strings = "")
clinical <- clinical[,c(1,2,5,4,10,7,9)]
names(clinical) <- c("sample","name","Gender","Age","Site","Stage","Type")
sample_twt <- read.table("sample_TWT.txt",header = T,stringsAsFactors = F,sep = "\t")
names(sample_twt) <- c("sample","TWT")
clinical <- merge(clinical,sample_twt,by="sample",all=T)
clinical$TWT[is.na(clinical$TWT)] <- "TWT"
clinical <- merge(clinical,wgd,by="sample")


clinical$Type[which(clinical$Type=="metastasis")] <-20
clinical$Type[which(clinical$Type=="primary")] <- 19
clinical$Type[which(clinical$Type=="NA")] <- 25



#Age

# clinical$Age[which(clinical$Age=="NA")] <- 25
# clinical$Age <- as.numeric(clinical$Age)
age <- median(na.omit(clinical$Age) )
print(age)
clinical$age[which(clinical$Age >=58)] <- ">58"#"70-79y"
clinical$age[which(clinical$Age <58)] <- "<58"
clinical$age[which(clinical$Age =="NA")] <- "NA"
#
clinical$Age[which(clinical$age =="NA")] <- 25#"70-79y"
clinical$Age[which(clinical$age %in% c(">58"))] <- 7
clinical$Age[which(clinical$age %in% c("<58"))] <- 4

# clinical$Age[which(clinical$Age >10 & clinical$Age<40 )] <-2# "15-39y"
# clinical$Age[which(clinical$Age<50  & clinical$Age>39)] <- 3#"40-49y"
# clinical$Age[which(clinical$Age<60 & clinical$Age>49)] <- 4#"50-59y"
# clinical$Age[which(clinical$Age<70  & clinical$Age>59)] <-5#"60-69y"
# clinical$Age[which(clinical$Age<80 & clinical$Age>69)] <- 6#"70-79y"
# clinical$Age[which(clinical$Age<90 & clinical$Age>79)] <- 7#"80-89y"

#Stage
# clinical$Stage[which(clinical$Stage=="Ib")] <- "I"
# clinical$Stage[which(clinical$Stage=="IIIC")] <- "III"
# clinical$Stage[which(clinical$Stage %in% c("IIc","II"))] <- "II"
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
#twt
clinical$TWT[which(clinical$TWT=="BRAF")] <- 34
clinical$TWT[which(clinical$TWT=="NRAS")] <- 35
clinical$TWT[which(clinical$TWT=="NF1")] <- 36
clinical$TWT[which(clinical$TWT=="TWT")] <- 37
############


###merge
data <- merge(df,clinical,by="sample")
data <- data[order(data$perMb,decreasing=F),]
data <- data[order(data$Age,data$Gender,data$Stage,decreasing=F),]
data <- data[order(data$Site,decreasing=F),]
sample_order <- as.data.frame(data$sample)
names(sample_order) <- "sample"
######################################################

# sig <- merge(sample_order,sbs,by="sample",sort=F)
sig <- sbs_sig
sig[2:ncol(sig)] <- sig[2:ncol(sig)]/ rowSums(sig[2:ncol(sig)])

# sig <- sig[order(sig[,1],sig[,2],sig[,3],sig[,4],sig[,5],decreasing = T),]
#####仅供排序

combin <- merge(data,sig,by="sample")
combin$uv <- combin$SBS7a+combin$SBS7b+combin$SBS7c+combin$SBS7d
combin <- combin[order(combin$Age,decreasing=F),]
combin <- combin[order(combin$Stage,decreasing=T),]
# clinical<-merge(sample_order,clinical,by=0,sort=FALSE,all=T)
combin<-combin[order(combin$uv,combin$SBS7a,combin$SBS7b,combin$SBS7c,combin$SBS7d,combin$SBS38,decreasing = T),]
# combin<-combin[order(combin[,"Site"],decreasing = F),]
library(dplyr)
combin <- select(combin,-uv)
sample_order <- as.data.frame(combin$sample)
names(sample_order) <- "sample"
############
sig <- sig[match(sample_order$sample,sig$sample),]
rownames(sig) <- sig$sample
sig <- as.matrix(sig[,c(-1)])
signatre.p <- t(sig)

#dbs_sig
dbs_sig <- merge(sample_order,dbs,by="sample",all = T,)
dbs_sig[is.na(dbs_sig)] <- 0
dbs_sig <- dbs_sig[match(sample_order$sample,dbs_sig$sample),]
rownames(dbs_sig) <- dbs_sig$sample
dbs_sig <- dbs_sig[,c(-1)]
dbs_sig$other <- 0

dbs_sig.t <- as.matrix(t(dbs_sig))
dbs_signatre.p <- t(t(dbs_sig.t)/colSums(dbs_sig.t))
dbs_signatre.p[is.na(dbs_signatre.p)] <- 0
dbs_signatre.p["other",] <-  (1- colSums(dbs_signatre.p))


###
id_sig <- merge(sample_order,id,by="sample",all=T)
id_sig <- id_sig[match(sample_order$sample,id$sample),]
rownames(id_sig) <- id_sig$sample
id_sig <- id_sig[,c(-1)]

id_sig.t <- as.matrix(t(id_sig))
id_signatre.p <- t(t(id_sig.t)/colSums(id_sig.t))
id_signatre.p[is.na(id_signatre.p)] <- 0
# id_signatre.p["other",] <- (1- colSums(id_signatre.p)  )


############
data <- data[match(sample_order$sample,data$sample),]
clin <- data[,c("Stage","Type","Age","Gender","WGD","Site","TWT")]
names(clin) <- c("Stage","Specimen type","Age","Gender","WGD","Primary site","driver")
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
col2[class==22]<-"#FF44AA"
#WGD "#555555","#99CC99"
col2[class==23]<-"#555555"
col2[class==24]<-"#99CC99"
#
col2[class==25]<-"#d3d3d3"
#twt
#"#FF33FF","#990066","#99FF99","#CC9999"
col2[class==34]<-"#FF33FF"
col2[class==35]<-"#990066"
col2[class==36]<-"#99FF99"
col2[class==37]<-"#CC9999"
cols <-  c("#ed1299", "#ddd53e", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92")

##########################################################################
pdf("china_signature_clin.pdf",width=15, height=9 )
layout(matrix(c(1:6), 3, 2, byrow = TRUE),  heights=c(1.5,2,4),widths = c(3,1))
###############
par(mar=c(1, 11.5, 3, 0), xpd=TRUE)
barplot(data$perMb, col="#0081a7", border="#0081a7", space=0, 
        names.arg=NULL, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=1.7, font=2, las=1);
axis(2, at=c(0,4,8,12,16), labels=c(0,4,8,12,16), las=1, cex.axis=1.5,font = 1,line=0);
mtext(side = 2, text = "Mutation per\n megabase", line = 3, cex=1.2,font = 1)

par(mar=c(5, 1, 6, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)

############
#SBS n=9
col1 <-c("#800f2f","#a4133c","#c9184a","#ff4d6d",
         "#66FF66","#CCCCFF","#7744FF","#FF8C00","#fdffb6",
         "#6c757d","#caffbf","#adb5bd","#CC99CC","#993399")
par(mar=c(0, 11.5, 1, 0), xpd=TRUE)
barplot(signatre.p,col =col1,space = 0,border=FALSE, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=2, font=2, las=1)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c("0","20","40","60","80","100"), las=1, cex.axis=1.5,font = 1,line=0);
mtext(side = 2, text = "Proportion of SBS(%)", line = 3, cex=1.2,font = 1)

par(mar=c(5, 0, 0, 8), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=col1[1:nrow(signatre.p)],legend=rownames(signatre.p), title="SBS96", title.adj = 0, border=FALSE, cex=1.7,horiz=F,ncol = 2,bty="n")
##


#dbs_sig n=6+grep
# col3 <-  c("#003C9D","#DDAA00","#00AA00","#FF3333", "#8000FF","#ff8000","grey")
# col3 <-c("#800f2f","#a4133c","#ff4d6d","#ff8fa3","#ddb892","#d5bdaf","#a26769","#6c757d","#adb5bd","#48bfe3","#bde0fe")
# par(mar=c(0, 11.5, 1, 0), xpd=TRUE)
# barplot(dbs_signatre.p,col =col3,space = 0.1,border=FALSE, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=2, cex.axis=1.7, font=2, las=1)
# axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c("0","20","40","60","80","100"), las=1, cex.axis=1.5,font = 1,line=-3);
# mtext(side = 2, text = "Proportion of DBS(%)", line = 2, cex=1.1,font = 1)
# 
# par(mar=c(5, 0, 0, 0), xpd=TRUE)
# plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
# legend(3,11, pch=15,  col=col3,legend=rownames(dbs_signatre.p), title="DBS78", title.adj = 0, border=FALSE, cex=2,horiz=F,ncol = 2,bty="n")


#id_sig
# col4 <-c("#ddb892","#ff4d6d","#bde0fe")
# par(mar=c(0, 11.5, 1, 0), xpd=TRUE)
# barplot(id_signatre.p,col =col4,space = 0.1,border=FALSE, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=2, font=2, las=1)
# axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c("0","20","40","60","80","100"), las=1, cex.axis=1.5,font = 1,line=-3);
# mtext(side = 2, text = "Proportion of IDS(%)", line = 2, cex=1.1,font = 1)
# 
# par(mar=c(5, 0, 0, 0), xpd=TRUE)
# plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
# legend(3,11, pch=15,  col=col4,legend=rownames(id_signatre.p), title="ID83", title.adj = 0, border=FALSE, cex=2,horiz=F,ncol = 2,bty="n")
#
# axis(2, at=c(0,10,20), labels=c(20,30,40),las=1, cex.axis=1.5,font = 1);
# mtext(side = 2, text = "blood_depth", line = 4.5, cex=1.8,font = 1)
#clin
par(mar=c(22, 14.2, 1, 2.8), xpd=TRUE)
color2D.matplot(class, cellcolors=col2, border = NA, main="", xlab="", ylab="", axes=FALSE,) 
axis(2,at=0.5:(nrow(class)-0.5),labels = rev(rownames(class)),las=2,cex.axis=1.7,font=1,line = F,tick = F)
for(i in (0:(ncol(class)))){
  segments(i,0,i,nrow(class),col = "black",lwd=0)
}
for(i in (0:(nrow(class)))){
  segments(0,i,ncol(class),i,col = "black",lwd=1)
}
# color2D.matplot(plot_df, cellcolors=col, border = "white", main="", xlab="", ylab="", axes=FALSE)
axis(1, at=c(1:length(data$name))-0.5,labels=data$name, line=FALSE, tick=FALSE, las=2, cex.axis=1,font=1,mgp=c(0, 1,1))

par(mar=c(0, 0, 0, 1), xpd=TRUE)
plot(x=5,y=11,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,16, pch=15,  col=c("#EDE0D4","#E6CCB2","#DDB892","#9C6644","#d3d3d3"),legend=c("I","II","III","IV","NA"), title="", title.adj = 0, border=FALSE, cex=1.5,horiz=T,bty="n")
legend(3,15.5, pch=15,  col=c("#03827f","#f94144"),legend=c("Primary","Metastasis"),title="", title.adj = 0, border=FALSE,  cex=1.5,horiz=T,bty="n")
legend(3,15, pch=15,  col=c("#ff8fa3","#590d22","#d3d3d3"),legend=c("<58",">=58","NA"), title="",title.adj = 0, border=FALSE, cex=1.5,horiz=T,bty="n")
legend(3,14.5, pch=15,  col=c("#fbb1bd","#bbd0ff"),legend=c("F","M"),title="", title.adj = 0, border=FALSE,  cex=1.5,horiz=T,bty="n")
# legend(3,8, pch=15, col=c("#48D1CC","#996666"),legend=c("upper","lower"), title="",title.adj = 0, border=FALSE, cex=2,horiz=T,bty="n")
legend(3,14, pch=15, col=c("#555555","#99CC99"),legend=c("non-WGD","WGD"), title="",title.adj = 0, border=FALSE, cex=1.5,horiz=F,ncol = 2,bty="n")
legend(3,13.5, pch=15, col=c("#CC9933","#D28EFF","#00A087FF","#6633FF"),legend=c("hand","heel","sole","other(foot)"), title="",title.adj = 0, border=FALSE, cex=1.5,horiz=F,ncol = 2,bty="n")
legend(3,12.5, pch=15, col=c("#FF33FF","#990066","#99FF99","#CC9999"),legend=c("BRAF","NRAS","NF1","Triple WT"), title="",title.adj = 0, border=FALSE, cex=1.5,horiz=F,ncol = 2,bty="n")

#
#Site
dev.off()

