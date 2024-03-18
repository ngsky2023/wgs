###add top genes on by frequency
rm(list=ls())
library(reshape2)
library(RColorBrewer)
library(tidyverse)
library(plotrix)
library(scales)
library(ggsci)

Frequency<-function(df=df){
  freq <-NULL
  a <- NULL
  for(i in as.character(unique(df$sample))){
    for (j in as.character(unique(df$type))) {
      tmp<-df[which(df$sample==i & df$type==j),]
      num<-sum(tmp$number)
      freq <- c(i,j,num)
      a <- as.data.frame(rbind(a,freq))
    }
    
  }
  names(a) <-  c("sample","type","number")
  return(a)
}

#read data
df<-read.table("mutation_burden.xls", sep="\t", quote="", header=T,stringsAsFactors = F)
df<- df[,c(1,5,6,7)]
df <- df[order(df$perMb,decreasing = F),]
#base_num
base_num <- read.table("SBS96.Samples.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
base_num <- melt(base_num)
base_num<-separate(data = base_num, col = Mutation.Types, into = c("A5", "Mutation.Types"), sep = "\\[")
base_num<-separate(data = base_num, col = Mutation.Types, into = c("MutationTypes", "type2"), sep = "\\]")
base_num <- base_num[,c(4,2,5)]
names(base_num) <- c("sample","type","number")
b <- Frequency(df=base_num)

# names(df) <- c("sample","perlMb")
sbs <- read.table("Decomposed_Solution_Activities_SBS96.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(sbs)[1] <- "sample"
# sbs <- sbs[,c(1,4,2,3,5,6)]
# names(sbs) <- c("sample","SBS7a","SBS5","SBS3-like","SBS38","SBS31")
#dbs_num
dbs_num <- read.table("DBS78.Samples.txt",sep="\t",quote = "",header = F,stringsAsFactors = F)
dbs_num <- as.data.frame(t(dbs_num))
colnames(dbs_num) <- dbs_num[1,]
dbs_num <- dbs_num[-1,]
dbs_num <- dbs_num[,c("Mutation Types","CC>AA","CC>TT","CC>Other","CT","GC","TC","TG","AC","TT","Other")]
names(dbs_num) <- c("sample","CC>AA","CC>TT","CC>Other","CT","GC","TC","TG","AC","TT","Other")



#DBS
dbs <- read.table("Decomposed_Solution_Activities_SBSDINUC.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(dbs)[1] <- "sample"

#ID_num
id_num <- read.table("ID83.Samples.txt",sep="\t",quote = "",header = F,stringsAsFactors = F)
id_num <- as.data.frame(t(id_num))
colnames(id_num) <- id_num[1,]
id_num <- id_num[-1,]
id_num <- id_num[,c("Mutation Types","Del-C","Del-T","Del-MH","Ins-C","Ins-T","Del-repeats","Ins-repeats")]
names(id_num) <- c("sample","Del-C","Del-T","Del-MH","Ins-C","Ins-T","Del-repeats","Ins-repeats")


#ID
id <- read.table("Decomposed_Solution_Activities_SBSINDEL.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(id)[1] <- "sample"

#clinical
#WGS
wgd <- read.table("01.purity_WGD.txt",sep="\t", quote="", header=T,stringsAsFactors = F)
wgd <- wgd[,c(1,4)]
names(wgd) <- c("sample","WGD")

clinical<-read.table("clinical_55.txt", sep="\t", quote="", header=TRUE,na.strings = "")
clinical <- clinical[,c(1,2,5,4,10,7,9)]
names(clinical) <- c("sample","name","Gender","Age","Site","Stage","Type")
clinical <- merge(clinical,wgd,by="sample")


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
clinical$Stage[which(clinical$Stage=="¢ô")] <- 11
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


###merge
data <- merge(df,clinical,by="sample")

data <- data[order(data$Age,data$Gender,data$Stage,decreasing=F),]
data <- data[order(data$perMb,decreasing=T,data$Site),]
# data <- data[order(data$Region,decreasing=F),]
sample_order <- as.data.frame(data[,1])
names(sample_order) <- "sample"

# data <- merge(sample_order,data,by="sample")

# sig <- merge(sample_order,sbs,by="sample",sort=F)
# sig <- sbs
# rownames(sig) <- sig$sample
# sig <- as.matrix(sig[,c(-1)])
# sig <- sig/rowSums(sig)
# sig <- sig[order(sig[,1],sig[,2],sig[,3],sig[,4],sig[,5],decreasing = T),]
# signatre.p <- as.matrix(t(sig))
# sample_order <- as.data.frame(rownames(sig))
# names(sample_order) <- "sample"

# data <- data[match(sample_order$sample,data$sample),]

# freq <- NrowSums()# freq <- NULL
# for(i in 1:nrow(sig.t)){
#   num <- 0;
#   for(j in 1:ncol(sig.t)){
#     num <- num+sig.t[i,j]
#   }
#   freq <- c(freq,num)
# }
# sig.t <- as.data.frame(sig.t)
# sig.t$Freq <- freq
# sig.t <- sig.t[order(sig.t[,"Freq"],decreasing = T),]
# sig.t <- subset(sig.t,select=-Freq)
# signatre.p <- t(sig)
###########
# base_num <- base_num[match(sample_order$sample,base_num$sample),]
# rownames(base_num) <- base_num$sample
# base_num <- base_num[,c(-1)]
# base_num.t <- as.matrix(t(base_num))
# base_num.p <- t(t(base_num.t)/colSums(base_num.t))
######
#dbs_num
dbs_num <- merge(sample_order,dbs_num,by="sample",all=T)
dbs_num[is.na(dbs_num)] <- 0
dbs_num <- dbs_num[match(sample_order$sample,dbs_num$sample),]
rownames(dbs_num) <- dbs_num$sample
dbs_num <- dbs_num[,-1]
dbs_num <- as.matrix(dbs_num)
dbs_num.t <- t(dbs_num)
#dbs_sig
dbs_sig <- merge(sample_order,dbs,by="sample",all = T,)
dbs_sig[is.na(dbs_sig)] <- 0
dbs_sig <- dbs_sig[match(sample_order$sample,dbs_sig$sample),]
rownames(dbs_sig) <- dbs_sig$sample
dbs_sig <- dbs_sig[,c(-1)]
dbs_sig$empty <- 0

dbs_sig.t <- as.matrix(t(dbs_sig))
dbs_signatre.p <- t(t(dbs_sig.t)/colSums(dbs_sig.t))
dbs_signatre.p[is.na(dbs_signatre.p)] <- 0
dbs_signatre.p["empty",] <-  (1- colSums(dbs_signatre.p))

####
id_num <- merge(sample_order,id_num,by="sample",all=T)
id_num[is.na(id_num)] <- 0
id_num <- id_num[match(sample_order$sample,id_num$sample),]
rownames(id_num) <- id_num$sample
id_num <- id_num[,-1]
id_num <- as.matrix(id_num)
id_num.t <- t(id_num)

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
clin <- data[,c("Stage","Type","Age","Gender","WGD","Site")]
names(clin) <- c("Stage","Specimen type","Age","Gender","WGD","Primary site")
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
##########################################################################
pdf("sig_dbs_id_barplot.pdf",width=28, height=19)
layout(matrix(c(1:14), 7, 2, byrow = TRUE),  heights=c(1.5,1.5,1.5,1.5,1.5,3.5),widths = c(3,1))
#TMB
par(mar=c(0, 11.5, 8, 3.1), xpd=TRUE)
barplot(data$perMb, col="#A9A9A9", border="#A9A9A9", space=0.1, 
        names.arg=NULL, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=2, cex.axis=3, font=2, las=1);
axis(2, at=c(0,4,8,12,16), labels=c(0,4,8,12,16), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "Mutation per \nmegabase", line = 4.5, cex=1.7,font = 1)

par(mar=c(5, 1, 10, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)



##dbs_num
# col5 <-  c("#003C9D","#DDAA00","#00AA00","#FF3333", "#8000FF","#ff8000","#FF44AA","#660000","#003333","#99FFFF")
col5 <-  c(pal_npg(alpha=0.5)(8),pal_npg()(8)[8],"grey")
par(mar=c(1, 11.5, 5, 3.1), xpd=TRUE)
barplot(dbs_num.t,col =col5,space = 0.1, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=1.7, font=2, las=1)
axis(2, at=c(0,100,200,300,400), labels=c(0,100,200,300,400), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "Number of DNVs", line = 5.5, cex=1.7,font = 1)

par(mar=c(5, 0, 0, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=col5,legend=rownames(dbs_num.t), title="", title.adj = 0, border=FALSE, cex=2.5,horiz=F,ncol = 3,bty="n")


#dbs_sig
# col3 <-  c("#003C9D","#DDAA00","#00AA00","#FF3333", "#8000FF","#ff8000","#FF44AA","#660000","#003333")
col3 <- c(pal_npg(alpha=0.8)(8),"grey")
par(mar=c(0, 11.5, 1, 3.1), xpd=TRUE)
barplot(dbs_signatre.p,col =col3,border="grey",space = 0.1, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=1.7, font=2, las=1)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c("0","20","40","60","80","100"), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "Proportion of DNV\nchanges(%)", line = 5.5, cex=1.7,font = 1)

par(mar=c(5, 0, 0, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=col3,legend=rownames(dbs_signatre.p), title="DBS78", title.adj = 0, border=FALSE, cex=2.5,horiz=F,ncol = 3,bty="n")

###id_num
# col5 <-  c("#003C9D","#DDAA00","#00AA00","#FF3333", "#8000FF","#ff8000","#FF44AA", rainbow(8)[5])
col5 <- pal_npg(alpha=0.5)(7)
par(mar=c(1, 11.5, 1, 3.1), xpd=TRUE)
barplot(id_num.t,col =col5,space = 0.1, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=1.7, font=2, las=1)
axis(2, at=c(0,500,1000,1500,2000), labels=c(0,500,1000,1500,2000), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "Number of In dels", line = 5.5, cex=1.7,font = 1)

par(mar=c(5, 0, 0, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=col5,legend=rownames(id_num.t), title="", title.adj = 0, border=FALSE, cex=2.5,horiz=F,ncol = 2,bty="n")

#id_sig
# col4 <-  c("#003C9D","#DDAA00","#00AA00","#FF3333", "#8000FF","#ff8000","#FF44AA","grey")
col4 <-pal_npg(alpha=0.8)(5)
par(mar=c(0, 11.5, 1, 3.1), xpd=TRUE)
barplot(id_signatre.p,col =col4,space = 0.1, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=3, cex.axis=3, font=2, las=1)
axis(2, at=c(0,0.2,0.4,0.6,0.8,1), labels=c("0","20","40","60","80","100"), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "Proportion of ID \nchanges(%)", line = 4.5, cex=1.7,font = 1)

par(mar=c(5, 0, 0, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=col4,legend=rownames(id_signatre.p), title="ID83", title.adj = 0, border=FALSE, cex=3,horiz=F,ncol = 3,bty="n")
#
# axis(2, at=c(0,10,20), labels=c(20,30,40),las=1, cex.axis=1.5,font = 1);
# mtext(side = 2, text = "blood_depth", line = 4.5, cex=1.8,font = 1)
#clin
par(mar=c(22.5, 16.8, 1, 8.2), xpd=TRUE)
color2D.matplot(class, cellcolors=col2, border = NA, main="", xlab="", ylab="", axes=FALSE,) 
axis(2,at=0.5:(nrow(class)-0.5),labels = rev(rownames(class)),las=2,cex.axis=2.5,font=1,line = F,tick = F)
for(i in (0:(ncol(class)))){
  segments(i,0,i,nrow(class),col = "black",lwd=0)
}
for(i in (0:(nrow(class)))){
  segments(0,i,ncol(class),i,col = "black",lwd=1)
}
# color2D.matplot(plot_df, cellcolors=col, border = "white", main="", xlab="", ylab="", axes=FALSE)
axis(1, at=c(1:length(data$name))-0.5,labels=data$name, line=FALSE, tick=FALSE, las=2, cex.axis=2,font=1,mgp=c(1, 3,1))

par(mar=c(0, 0, 8, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,13.5,pch=15,  col=c("#EDE0D4","#E6CCB2","#DDB892","#9C6644","#d3d3d3"),legend=c("I","II","III","IV","NA"), title="", title.adj = 0, border=FALSE, cex=3,horiz=T,bty="n")
legend(3,13, pch=15,  col=c("#03827f","#f94144"),legend=c("Primary","Metastasis"),title="", title.adj = 0, border=FALSE,  cex=3,horiz=T,bty="n")
legend(3,12.5, pch=15,  col=c("#ff8fa3","#590d22","#d3d3d3"),legend=c("<58",">=58","NA"), title="",title.adj = 0, border=FALSE, cex=3,horiz=T,bty="n")
legend(3,12, pch=15,  col=c("#fbb1bd","#bbd0ff"),legend=c("F","M"),title="", title.adj = 0, border=FALSE,  cex=3,horiz=T,bty="n")
legend(3,11.5, pch=15, col=c("#555555","#99CC99"),legend=c("non-WGD","WGD"), title="",title.adj = 0, border=FALSE, cex=3,horiz=F,ncol = 2,bty="n")
legend(3,11, pch=15, col=c("#CC9933","#D28EFF","#00A087FF","#6633FF"),legend=c("hand","heel","sole","other(foot)"), title="",title.adj = 0, border=FALSE, cex=3,horiz=F,ncol = 2,bty="n")
dev.off()
#################################

