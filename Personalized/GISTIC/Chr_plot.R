#!usr/Rscript
library(plotrix)

F<-read.table("result/wgs_Chart.txt", sep="\t", quote="", header=TRUE)
chr<-read.table("/share/work1/hanwj4457/softwore/GISTIC2/gisticplot/ChrInfo_ZhonghuaAll.txt", sep="\t", quote="", header=TRUE)
col1<-ifelse(F$LogQ_Amp>=1.30102999566398, "#E41A1C", "#FBB4AE")
col2<-ifelse(F$LogQ_Del>=1.30102999566398, "#377EB8", "#B3CDE3")
col3<-paste("#",F$ColFreq_Amp,sep="")
col4<-paste("#",F$ColFreq_Del,sep="")
cl<-chr$MaxIndex
ct<-chr$MedianIndex
c<-chr$Chr
ylimG=2
ylimF=1
h<-ifelse(c%%2>0, ylimG-0.03, ylimG-0.07)  ## 设定将“chr N”标示位置高度错开，以免影响相邻标示

pdf("result/Chr_plot.pdf", width=12, height=6)
par(mar=c(1, 5, 1, 2.1))
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), heights=c(2,1,2))

plot(F$Index, F$G_Amp, type="h", col=col1, xaxt="n", frame.plot=FALSE, xlim=c(1,chr$MaxIndex[22]), ylim=c(0,ylimG), 
     xlab="", ylab="G Score Amp", cex.lab=1.5, cex.axis=1, las=1);
text(ct, h, c, cex=1, font=4); abline(v=c(1,cl), col="#999999", lwd=1, lty=3, xpd=T)

plot(F$Index, F$Freq_Amp, type="h", col=col3, xaxt="n", frame.plot=FALSE, xlim=c(1,chr$MaxIndex[22]), ylim=c(0-ylimF,ylimF), 
     xlab="", ylab="Freq of CNA", cex.lab=1.5, cex.axis=1, las=1);
points(F$Index, F$Freq_Del, type="h", col=col4, xaxt="n");abline(v=c(1,cl), col="#999999", lwd=1, lty=3)

plot(F$Index, F$G_Del, type="h", col=col2, xaxt="n", frame.plot=FALSE, xlim=c(1,chr$MaxIndex[22]), ylim=c(0-ylimG,0), 
     xlab="", ylab="G Score Del", cex.lab=1.5, cex.axis=1, las=1);
abline(v=c(1,cl), col="#999999", lwd=1, lty=3)
dev.off()
