library(plotrix)
library(RColorBrewer)
library(argparser)

argv <- arg_parser("")
argv <- add_argument(argv,"--file1", help="F1_MutGene_Matrix.txt")
argv <- add_argument(argv,"--file2", help="F2_MutGene_Stats.txt")
argv <- add_argument(argv,"--file3", help="F3_NumMut_Sorted.txt")
argv <- add_argument(argv,"--out", help="out")
argv <- parse_args(argv)
F1_MutGene_Matrix <- argv$file1
F2_MutGene_Stats <- argv$file2
F3_NumMut_Sorted <- argv$file3
out <- argv$out

F1<-read.table(F1_MutGene_Matrix, sep="\t", quote="", header=TRUE, row.names=1)
F2<-read.table(F2_MutGene_Stats, sep="\t", quote="", header=TRUE, row.names=1)
F3<-read.table(F3_NumMut_Sorted, sep="\t", quote="", header=TRUE, row.names=1)

mutgene<-as.matrix(F1)
freq=0-F2$MutFreq

col<-matrix(rep("#000000",length(mutgene)),nrow=nrow(mutgene))
col[mutgene==0]<-"lightgrey"
col[mutgene==1]<-"#2B406D"
col[mutgene==2]<-"#5F7133"
col[mutgene==3]<-brewer.pal(9,"Set1")[8]
col[mutgene==4]<-brewer.pal(9,"Set1")[1]
col[mutgene==5]<-"#7030A0"
col[mutgene==6]<-"black"

F3$NumMut_MB[F3$NumMut_MB>20] = 20

loc1<-barplot(F3$NumMut_MB, space=0.4, plot=FALSE)
loc2<-barplot(rev(freq), horiz=TRUE, plot=FALSE)

pdf(out,width=10, height=7)
layout(matrix(c(1,1,1,2,3,4,5,5,5),3, 3, byrow = TRUE), widths=c(8,1,1), heights=c(1,2.5,0.3))

par(mar=c(0, 4.4, 1, 13), xpd=TRUE)
barplot(F3$NumMut_MB, col="#658B8E", border="#658B8E", space=0.2,	#col and border set color
names.arg="", ylab="TMB", xlab="", main="", cex.lab=1, cex.axis=1, font=2, las=1, ylim=c(0,20));

par(mar=c(1, 6.5, 1, 0) , xpd=TRUE)
color2D.matplot(mutgene, cellcolors=col, border = NA, lwd=20, main="", xlab="", ylab="", axes=FALSE)
axis(2,at=0.5:(nrow(F1)-0.5), labels = rev(row.names(F2)),las=1, cex.axis=1.4, font=4, line=FALSE,tick=FALSE)
for(i in (2:ncol(F1)-1)){
	segments(i,0,i,nrow(F1),col = "white",lwd=1)
}
for(i in (2:nrow(F1)-1)){
	segments(0,i,ncol(F1),i,col = "white",lwd=1)
}

par(mar= c(0, 0.2, 0, 0), xpd=TRUE)
barplot(rev(freq), horiz=TRUE, xaxt="n", xlim=c(-0.8, 0), xlab="", col="#5A7BA4", border=NA,space=0.25) 
axis(1, at=c(-0.8,-0.6, -0.4, -0.2, 0), labels=c("0.8","0.6", "0.4", "0.2", "0"), las=1, cex.axis=1,font = 2); 
mtext(side = 1, text = "Frequency", line = 3, cex=1,font = 1)

par(mar=c(0, 0, 0, 2), xpd=TRUE)
barplot(rev(F2$MinusLogP), horiz=TRUE, xlab="", col="#B36348", border=NA, cex.axis=1,font=1, space=0.25); 
mtext(side = 1, text = "-LogP", line = 3, cex=1,font = 1)

par(mar=c(0,12,0,24),xpd=TRUE)
plot(x=5,y=8,cex=0,xlab="",ylab="",xaxt="n",yaxt="n",frame.plot=FALSE)
legend("top",pch=15,col=c("lightgrey","#2B406D","#5F7133",border=brewer.pal(9,"Set1")[8],border=brewer.pal(9,"Set1")[1],"#7030A0","black"), border=FALSE,legend=c("No mutation", "Missense", "Splice-site", "Stop-lost", "Stop-gain", "In-frame indel", "Frame shift indel" ),cex=1,ncol=7,bty='n')
dev.off()
#################################################################################################################################################################
