#!/usr/bin/Rscript
args <-(commandArgs(TRUE))
if(args[1]=='-h' | args[1]=='--help'){print("Usage: ~/LRR_Calculation.R <Input text file of ExomeSomatic> <Read depth (.sample_interval_summary) from GATK DepthOfCoverage> <Output LRR result>")} else {
Input <- read.table(args[1], header=FALSE, sep="\t", col.names = paste("V",seq_len(10),sep=""), fill = TRUE)
Depth <- read.table(args[2], header=TRUE, row.names=1, sep="\t")
LRR_ALL<-0;headers<-0;j=1;
for (i in seq(1,nrow(Input)-1,by=2)){
  NID=paste(as.character(Input[i,2]),"_mean_cvg",sep="")
  TID=paste(as.character(Input[i+1,2]),"_mean_cvg",sep="")
  Nmean=mean(Depth[,NID])
  Tmean=mean(Depth[,TID])
  LRR=as.matrix(log((Depth[,TID]/Tmean)/(Depth[,NID]/Nmean),2),byrow=FALSE)
  LRR_ALL=cbind(LRR_ALL,LRR)
  rm(NID,TID,Nmean,Tmean,LRR)
  headers[j]=as.character(Input[i+1,2])
  j=j+1
}
LRR_ALL<-cbind(rownames(Depth),LRR_ALL[,-1])
write.table(LRR_ALL, file=args[3], 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=c("Target",headers))
}
