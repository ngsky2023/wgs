#!/usr/bin/Rscript
args <-(commandArgs(TRUE))
if(args[1]=='-h' | args[1]=='--help'){print("Usage: ~/Ranking.R <Formated tumor depth data> <Formated normal depth data> <Output file name>")} else {
Tumor <- read.table(args[1], header=TRUE, row.names=1, sep="\t")
Normal <- read.table(args[2], header=TRUE, row.names=1, sep="\t")
if(ncol(Tumor)==2){
	TumorNorm<-scale(Tumor[,2],center=FALSE,scale=mean(Tumor[,2]))
}else if(ncol(Tumor)>2){
	TumorNorm<-scale(Tumor[,2:ncol(Tumor)],center=FALSE,scale=colMeans(Tumor[,2:ncol(Tumor)]))
}
NormalNorm<-scale(Normal[,2:ncol(Normal)],center=FALSE,scale=colMeans(Normal[,2:ncol(Normal)]))

sink(args[3])
cat("Target\tGene",sep="")
for (col in 2:ncol(Tumor)){
	cat("\t",colnames(as.matrix(Tumor))[col],"_Normalized\t",colnames(as.matrix(Tumor))[col],"_RankNormal",sep="")
}
cat("\n",sep="")

for (row in 1:nrow(TumorNorm)){
	cat(rownames(as.matrix(Tumor))[row],"\t",as.matrix(Tumor)[row,1],sep="")
	for (col in 1:ncol(TumorNorm)){
		quantile=0
		quantile=as.vector(rank(c(TumorNorm[row,col],NormalNorm[row,])))[1]/length(c(TumorNorm[row,col],NormalNorm[row,]))
		cat("\t",TumorNorm[row,col],"\t",quantile,sep="")
	}
	cat("\n",sep="")
}
sink()
}
