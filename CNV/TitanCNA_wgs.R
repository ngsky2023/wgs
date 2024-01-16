library('TitanCNA')
library(doMC)
registerDoMC()
options(cores=4)
args = commandArgs(trailingOnly=TRUE)

numClusters <- 2

id <- as.character(args[1])
infile <- paste(id, "AlleleReadCounts.txt", sep=".")
tumWig <- paste(id, "tum_reads.wig", sep=".")
normWig <- paste(id, "norm_reads.wig", sep=".")
gc <- paste(id, "gc.wig", sep=".")
map <- paste(id, "map.wig", sep=".")

data <- loadAlleleCounts(infile)
params <- loadDefaultParameters(copyNumber=5, numberClonalClusters=numClusters, symmetric=TRUE, data=data)

params$normalParams$n_0 <- 0.5
params$ploidyParams$phi_0 <- 2

#cnData <- correctReadDepth(tumWig, normWig, gc, map, genomeStyle = "UCSC")
cnData <- correctReadDepth(tumWig, normWig, gc, map)

logR <- getPositionOverlap(data$chr,data$posn,cnData)
data$logR <- log(2^logR)
rm(logR,cnData)

data <- filterData(data,c(1:22,"X"),minDepth=10,maxDepth=2000)

mScore <- as.data.frame(wigToRangedData(map))
mScore <- getPositionOverlap(data$chr,data$posn,mScore[,-4])
data <- filterData(data,c(1:22,"X"),minDepth=10,maxDepth=2000,map=mScore,mapThres=0.8)

convergeParams <- runEMclonalCN(data,gParams=params$genotypeParams,nParams=params$normalParams,
                                pParams=params$ploidyParams,sParams=params$cellPrevParams,
                                maxiter=20,maxiterUpdate=1500,txnExpLen=1e15,txnZstrength=1e5,
                                useOutlierState=FALSE,
                                normalEstimateMethod="map",estimateS=TRUE,estimatePloidy=TRUE)

optimalPath <- viterbiClonalCN(data,convergeParams)

outfile <- paste(id, "titan.txt",sep=".")
results <- outputTitanResults(data,convergeParams,optimalPath,
                                     filename=outfile,posteriorProbs=FALSE,subcloneProfiles=TRUE)
outparam <- paste(id, "params.txt",sep=".")
outputModelParameters(convergeParams,results,outparam)


#norm <- convergeParams$n[length(convergeParams$n)]
#ploidy <- convergeParams$phi[length(convergeParams$phi)]
#png(paste(id, chr, 'TitanCNA.png', sep="."),width=1200,height=1000,res=100,type="cairo")
#par(mfrow=c(3,1))
#plotCNlogRByChr(results, chr, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-4,6),cex=0.5,main=paste('Chr',chr,sep=" "))
#plotAllelicRatio(results, chr, geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5,main=paste('Chr',chr,sep=" "))
#plotClonalFrequency(results, chr, normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5,main=paste('Chr',chr,sep=" "))
#if (as.numeric(numClusters) <= 2){ 
#	plotSubcloneProfiles(results, chr, cex = 2, spacing=6, main=paste('Chr',chr,sep=" "))
#}
#dev.off()



