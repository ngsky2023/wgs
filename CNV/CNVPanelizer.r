#!/share/public/software/R-3.2.5/bin/Rscript
args <-(commandArgs(TRUE))
if(args[1]=='-h' | args[1]=='--help'){print("Usage: ~/CNVPanelizer.r <Output directory> <BED file> <BAMList tumor> <BAMList Normal>")} else {
library(CNVPanelizer)
out_wd=as.character(args[1])
setwd(out_wd)
### loading annotation files ###
bedFilepath <- as.character(args[2])
amplColumnNumber <- 4
genomicRangesFromBed <- BedToGenomicRanges(bedFilepath, ampliconColumn = amplColumnNumber, split = "_")
metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
geneNames = metadataFromGenomicRanges["geneNames"][, 1]
ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]

### loading data files ###
samp<-read.table(args[3], header=FALSE, sep="", quote="")
sampleFilenames <- as.character(samp$V1)
ref<-read.table(args[4], header=FALSE, sep="", quote="")
referenceFilenames <- as.character(ref$V1)
sampleNames=0;referenceNames=0;
for(i in 1:length(sampleFilenames)){
  sampleNames[i]=strsplit(sampleFilenames,split="/")[[i]][length(strsplit(sampleFilenames,split="/")[[i]])-2]
}
for(i in 1:length(referenceFilenames)){
  referenceNames[i]=strsplit(referenceFilenames,split="/")[[i]][length(strsplit(referenceFilenames,split="/")[[i]])-2]
}

### counting reads ###
removePcrDuplicates <- FALSE
sampleReadCounts <- ReadCountsFromBam(sampleFilenames, genomicRangesFromBed, sampleNames = sampleNames, ampliconNames = ampliconNames, removeDup = removePcrDuplicates)
referenceReadCounts <- ReadCountsFromBam(referenceFilenames, genomicRangesFromBed, sampleNames = referenceNames, ampliconNames = ampliconNames, removeDup = removePcrDuplicates)

### normalization ###
normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts, referenceReadCounts, ampliconNames = ampliconNames)
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]

### bootstrap based CNV ###
replicates <- 10000
bootList <- BootList(geneNames, samplesNormalizedReadCounts, referenceNormalizedReadCounts, replicates = replicates)
backgroundNoise <- Background(geneNames, samplesNormalizedReadCounts, referenceNormalizedReadCounts, bootList, replicates = replicates, significanceLevel = 0.05, robust = TRUE)

### report ###
reportTables <- ReportTables(geneNames, samplesNormalizedReadCounts, referenceNormalizedReadCounts, bootList, backgroundNoise)
tab<-cbind(GeneID=rownames(data.frame(reportTables)), data.frame(reportTables, row.names=NULL))
write.table(tab, file = "LC_CRC_CNVPanelizer_Output.txt", append = FALSE, quote = FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
PlotBootstrapDistributions(bootList, reportTables, outputFolder = out_wd, sampleNames = NULL, save = TRUE, scale = 7)
}
