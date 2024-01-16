library('TitanCNA')
library(doMC)
registerDoMC()
options(cores=4)
args = commandArgs(trailingOnly=TRUE)


id <- as.character(args[1])
chr <- args[2]
bam <- as.character(args[3])
bai <- paste(bam, 'bai', sep=".")
vcf <- paste(id, chr, 'normHetSNPs.vcf', sep=".")
infile <- paste(id, chr, "AlleleReadCounts.txt", sep=".")

extractAlleleReadCounts(bam, bai, vcf, outputFilename = infile)
