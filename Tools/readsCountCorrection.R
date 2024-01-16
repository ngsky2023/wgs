library('TitanCNA')
#library(doMC)
#registerDoMC()
require(HMMcopy)
args = commandArgs(trailingOnly=TRUE)

normWig <- args[1]
gcWig <- args[2]
mapWig <- args[3]
out <- args[4]

setGenomeStyle <- function(x, genomeStyle = "NCBI", species = "Homo_sapiens"){
	#chrs <- genomeStyles(species)[c("NCBI","UCSC")]
	if (!genomeStyle %in% seqlevelsStyle(as.character(x))){
    	x <- suppressWarnings(mapSeqlevels(as.character(x), 
    					genomeStyle, drop = FALSE)[1,])
    }
    
    autoSexMChr <- extractSeqlevelsByGroup(species = species, 
    				style = genomeStyle, group = "all")
    x <- x[x %in% autoSexMChr]
    return(x)
}

genomeStyle = "NCBI"
message("Reading GC and mappability files")
gc <- wigToRangedData(gcWig)
map <- wigToRangedData(mapWig)

### LOAD TUMOUR AND NORMAL FILES ###
message("Loading normal file:", normWig)
normal_reads <- wigToRangedData(normWig)

### set the genomeStyle: NCBI or UCSC
#require(GenomeInfoDb)
names(gc) <- setGenomeStyle(names(gc), genomeStyle)
names(map) <- setGenomeStyle(names(map), genomeStyle)
names(normal_reads) <- setGenomeStyle(names(normal_reads), genomeStyle)

### make sure tumour wig and gc/map wigs have same
### chromosomes
gc <- gc[gc$space %in% normal_reads$space, ]
map <- map[map$space %in% normal_reads$space, ]
samplesize <- 50000

### add GC and Map data to IRanges objects ###
normal_reads$gc <- gc$value
normal_reads$map <- map$value
colnames(normal_reads) <- c("reads", "gc", "map")

### CORRECT NORMAL DATA FOR GC CONTENT AND
### MAPPABILITY BIASES ###
message("Correcting Normal")
normal_copy <- correctReadcount(normal_reads, samplesize = samplesize)
 



write.table(normal_copy, file=out, sep="\t", col.names=T, row.names=T, quote=F)

