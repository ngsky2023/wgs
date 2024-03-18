args=commandArgs(trailingOnly=T)
maf=args[1]
gene=args[2]

library(maftools)
library(stringr)
laml = read.maf(maf = maf)
labpos=sort(as.numeric(unique(str_extract(unique(laml@data$HGVSp_Short), "[[:digit:]]+"))))
pdf(paste(gene,"_mutation.pdf",sep=""),height=2.5,width=8)
lollipopPlot(maf = laml, gene = gene , AACol = 'HGVSp_Short', labelPos=labpos ,labPosAngle = 0,showMutationRate = FALSE,axisTextSize = c(1,1),defaultYaxis=TRUE,labelOnlyUniqueDoamins=TRUE)
dev.off()
