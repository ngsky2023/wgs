#!/share/public/software/R/3.5.3/lib64/R/bin/Rscript
library(maftools)

laml.gistic = readGistic(gisticAllLesionsFile = 'all_lesions.conf_95.txt',
                         gisticAmpGenesFile = 'amp_genes.conf_95.txt', 
                         gisticDelGenesFile = 'del_genes.conf_95.txt', 
                         gisticScoresFile = 'scores.gistic', 
                         isTCGA = F)
pdf("gene_bar.pdf",height = 3,width = 7)
#gisticChromPlot(gistic = laml.gistic, ref.build = 'hg19',cytobandTxtSize=0.8,cytobandOffset=0.12,txtSize=0.6)
gisticChromPlot(gistic = laml.gistic, ref.build = 'hg19', markBands = "all")
gisticBubblePlot(gistic = laml.gistic)
dev.off()
