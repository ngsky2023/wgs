#!/share/public/software/R/3.5.3/lib64/R/bin/Rscript
library(maftools)

laml.gistic = readGistic(gisticAllLesionsFile = 'all_lesions.conf_95.txt',
                         gisticAmpGenesFile = 'amp_genes.conf_95.txt', 
                         gisticDelGenesFile = 'del_genes.conf_95.txt', 
                         gisticScoresFile = 'scores.gistic', 
                         isTCGA = F)
pdf("gene_bar.pdf",height = 5,width = 10)
gisticChromPlot(gistic = laml.gistic, ref.build = 'hg19',cytobandTxtSize=0.4,
                cytobandOffset=0.12,txtSize=0.6,markBands = "all",fdrCutOff = 0.05)
gisticBubblePlot(gistic = laml.gistic,txtSize = 0.3)
dev.off()
