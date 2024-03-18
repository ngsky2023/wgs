rm(list=ls())
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(grid)
library(cowplot)

clin <- read.delim2("clinical_55.txt",header = T,stringsAsFactors = F)
clin <- clin[,c(1,2,7)]
names(clin) <- c("sample","name","Site")
link <- read.table("tra.xls",header=T,stringsAsFactors = F)
names(link)[1] <-"sample" 
link <- link[which(link$TYPE=="BND"),]
link <- merge(clin,link,by="sample")
# hot_chr <- c("chr8","chr16","chr11","chr12","chr19","chr22","chr6","chr5","chr7","chr5","chr17","chr18")
# link <- link[c(which(link$CHROM_A %in% hot_chr & link$CHROM_B %in% hot_chr)),]
sample <- unique(link$name)
####
# for(i in 1:length(sample)){
  # p_list <- lapply(1:length(sample),function(x){
  #   sam <- sample[x]
  #   link_data <- link[grep(sam,link$name),c(2,4,5,6,7,8,9)]
  #   colnames(link_data) <- c("ID","chr1","start1","end1","chr2","start2","end2")
  #   # link_data$ID <- paste0("ID",1:nrow(link_data))
  #   bed1 <- link_data[,c(2,3,4)];colnames(bed1) <- c('chr','start','end');bed1$value <- rep(0,nrow(bed1))
  #   bed2 <- link_data[,c(5,6,7)];colnames(bed2) <- c('chr','start','end');bed2$value <- rep(0,nrow(bed2))
  #   #基因组的展示
  #   pdf(paste0(sam,'_circos.pdf'),width=5,height=5)
  #   circos.initializeWithIdeogram(axis.labels.cex = 0.2, labels.cex = 0.5)
  #   #
  #   circos.genomicLink(bed1, bed2,col = rand_color(nrow(bed1),transparency = 0.01))
  #   
  #   title = Legend(at = sam, type = "points", ncol=1,
  #                  legend_gp = gpar(col = 'transparent'),background = "transparent",
  #                  labels_gp = gpar(fontsize = 6),size = unit(1, "mm"))
  #   draw(title, x = unit(32, "mm"), y = unit(65, "mm"), just = c("top"))
  #   dev.off()
  # })
# }
############################################################################################################
library(scales)
 col <- as.data.frame(clin$sample)
  col$cols <- rand_color(length(unique(clin$name)),transparency = 0.01)
  show_col(col$cols)
  names(col) <- c("sample","cols")
  link <- merge(link,col)
  # p_list <- lapply(1:length(sample),function(x){
  #   sam <- sample[x]
  #   link_data <- link[grep(sam,link$name),c(2,4,5,6,7,8,9)]
    # hot_chr <- c("chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr2","chr20","chr21","chr22","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chrX")
    # hot_chr_paire <- list(c("chr5","chr1"),c("chr15","chr22"),c("chr5","chr11"),c("chr5","chr7"),c("chr8","chr9"),c("chr16","chr17"),c("chr4","chr11"),c("chr6","chr11"),
    #                       c("chr10","chr16"),c("chr3","chr11"),c("chr5","chr12"),c("chr9","chr18"),c("chr12","chr15"),c("chr2","chr4"),
    #                       c("chr5","chr14"),c("chr10","chrX"),c("chr3","chr19"),c("chr10","chr17"),c("chr13","chr17"),c("chr13","chr22"),
    #                       c("chr10","chr21"),c("chr11","chr19"),c("chr20","chr21"),c("chr4","chrX"))
    hot_chr_paire <- list(c("chr5","chr15"),c("chr5","chr22"),c("chr8","chr11"),c("chr11","chr22"),c("chr5","chr1"),c("chr15","chr22"),c("chr5","chr11"))
    for (hot_chr in hot_chr_paire ) {
      link2 <- link[c(which(link$CHROM_A %in% hot_chr & link$CHROM_B %in% hot_chr)),]
      link_data <- link2[,c(2,4,5,6,7,8,9,11)]
      colnames(link_data) <- c("ID","chr1","start1","end1","chr2","start2","end2","color")
      # link_data$ID <- paste0("ID",1:nrow(link_data))
      bed1 <- link_data[,c(2,3,4)];colnames(bed1) <- c('chr','start','end');bed1$value <- rep(0,nrow(bed1))
      bed2 <- link_data[,c(5,6,7)];colnames(bed2) <- c('chr','start','end');bed2$value <- rep(0,nrow(bed2))
      #基因组的展示
      pdf(paste0('AM_',hot_chr[1],"-",hot_chr[2],'_circos.pdf'),width=5,height=5)
      circos.initializeWithIdeogram(chromosome.index =c(hot_chr[1],hot_chr[2]),axis.labels.cex = 0.2, labels.cex = 0.5)
      #
      # circos.genomicLink(bed1, bed2,col = rand_color(nrow(bed1),transparency = 0.01))
      circos.genomicLink(bed1, bed2,col =link_data$color)
      title = Legend(at = paste(hot_chr[1],hot_chr[2],sep = "-"), type = "points", ncol=1,
                     legend_gp = gpar(col = 'transparent'),background = "transparent",
                     labels_gp = gpar(fontsize = 6),size = unit(1, "mm"))
      draw(title, x = unit(32, "mm"), y = unit(65, "mm"), just = c("top"))
      dev.off()
    }
    
  # })