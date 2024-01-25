cnv_genes <- read.table('result/CNV_region_genelist.txt',header=T,stringsAsFactors=F)
# COSMIC <- read.table("/share/work1/hanwj4457/softwore/GISTIC2/COSMIC_driver/COSMIC.driver.tier1",header=F,stringsAsFactors=F)
COSMIC <- read.table("cancer_gene_census_OncTSG_name.txt",header=F,stringsAsFactors=F)

cnv_genes <- cnv_genes[which(cnv_genes$Gene %in% COSMIC[,1]),]

geneinfo <- read.table('gene.info',header=F,stringsAsFactors=F)
geneinfo <- geneinfo[which(geneinfo[,5] %in% cnv_genes$Gene),]

indexinfo <- read.table('Markers_ZhongHuaAll_Index.txt',header=T,stringsAsFactors=F)

index <- unlist(apply(geneinfo,1,function(x){
	chr=x[1];start=as.numeric(x[2]);end=as.numeric(x[3])
	index_temp <- indexinfo[indexinfo$Chr==chr,]
	dis_strat <- abs(index_temp$Pos-start);names(dis_strat) <- index_temp$Name
	dis_end <- abs(index_temp$Pos-end);names(dis_end) <- index_temp$Name
	min <- c(dis_strat[which.min(dis_strat)],dis_end[which.min(dis_end)])
	index <- index_temp$Index[index_temp$Name==names(min)[which.min(min)]]
	return(index)
	}))
names(index) <- geneinfo[,5]

driver_index <- data.frame(chr=geneinfo[,1],Gene=geneinfo[,5],index=index[1:nrow(geneinfo)],stringsAsFactors=F)
driver_index <- merge(driver_index,cnv_genes,by='Gene')
write.table(driver_index,file="result/CNV_driver_index.txt",quote=F,sep="\t",row.names=F)
