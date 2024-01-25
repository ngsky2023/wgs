args<-commandArgs(T)
amp_genes<-NULL
if(file.exists('result/amp_genes.conf_95.txt')){
  amp <- read.delim('result/amp_genes.conf_95.txt',header=F,stringsAsFactors=F,row.names=NULL)
  
  for(i in 2:ncol(amp)){
      x <- amp[,i]
      genes <- x[5:length(x)]
      genes <- as.data.frame(genes)
      names(genes) <- "Gene"
      genes$chr <- x[1]
      genes$q <- x[2]
      genes$residual_q <- x[3]
      genes$peak <- x[4]
      amp_genes<-as.data.frame(rbind(amp_genes,genes))
  }
  amp_genes[amp_genes==""] <- NA
  amp_genes <- na.omit(amp_genes)
  amp_genes$Type<-"Amp"
}

del_genes<-NULL
if(file.exists('del_genes.conf_95.txt')){
  del <- read.delim('del_genes.conf_95.txt',header=F,stringsAsFactors=F,row.names=NULL)
  if(ncol(del)==2){
    for(i in 2:ncol(del)){
      x <- del[,i]
      genes <- x[5:length(x)]
      genes <- as.data.frame(genes)
      names(genes) <- "Gene"
      genes$chr <- x[1]
      genes$q <- x[2]
      genes$residual_q <- x[3]
      genes$peak <- x[4]
      del_genes<-as.data.frame(rbind(del_genes,genes))
    }
    del_genes[del_genes==""] <- NA
    del_genes <- na.omit(del_genes)
    del_genes$Type<-"Del"
  }
}

CNV <- as.data.frame(rbind(amp_genes,del_genes))
#head(CNV)
CNV$q<-as.numeric(CNV$q)
threshold<-as.numeric(args[1])
CNV <- CNV[CNV$q<threshold,]
# cosmic<-read.table("/share/work1/hanwj4457/database/cancergene/cancer_gene_census_OncTSG.txt",sep='\t',header = T,stringsAsFactors = F)
# CNV<- merge(CNV,cosmic,by=1,sort=F)

write.table(CNV,file="result/genes.conf_95.txt",quote=F,sep="\t",row.names=F)

