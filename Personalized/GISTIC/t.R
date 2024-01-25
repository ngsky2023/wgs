amp <- read.delim2("result/amp_genes.conf_95.txt",header = F,quote = "",stringsAsFactors = F,sep="\t")
amp.t <- na.omit(t(amp))

write.table(amp.t,"result/amp_genes.conf_95.t.txt",sep = "\t",row.names=FALSE,col.names = FALSE,quote=FALSE)

amp <- read.delim2("result/del_genes.conf_95.txt",header = F,quote = "",stringsAsFactors = F,sep="\t")
amp.t <- na.omit(t(amp))

write.table(amp.t,"result/del_genes.conf_95.t.txt",sep = "\t",row.names=FALSE,col.names = FALSE,quote=FALSE)

