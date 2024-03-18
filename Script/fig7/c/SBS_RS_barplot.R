###add top genes on by frequency
rm(list=ls())
library(reshape2)
library(RColorBrewer)
#library(plotrix)
library(ggsci)
library(scales)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(grid)

library(cowplot)
library(ggpubr)

library(ggsignif)
##################################RS
RS_sig <- read.delim("sv_5_samples.xls",header = T,stringsAsFactors = F)
colnames(RS_sig)<-paste0("RS",c(1:(ncol(RS_sig))))
# RS_sig$Sample <- rownames(RS_sig)
# cnv_sig$Sample <- rownames(cnv_sig)
sig4 <-RS_sig
# sig4 <- sig4[match(sample_order1$sample,rownames(sig4)),]
# rownames(sig4) <- sig4[,1]
# sig4 <- sig4[,-1]
sig4 <- as.matrix(t(sig4))
result8 <- pheatmap(sig4,
                    scale="column",cluster_rows = TRUE,
                    cutree_cols=1,cutree_rows = 1,show_colnames = TRUE,
                    color = colorRampPalette(colors = c("#fad2e1","white","#EE7E32"))(20))#"#fff0f3","#c9184a"
sample_order7 <- as.data.frame(colnames(sig4)[result8$tree_col$order])
colnames(sample_order7) <- "sample"
sample_order <- sample_order7
# write.table(file = "sample_order.xls",sample_order,sep = "\t",quote = F)
# 
############
p <- as.data.frame(cutree(result8$tree_col,k=6))
colnames(p) <- c("cluster")
p$sample <- rownames(p)
p <- p[match(sample_order$sample,p$sample),]
cluster <- p
print(cluster$cluster)
cluster$cluster[which(cluster$cluster %in% c(2))] <- "A"
cluster$cluster[which(cluster$cluster %in% c(4))] <- "B"
cluster$cluster[which(cluster$cluster %in% c(6))] <- "C"
cluster$cluster[which(cluster$cluster %in% c(1))] <- "D"
cluster$cluster[which(cluster$cluster %in% c(5))] <- "E"
cluster$cluster[which(cluster$cluster %in% c(3))] <- "F"
RS_cluster <- as.data.frame(cluster[,c("sample","cluster")])
names(RS_cluster) <- c("sample","RS_cluster")
write.table(RS_cluster,file = "RS_cluster.xls",row.names = F,quote = F,sep = "\t")
print(RS_cluster)

annotation_col = as.data.frame(RS_cluster[,c("RS_cluster")])  #增加Time，CellType分组信息
names(annotation_col) <- c("RS_cluster")
ann_colors = list(RS_cluster=c(A="#E64B35FF",B="#F39B7FFF",C="#4DBBD5FF",D="#00A087FF",E="#03827f",F="#f94144"))
rownames(annotation_col) <- RS_cluster$sample
resultRs <- pheatmap(sig4,annotation_col = annotation_col,annotation_colors = ann_colors,
                     scale="column",cluster_rows = TRUE,
                     cutree_cols=1,cutree_rows = 1,show_colnames = TRUE)
ggsave("RS_cluster.pdf",resultRs,height = 5,width = 10)
##################################SBS####################################
sbs_sig <- read.table("Decomposed_Solution_Activities_SBS96.txt",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(sbs_sig)[1] <- "sample"
sbs_sig$drug_mutagenesis <- sbs_sig$SBS17a+sbs_sig$SBS17b+sbs_sig$SBS31+sbs_sig$SBS32+sbs_sig$SBS35
#Unknown
# sbs_sig$Unknown <- sbs_sig$SBS33+sbs_sig$SBS34+sbs_sig$SBS19+sbs_sig$SBS54 #新的结果没有sbs19
sbs_sig$Unknown <- sbs_sig$SBS33+sbs_sig$SBS34+sbs_sig$SBS54
#
sbs_sig <- sbs_sig[,c("sample","SBS7a","SBS7b","SBS7c","SBS7d","SBS38","SBS2",
                      "SBS13","drug_mutagenesis","SBS9","Unknown","SBS36","SBS1","SBS5","SBS40")]

# #UV_light_exposure(SBS7a,7b,7c,7d)
# sbs_sig$UV_light_exposure <- sbs_sig$SBS7a+sbs_sig$SBS7b+sbs_sig$SBS7c+sbs_sig$SBS7d
# ##UV_light_exposure_indirect_effect(SBS38)
# sbs_sig$UV_light_exposure_indirect_effect <- sbs_sig$SBS38
# # #
# sbs_sig$aging <- sbs_sig$SBS1+sbs_sig$SBS5+sbs_sig$SBS40
# 
# 
# #APOBEC_activity(SBS2 and SBS13)
# sbs_sig$APOBEC_activity <- sbs_sig$SBS2+sbs_sig$SBS13
# # #Polymerase eta somatic hypermutation
# sbs_sig$Polymerase <- sbs_sig$SBS9
# # #BER deficiency
# sbs_sig$BER_deficiency <- sbs_sig$SBS36
# ## drug mutagenesis
# sbs_sig$drug_mutagenesis <- sbs_sig$SBS17a+sbs_sig$SBS17b+sbs_sig$SBS31+sbs_sig$SBS32+sbs_sig$SBS35
# #Unknown
# sbs_sig$Unknown <- sbs_sig$SBS33+sbs_sig$SBS34++sbs_sig$SBS54


# sbs_sig <- sbs_sig[,c("sample","UV_light_exposure","UV_light_exposure_indirect_effect","APOBEC_activity","Polymerase","BER_deficiency","drug_mutagenesis","Unknown","aging")]

#计算每个样本突变的百分比
sbs_sig[2:ncol(sbs_sig)] <- sbs_sig[2:ncol(sbs_sig)]/rowSums(sbs_sig[2:ncol(sbs_sig)])
sig1 <-sbs_sig
rownames(sig1) <- sig1[,1]     
sig1 <- sig1[,-1]
sig1 <- as.matrix(t(sig1))
sig1_scale <- round(t(apply(sig1, 1, scale)),2)
colnames(sig1_scale) <- colnames(sig1)

# sig8 <- sig1_scale[,match(sample_order7$sample,colnames(sig1_scale))]
# 
# 
# resultSBS <- pheatmap(sig8,scale="column",treeheight_row=10,treeheight_col=3,cluster_rows = FALSE,cluster_cols = FALSE,
#                       cutree_cols=1,cutree_rows = 1,show_colnames = TRUE,color = colorRampPalette(colors = c("#fad2e1","white","#EE7E32"))(20))
# ggsave("RS_SBS.pdf",resultSBS,height = 3,width = 8)
data <- merge(RS_cluster,sbs_sig,by='sample')


data.m <- melt(data = data)
names(data.m)[3:4] <- c("sbs","num")
# data.m <- aggregate(x=data.m$count, by=list(data.m$region,data.m$sig),sum)
col1 <-c("#800f2f","#a4133c","#c9184a","#ff4d6d",
         "#66FF66","#CCCCFF","#7744FF","#FF8C00","#fdffb6",
         "#6c757d","#caffbf","#adb5bd","#CC99CC","#993399")
#############
variables <- c("RS_cluster")
p_list <- lapply(1:length(variables),function(x){
  clin <- variables[x]
  data <- data.m[!is.na(data.m[,clin]),]
  plotdata <- data %>% group_by(RS_cluster,sbs) %>% summarise(sum(num))####修改
  ######################
  plotdata <- as.data.frame(plotdata)
  colnames(plotdata) <- c(clin,"sbs","num")
  p <- ggplot(plotdata,aes(x=as.character(plotdata[,clin]),y=num,fill="sbs"))+
    geom_bar(aes(fill=sbs),stat='identity',position = "fill")+
    # geom_text(aes(fill=sbs,label = paste0(format(data.m$Freq*100,digits = 3),"%"),size=8),position=position_fill(vjust =0.5),size = 2.7)+
    scale_fill_manual(values =col1 )+
    theme_bw()+
    labs(y="The percentage of SBS",x="",fill="sbs")+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line=element_line(size=1),
          axis.text = element_text(color='black',size=8))
})
p_leg <- as_ggplot(get_legend(p_list[[1]]))
p_list <- lapply(p_list,function(x) x<-x+theme(legend.position = "none") )

pdf("AM_RS_cluster_SBS_barplot.pdf",width=length(variables)*4+2,height=4)
p_list[[length(variables)+1]] <- p_leg
plot_grid(plotlist = p_list,ncol=length(variables)+1,rel_widths = c(rep(1,length(variables)-1),1.8,1))
dev.off()


