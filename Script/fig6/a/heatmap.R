#fishertest¨¦????¨¨??????????¨¦???????¡§o?????¡ã??¦Ì??????¡°?????

library(ggplot2)
library(ggrepel)
library(corrplot)
library(reshape2)

kf<-read.table("tra_chr_join_samplecount.xls",sep = "\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)

df<-kf
df[is.na(df)]<-0
data<-NULL

get_upper_tri <- function(df){
  #df[lower.tri(df)]<- NA
  df[upper.tri(df)]<- NA
  return(df)
}
upper_tri <- get_upper_tri(df)


melted_cormat <- melt(as.matrix(upper_tri), na.rm = TRUE)
melted_cormat$Var1 <- factor(melted_cormat$Var1,levels = rev(as.character(unique(melted_cormat$Var1))))

p<-ggplot(data = melted_cormat,aes(Var2,Var1,fill=value)) + geom_tile() +
  labs(x="",y="") +
  # scale_fill_gradient2(high = "red",mid='white',midpoint = 0.085,low = "#006633")+
  # geom_text(aes(label=paste0(round((value/55*100),0),'%')))+
  scale_fill_gradient2(high = "red",mid='white',midpoint = 5,low = "#006633")+
  geom_text(aes(label=value))+
  ggtitle("No.of sample in tra")+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black",angle = 90,vjust=0.5, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9,0.4),
        legend.title = element_blank(),
        legend.text = element_text(angle = 0))
p

ggsave("tra_sample_count.pdf",width = 8,height = 8)
