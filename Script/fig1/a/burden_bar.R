###add top genes on by frequency
rm(list=ls())
library(reshape2)
library(RColorBrewer)
library(plotrix)
library(scales)
##################################func####################################
Plot_func<-function(gene=gene,sample=sample,df=df,topgene=topgene,mut_col=mut_col,cex1=2,cex2=1,cex3=0.7,xdist=0,ydist=0.2){
  y<-topgene
  plot(x=nrow(sample)+1,y=topgene,xlim = c(0,nrow(sample)-1.5),ylim=c(0,topgene),cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
  for(i in as.character(gene$gene)){
    x<-0
    for(j in as.character(sample$Sample)){
      tmp<-df[which(df$gene==i & df$Sample==j),]
      tmp<-tmp[order(tmp[,"class"],decreasing=TRUE),]
      if(nrow(tmp)==0){
        points(x=x,y=y,pch=15,col="#EDEDED",cex=cex1)
      }else if(nrow(tmp)>0 & nrow(tmp)<=2){
        for(k in 1:nrow(tmp)){
          if(k==1){
            points(x=x,y=y,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[1,"type"]),"col"]),cex=cex1)
          }else if(k==2){
            points(x=x-xdist,y=y,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[2,"type"]),"col"]),cex=cex2)
          }
        }
      }else if(nrow(tmp)==3){
        for(k in 1:nrow(tmp)){
          if(k==1){
            points(x=x,y=y,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[1,"type"]),"col"]),cex=cex1)
          }else if(k==2){
            points(x=x-xdist,y=y+0.2,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[2,"type"]),"col"]),cex=cex3)
          }else if(k==3){
            points(x=x-xdist,y=y-0.2,pch=15,col=as.character(mut_col[which(mut_col[,"mut"]==tmp[3,"type"]),"col"]),cex=cex3)
          }
        }
      }
      x<-x+1
    }
    #print(nrow(sample)+1)
    #print(y)
    points(x=nrow(sample),y=y,pch=15,col="white",cex=cex1)
    y<-y-1
  }
  for(i in 0:nrow(sample)){
    points(x=i,y=0,pch=15,col="white",cex=cex1)
  }
}

###pick topgene by gene frequency
Pick_topgene<-function(df,topgene=topgene){
  df1<-df[,c("sample","sbs")]
  df1<-df1[!duplicated(df1),]
  gene<-as.data.frame(table(df1$gene))
  gene<-gene[order(gene[,2],decreasing = TRUE),]
  print(head(gene,n=50L))
  gene<-gene[c(1:topgene),]
  colnames(gene)<-c("sbs","Freq")
  return(gene)
}

###pick subset data from df by topgenes
Subdata<-function(df,gene){
  df<-merge(df,gene,by="gene",sort=FALSE)
  df<-df[,c("gene","Sample","type")]
  df$type<-as.character(df$type)
  df$type[df$type=="frameshift insertion"]<-"Frameshift indel"
  df$type[df$type=="Frameshift indel"]<-"Frameshift indel"
  df$type[df$type=="frameshift substitution"]<-"Frameshift indel"
  df$type[df$type=="nonframeshift insertion"]<-"In-frame indel"
  df$type[df$type=="InFrame indel"]<-"In-frame indel"
  df$type[df$type=="nonframeshift substitution"]<-"In-frame indel"
  df$type[df$type=="splicing"]<-"Splice-site"
  df$type[df$type=="amplification"]<-"Amplification"
  df$type[df$type=="CCDC6--RET"]<-"fusion"
  df$type[df$type=="EML4--ALK"]<-"fusion"
  df$type[df$type=="EZR--ROS1"]<-"fusion"
  df$type[df$type=="KIF5B--RET"]<-"fusion"
  
  
  
  
  df<-df[!duplicated(df),]
  
  ###indicate mutation types by numbers
  df$class[df$type=="Missense"]<-1
  df$class[df$type=="Splice-site"]<-2
  df$class[df$type=="fusion"]<-3
  df$class[df$type=="stopgain"]<-4
  df$class[df$type=="In-frame indel"]<-5
  df$class[df$type=="Frameshift indel"]<-6
  # df$class[df$type=="Amplification"]<-7
  # df$class[df$type=="deletion"]<-8
  return(df)
}

###sample order in plot
Sample_order<-function(df,gene){
  ##for genes of each sample with more than one mutation type, retain the greatest "class"
  df<-df[order(df[,"gene"],df[,"Sample"],df[,"class"]),]
  df$filter<-paste(df$gene,df$Sample,sep="-")
  df_uniq_type<-df[!duplicated(df$filter,fromLast=TRUE),]
  df_uniq_type<-df_uniq_type[,c("gene","Sample","class")]
  ##order each line by mutation type according to class number
  df_uniq_type_2<-dcast(df_uniq_type,gene~Sample)
  df_uniq_type_2<-merge(gene,df_uniq_type_2,by="gene",sort = FALSE)
  rownames(df_uniq_type_2)<-df_uniq_type_2[,1]
  df_uniq_type_2<-df_uniq_type_2[,-c(1,2)]
  df_uniq_type_2[is.na(df_uniq_type_2)]<-0
  order <- "df_uniq_type_2 <- df_uniq_type_2[,order("
  for(i in 1:nrow(df_uniq_type_2)){
    order <- paste(order,"df_uniq_type_2[",i,",],",sep = "")
  }   
  order <- paste(order,"decreasing = TRUE)]")
  print(order)
  df_uniq_type_2 <- eval(parse(text=order))
  sample<-as.data.frame(colnames(df_uniq_type_2))
  names(sample)<-"Sample"
  return(sample)
}

###frequency
Frequency<-function(df=df,gene=gene){
  freq <-NULL
  for(i in as.character(gene$gene)){
    tmp<-df[which(df$gene==i),]
    sample<-unique(tmp$Sample)
    num<-length(sample)
    freq <- c(freq,num)
  }
  return(freq)
}

###plot matrix
Matplot_function<-function(df=df){
  df[is.na(df)]<-"NN"
  num<-2
  legend_at<-c(0,1)
  legend_label<-c("No","Yes")
  for(i in 1:ncol(df)){
    unit<-sort(unique(df[,i]))
    if(length(unit)<=3){
      for(j in 1:length(unit)){
        if(unit[j]!="NN" & unit[j]!="0" & unit[j]!="1" & unit[j]!=0 & unit[j]!=1){
          df[,i][df[,i]==unit[j]]<-num
          legend_at<-c(legend_at,num)
          legend_label<-c(legend_label,unit[j])
          num<-num+1
        }else{
          next
        }
      }
    }else if(length(unit)>3){
      for(j in 1:length(unit)){
        if(unit[j]!="NN"){
          df[,i][df[,i]==unit[j]]<-num
          legend_at<-c(legend_at,num)
          legend_label<-c(legend_label,unit[j])
          num<-num+1
        }else{
          next
        }
      }
    }
  }
  df[df=="NN"]<-max(legend_at)+1
  legend_at<-c(legend_at,max(legend_at)+1)
  legend_label<-c(legend_label,"NA")
  
  df<-as.data.frame(lapply(df,as.numeric))
  df_plot<-t(df)
  col<-matrix(rep("#000000",length(df_plot)),nrow=nrow(df_plot))
  # color<-c("red","white",brewer.pal(8,"Set1")[3:7],"#FFCC99","#99CC99","#FFFFCC","#FFCCCC","#99CC99","#FF99CC","#CCFFFF","#66CCCC","#009999","#336666","#CCCCFF","#9999CC","#996699","#663366","#330033","#CC9966","#666666","#CC9999","#CCCC99","#000000","#CCCCCC")
  
  color<-c("red","white","#99CCFF","#6699CC","#006699","#336699","#003366","#FFCC99","#99CC99","#FFFFCC","#FFCCCC","#99CC99","#FF99CC","#CCFFFF","#66CCCC","#009999","#336666","#CCCCFF","#9999CC","#996699","#663366","#330033","#CC9966","#666666","#CC9999","#CCCC99","#000000","#CCCCCC")
  # color<-c(color[1:(length(legend_at)-1)],"white")
  for(i in legend_at){
    col[df_plot==i]<-color[i+1]
  }
  return(list(plot_df=df_plot,col=col,legend_color=color,legend_label=legend_label))
}



plot_clin<-function(data=clinical,clin_col=clin_col,cex=cex){
  x<-1:nrow(clinical)
  y<-1:ncol(clinical)
  plot(x=length(x),y=length(y),xlim = c(1,length(x)),ylim=c(1,length(y)),cex=0,xlab = "",ylab = "",xaxt='n',yaxt='n',frame.plot=F)
  for(i in y){
    for(j in x){
      points(x=j,y=i,pch=15,col=as.character(clin_col[match(clinical[j,i],clin_col[,"clin"]),"col2"]),cex=cex)
    }
    points(x=nrow(clinical)+1,y=i,pch=15,col="white",cex=cex)
  }
  
}

#read data
df<-read.table("mutation_burden.xls", sep="\t", quote="", header=T,stringsAsFactors = F)
df<- df[,c(1,5)]
df <- df[order(df$perMb,decreasing = T),]
names(df) <- c("sample","Mutation per megabase")
sv <- read.table("sv_filter_BND_number.xls",sep="\t", quote="", header=T,stringsAsFactors = F)
names(sv) <- c("sample","Translocation","Deletion","Duplication","Inversion")
#cnv
cnv <- read.table("cnv_type_length.xls",sep="\t",quote = "",header = T,stringsAsFactors = F)
names(cnv) <- c("sample","Deletion","Loss CN1","Copy Neutral LOH","CN2","CN3-5","CN>=6")
cnv <- cnv[,c("sample","Deletion","Loss CN1","Copy Neutral LOH","CN>=6","CN2","CN3-5")]
#WGS
wgd <- read.table("01.purity_WGD.txt",sep="\t", quote="", header=T,stringsAsFactors = F)
wgd <- wgd[,c(1,4)]
names(wgd) <- c("sample","WGD")
#clinical

clinical<-read.table("clinical_55.txt", sep="\t", quote="", header=TRUE,na.strings = "")
clinical <- clinical[,c(1,2,5,4,10,7,9)]
names(clinical) <- c("sample","name","Gender","Age","Site","Stage","Type")
clinical <- merge(clinical,wgd,by="sample")


clinical$Type[which(clinical$Type=="metastasis")] <-20
clinical$Type[which(clinical$Type=="primary")] <- 19
clinical$Type[which(clinical$Type=="NA")] <- 25



#Age

age <- median(na.omit(clinical$Age) )
print(age)
clinical$age[which(clinical$Age >=58)] <- ">58"#"70-79y"
clinical$age[which(clinical$Age <58)] <- "<58"
clinical$age[which(clinical$Age =="NA")] <- "NA"
#
clinical$Age[which(clinical$age =="NA")] <- 25#"70-79y"
clinical$Age[which(clinical$age %in% c(">58"))] <- 7
clinical$Age[which(clinical$age %in% c("<58"))] <- 4
#Stage
clinical$Stage[which(clinical$Stage=="Ib")] <- "I"
clinical$Stage[which(clinical$Stage=="IIIC")] <- "III"
clinical$Stage[which(clinical$Stage=="IIc")] <- "II"
clinical$Stage[which(clinical$Stage=="IIa")] <- "II"
clinical$Stage[which(clinical$Stage=="IIIB")] <- "III"
clinical$Stage[which(clinical$Stage=="IIb")] <- "II"
#
clinical$Stage[which(clinical$Stage=="NA")] <-25
clinical$Stage[which(clinical$Stage=="I")] <- 8#
clinical$Stage[which(clinical$Stage=="II")] <- 9#
clinical$Stage[which(clinical$Stage=="III")] <- 10#
clinical$Stage[which(clinical$Stage=="IV")] <- 11#
clinical$Stage[which(clinical$Stage=="Ⅳ")] <- 11
#Gender
clinical$Gender[which(clinical$Gender=="F")] <- 13
clinical$Gender[which(clinical$Gender=="M")] <- 14
clinical$Gender[which(clinical$Gender=="NA")] <- 25
# clinical$Gender <- factor(clinical$Gender,levels = c("F","M"))
#Site
clinical$Site[which(clinical$Site=="NA")] <- 25

clinical$Site[which(clinical$Site=="hand")] <- 15
clinical$Site[which(clinical$Site=="heel")] <- 16
clinical$Site[which(clinical$Site=="sole")] <- 17
clinical$Site[which(clinical$Site=="other")] <- 18
#wgd

clinical$WGD[which(clinical$WGD==FALSE)] <- 23
clinical$WGD[which(clinical$WGD==TRUE)] <- 24
#
############


###merge
data <- merge(df,clinical,by="sample")
data <- data[order(data$Age,data$Gender,data$Stage,decreasing=F),]
data <- data[order(data$Site,decreasing=F),]
# data <- data[order(data$Region,decreasing=F),]
data <- data[order(data$`Mutation per megabase`,decreasing=T),]
# sample_order <- as.data.frame(data[,1])
# names(sample_order) <- "sample"
##?如何让每个样本只显示4个
# data <- merge(sample_order,data,by="sample")

# sig <- merge(sample_order,sbs,by="sample",sort=F)
sample_order <- as.data.frame(data$sample)
names(sample_order) <- "sample"
# data <- data[match(sample_order$sample,data$sample),]

# freq <- NrowSums()# freq <- NULL
# for(i in 1:nrow(sig.t)){
#   num <- 0;
#   for(j in 1:ncol(sig.t)){
#     num <- num+sig.t[i,j]
#   }
#   freq <- c(freq,num)
# }
# sig.t <- as.data.frame(sig.t)
# sig.t$Freq <- freq
# sig.t <- sig.t[order(sig.t[,"Freq"],decreasing = T),]
# sig.t <- subset(sig.t,select=-Freq)

###sv
sv <- merge(sample_order,sv,by="sample",all=T)
sv <- sv[match(sample_order$sample,sv$sample),]
rownames(sv) <- sv$sample
sv <- sv[,c(-1)]

sv.t <- as.matrix(t(sv))
# sv.p <- t(t(sv.t)/colSums(sv.t))
# sv.p[is.na(sv.p)] <- 0
# svnatre.p["other",] <- (1- colSums(svnatre.p)  )
#cnv
cnv <- merge(sample_order,cnv,by="sample",all=T)
cnv <- cnv[match(sample_order$sample,cnv$sample),]
rownames(cnv) <- cnv$sample
cnv <- cnv[,c(-1)]

cnv.t <- as.matrix(t(cnv))
cnv.p <- t(t(cnv.t)/colSums(cnv.t))
cnv.p[is.na(cnv.p)] <- 0
# cnv.p["other",] <- (1- colSums(cnv.p)  )


############
clin <- data[,c("Stage","Type","Age","Gender","WGD","Site")]
names(clin) <- c("Stage","Specimen type","Age","Gender","WGD","Primary site")
for(i in 1:ncol(clin)){
  clin[,i]<-as.numeric(clin[,i])
}

class<-t(clin)
col2<-matrix(rep("#000000",length(class)),nrow=nrow(class))
col2[class==1]<-"#5371A3"
#Age"#fff0f3","#ffccd5","#ff8fa3","#ff4d6d","#c9184a","#590d22"
col2[class==2]<-"#fff0f3"
col2[class==3]<-"#ffccd5"
col2[class==4]<-"#ff8fa3"
col2[class==5]<-"#ff4d6d"
col2[class==6]<-"#c9184a"
col2[class==7]<-"#590d22"
#Stage"#EDE0D4","#E6CCB2","#DDB892","#9C6644"
col2[class==8]<-"#EDE0D4"
col2[class==9]<-"#E6CCB2"
col2[class==10]<-"#DDB892"
col2[class==11]<-"#9C6644"
#
col2[class==12]<-"lightgrey"
#sex"#fbb1bd","#bbd0ff"
col2[class==13]<-"#fbb1bd"
col2[class==14]<-"#bbd0ff"
#site
shows=c("#CC9933","#D28EFF","#00A087FF","#6633FF")
show_col(shows)
col2[class==15]<-"#CC9933" #hand
col2[class==16]<-"#D28EFF"#heel
col2[class==17]<-"#00A087FF"#sole
col2[class==18]<-"#6633FF"#other
#type"#03827f","#f94144"
col2[class==19]<-"#03827f"
col2[class==20]<-"#f94144"
#Region "#48D1CC","#996666"
col2[class==21]<-"#48D1CC"
col2[class==22]<-"#996666"
#WGD "#555555","#99CC99"
col2[class==23]<-"#555555"
col2[class==24]<-"#99CC99"
#
col2[class==25]<-"#d3d3d3"
cols <-  c("#ed1299", "#ddd53e", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92")
library(scales)
##########################################################################
pdf("china_AM_burden_clin.pdf",width=21, height=16 )
layout(matrix(c(1:8), 4, 2, byrow = TRUE),  heights=c(1.5,1.5,1.5,2.5),widths = c(3,1))
#TMB
par(mar=c(1, 11.5, 8, 3.1), xpd=TRUE)
barplot(data$`Mutation per megabase`, col="#2927c4", border="#2927c4", space=0.1, 
        names.arg=NULL,xaxt='n',yaxt='n',ylab="", xlab="", main="", cex.lab=1.4, cex.axis=2, font=2, las=1);
axis(2, at=c(0,4,8,12,16), labels=c(0,4,8,12,">=16"), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "Mutation per \nmegabase", line = 5, cex=2,font = 1)

par(mar=c(5, 1, 10, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)

###sv
col5 <-  c("#003C9D","#DDAA00","#00AA00","#FF3333")
par(mar=c(0, 11.5, 2, 3.1), xpd=TRUE)
barplot(sv.t,col =col5,border="black",space = 0.1, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=2, font=2, las=1)
axis(2, at=c(0,500,1000,1500), labels=c(0,500,1000,1500), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "SV count", line = 5, cex=2,font = 1)

par(mar=c(10, 0, 10, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=col5,legend=rownames(sv.t), title="", title.adj = 0, border=FALSE, cex=2,horiz=F,ncol = 1,bty="n")

#cnv
# col1 <-  c("#003C9D","#DDAA00","#00AA00","#FF3333", "#8000FF","#ff8000","#FF44AA","#660000","#003333","#99FFFF")
col1 <-  c("#003C9D","#DDAA00","#00AA00", "#ff8000","#8000FF","#FF3333")
par(mar=c(0, 11.5, 1, 3.1), xpd=TRUE)
barplot(cnv.p,col =col1,space = 0.1, xaxt='n',yaxt='n', ylab="", xlab="", main="", cex.lab=1.4, cex.axis=2, font=2, las=1)
axis(2, at=c(0,0.5,1), labels=c(0,50,100), las=1, cex.axis=2,font = 1);
mtext(side = 2, text = "%genome", line = 5, cex=2,font = 1)

par(mar=c(10, 0, 10, 0), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11, pch=15,  col=col1[1:nrow(cnv.p)],legend=rownames(cnv.p), title="cnv", title.adj = 0, border=FALSE, cex=2,horiz=F,ncol = 1,bty="n")
###
#clin
par(mar=c(22.5, 15, 1, 7), xpd=TRUE)
color2D.matplot(class, cellcolors=col2, border = NA, main="", xlab="", ylab="", axes=FALSE) 
axis(2,at=0.5:(nrow(class)-0.5),labels = rev(rownames(class)),las=2,cex.axis=2,font=1,line = F,tick = F)
for(i in (1:(ncol(class)-1))){
  segments(i,0,i,nrow(class),col = "white",lwd=0)
}
for(i in (1:(nrow(class)-1))){
  segments(0,i,ncol(class),i,col = "white",lwd=1)
}
# color2D.matplot(plot_df, cellcolors=col, border = "white", main="", xlab="", ylab="", axes=FALSE)
axis(1, at=c(1:length(data$name))-0.5,labels=data$name, line=FALSE, tick=FALSE, las=2, cex.axis=2,font=1,mgp=c(0, 1,1))

par(mar=c(0, 0, 1, 1), xpd=TRUE)
plot(x=5,y=8,cex=0,xlab = "",ylab = "",xaxt="n",yaxt="n",frame.plot=FALSE)
legend(3,11.5, pch=15,  col=c("#EDE0D4","#E6CCB2","#DDB892","#9C6644","#d3d3d3"),legend=c("I","II","III","IV","NA"), title="", title.adj = 0, border=FALSE, cex=2,horiz=T,bty="n")
legend(3,11, pch=15,  col=c("#03827f","#f94144"),legend=c("Primary","Metastasis"),title="", title.adj = 0, border=FALSE,  cex=2,horiz=T,bty="n")
legend(3,10.5, pch=15,  col=c("#ff8fa3","#590d22","#d3d3d3"),legend=c("<58",">=58","NA"), title="",title.adj = 0, border=FALSE, cex=2,horiz=T,bty="n")
legend(3,10, pch=15,  col=c("#fbb1bd","#bbd0ff"),legend=c("F","M"),title="", title.adj = 0, border=FALSE,  cex=2,horiz=T,bty="n")
legend(3,9.5, pch=15, col=c("#555555","#99CC99"),legend=c("non-WGD","WGD"), title="",title.adj = 0, border=FALSE, cex=2,horiz=F,ncol = 2,bty="n")
legend(3,9, pch=15, col=c("#CC9933","#D28EFF","#00A087FF","#6633FF"),legend=c("hand","heel","sole","other(foot)"), title="",title.adj = 0, border=FALSE, cex=2,horiz=F,ncol = 2,bty="n")


dev.off()
#################################

