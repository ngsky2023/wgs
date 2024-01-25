# args<-commandArgs(T)

# amp_level <- read.delim(args[1],header=T,stringsAsFactors=F)
amp_level <- read.delim("result/broad_significance_results.txt",header=T,stringsAsFactors=F)

amp_border_col <- sapply(amp_level$`Amp.q.value`,function(x){
	if(x<0.1){
		col <- 'black'
	} else {
		col <- NA
	}
})
del_border_col <- sapply(amp_level$`Del.q.value`,function(x){
	if(x<0.1){
		col <- 'black'
	} else {
		col <- NA
	}
})

pdf("result/Amp_level.pdf",width=8,height=3.5)
layout(matrix(c(1,2,3),ncol=1),heights=c(1,0.1,1))
par(mar=c(0.5,4,2,1),xpd=TRUE,lwd=2)
barplot(rep(1,nrow(amp_level)),col=rep(c('grey90',"grey90",'white','white'),11),space=0,xaxt='n',yaxt='n',ylim=c(0,1),border=NA)
barplot(amp_level$`Amp.frequency`,col='#FF6666',border=NA,mgp=c(1,0.5,0.5),las=2,xaxt='n',add=TRUE,space=0)
barplot(amp_level$`Amp.frequency`,col='transparent',border=amp_border_col,mgp=c(1,0.5,0.5),las=2,xaxt='n',add=TRUE,space=0,yaxt='n')
mtext("Amp Frequency",side=2,line=2.5,cex=0.7)

par(mar=c(0,4,0,1),xpd=TRUE)
x <- barplot(rep(1,nrow(amp_level)),col=rep(c('grey90','white'),22),space=0.05,xaxt='n',yaxt='n')
text(x=x-0.05,y=0.5,amp_level$Arm,xpd=TRUE,cex=0.9)

par(mar=c(3,4,0.5,1),xpd=TRUE,lwd=2)
barplot(rep(-1,nrow(amp_level)),col=rep(c('grey90',"grey90",'white','white'),11),space=0,xaxt='n',yaxt='n',ylim=c(-1,0),border=NA)
barplot(-amp_level$`Del.frequency`,col='#0066CC',border=NA,mgp=c(1,0.5,0.5),las=2,add=TRUE,space=0,yaxt='n')
barplot(-amp_level$`Del.frequency`,col='transparent',border=del_border_col,mgp=c(1,0.5,0.5),las=2,xaxt='n',add=TRUE,space=0,yaxt='n')
axis(side=2,at=seq(-1,0,0.2),labels=seq(1,0,-0.2),cex=0.6,las=2,mgp=c(1,0.5,0.5))
mtext("Del Frequency",side=2,line=2.5,cex=0.7)

dev.off()

