args = commandArgs(trailingOnly=TRUE)
infile <- args[1]
pdffile <- args[2]

data <- scan(infile)
pdf(pdffile,width=8)
maxins <- 600
if(range(data)[2] > maxins){
    	plot(table(data), xlab="Insert Size", ylab="Reads Number",xlim=c(0,maxins))
    	data <- data[data<=maxins]
		hist(data, xlim=c(0,maxins))
}else{
    	plot(table(data), xlab="Insert Size", ylab="Reads Number")
    	data <- data[data<=maxins]
		hist(data)
}
dev.off()
