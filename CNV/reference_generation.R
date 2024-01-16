
Args <- commandArgs(TRUE)
if (length(Args) != 5) {
        cat("usage: <infile prefix> <outdir> <nBlock> <nStep> <library>\n", file=stderr())
        quit(status=1)
}

getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

program_dir = getScriptPath()
parameter = list()

input_dir = Args[1]
parameter$output_reference = Args[2]
CHRs_block = paste(program_dir,'20K.CHRs.block',sep="/")
parameter$nBlock = as.numeric(Args[3])
parameter$nStep = as.numeric(Args[4])
parameter$library = as.logical(Args[5])

parameter$library = TRUE

parameter$Ncores = 3
parameter$file_GC = list.files(input_dir,".20K.GC")
parameter$file_20K = gsub('.GC','',parameter$file_GC)


parameter$no_na_ratio = 0.4 
parameter$chromosome_ratio = c(8.3728,9.0666,7.5199,7.1003,6.7485,6.3726,5.5302,5.4967,4.1752,4.8609,4.9726,4.8822,3.7005,3.328,2.9204,
                        2.77,2.6975,2.9284,1.7121,2.3308,1.295,1.2189,2.5806,0.1822,5.1116,0.0043)

library('parallel')

read_data <- function(filepath_20K,parameter){
  
  library(stringr)
  
  remove_low_coverage <- function(da_rd,da_gc,da_li,nBlock,nStep){
    
    if (nBlock != 1) {
	.libPaths("/share/work1/wangrr/local/Rlib/")
      library(zoo)
      da_rd <- t(rollapply(t(da_rd), width=nBlock, by=nStep,FUN=sum, partial=T, align="left"))
      da_gc <- t(rollapply(t(da_gc), width=nBlock, by=nStep, FUN=sum, partial=T, align="left"))
      da_li <- t(rollapply(t(da_li), width=nBlock, by=nStep, FUN=sum, partial=T, align="left"))
    }
    
    valide_cut = limit_cutoff * (window_len * (2.1/3)) * nBlock
    usefull <- da_li <= valide_cut
    
    da_rd[da_rd==0] <- NA
    da_rd[usefull] <- NA           
    da_gc <- round(da_gc/(da_rd*36)*100, 3)
    da_gc[!is.na(da_gc) & (da_gc > 70 | da_gc < 25)] <- NA  
    da_rd[is.na(da_gc)] <- NA  
    
    return(list(da_rd = da_rd, da_gc = da_gc))
  }
  
  smooth_spline_Fit <- function(rc, gc, deter){
    
    
    rc2 <- as.vector(rc)
    rc2 <- rc[!is.na(rc2)]
    gc2 <- as.vector(gc)
    gc2 <- gc[!is.na(gc2)]
    
    rc <- as.vector(rc[deter,])
    rc <- rc[!is.na(rc)]
    gc <- as.vector(gc[deter,])
    gc <- gc[!is.na(gc)]
    
    gc = gc[rc < 200 * nStep]
    rc = rc[rc < 200 * nStep]
    
    ## normal smooth spline
    
    smooth_model = smooth.spline(rc ~ gc,spar = 0.8)
    
    return(smooth_model)
  }
  
  smooth_spline_gc_correction <- function(da_rd,da_gc,deter){
    
    RCgcCorrection.tatol.reads <- sum(da_rd)
    normal_count <- sum(as.numeric(da_rd))  
    
    MedianRC <- median(da_rd[1:22,], na.rm=T)
    smooth_model <- smooth_spline_Fit(da_rd, da_gc,deter)
    gc_vector = as.vector(da_gc)
    Fit_vector <- predict(smooth_model, gc_vector[!is.na(gc_vector)])$y
    gc_vector[!is.na(gc_vector)] = Fit_vector
    
    #RCgcCorrection = da_rd 
    
    RCgcCorrection <- ( MedianRC / gc_vector ) * da_rd
    RCgcCorrection[!is.na(RCgcCorrection) & RCgcCorrection <= 0] <- NA
    
    return(RCgcCorrection)
  }
  
  window_len <- 20000
  lib_reads_count <- 2176351405
  limit_cutoff <- 0.2
  nBlock = parameter$nBlock
  nStep = parameter$nStep
  
  filepath_20GC <- str_replace(filepath_20K,"20K", "20K.GC")
  data = as.matrix(read.table(filepath_20K, row.names=1, stringsAsFactors=F, header=F))
  dataGC = as.matrix(read.table(filepath_20GC, row.names=1, stringsAsFactors=F, header=F))
  sampleName = strsplit(basename(filepath_20K),split='\\.')[[1]][1]
  
  
  ratio = round(sum(data[24,],na.rm=T)/sum(data[1:22,],na.rm=T)*100,3)
  sex = ifelse(ratio>0.08,'M','F')
  
  if(parameter$library){
    ori = remove_low_coverage(data,dataGC,parameter$da_li,nBlock,nStep)
  }else{
    parameter$da_li = matrix(rep(20000,12500*24),nrow=24)
    cutoff = floor(mean(data[data>0])*0.1)
    parameter$da_li[data <= cutoff] = 0
    ori = remove_low_coverage(data,dataGC,parameter$da_li,nBlock,nStep)
  }
  
  data = ori$da_rd
  dataGC = ori$da_gc
  
  if(sex == 'M'){
    chromosome_ratio = parameter$chromosome_ratio[1:24]
  }else{
    chromosome_ratio = c(parameter$chromosome_ratio[1:22],parameter$chromosome_ratio[25:26])
  }
  
  bili = 5000000 / sum(data[1:22,],na.rm=T)
  data = data * bili
  
  chr_ratio = apply(data,1,sum,na.rm=T) / sum(data[1:22,],na.rm=T) * 100
  index = 1:22
  deter = index[abs(chr_ratio[1:22] - chromosome_ratio[1:22]) / chromosome_ratio[1:22] < 0.4]
  data = smooth_spline_gc_correction(data,dataGC,deter)
  
  bili = 5000000 / sum(data[1:22,],na.rm=T)
  data = data * bili
  
  ratio = apply(data,1,sum,na.rm=T)/5000000 * 100
  deter = sapply(1:24,function(x) ifelse(abs(ratio[x]-chromosome_ratio[x])/chromosome_ratio[x]<0.2,TRUE,FALSE))
  
  return(list(data = data, sex = sex,deter = deter,ratio = ratio))
}

select_use_bin <- function(alldata){
  
  reference_all = array(0,dim=c(26,12500/parameter$nStep))
  
  for(i in 1:26){
    for(j in 1:12500/parameter$nStep){
      
      if(i == 23 | i == 24) {
        temp = sapply(alldata,function(x) return(ifelse(x$sex=='M',x$data[i,j],NA)))
      }else if(i == 25 | i == 26) {
        temp = sapply(alldata,function(x) return(ifelse(x$sex=='F',x$data[i-2,j],NA)))
      }else{
        temp = sapply(alldata,function(x) return(x$data[i,j]))
      }
      reference_all[i,j] = sum(temp,na.rm=T)
    }
  }
  
  reference_all[reference_all < 1] = NA
  cut_count = apply(reference_all,1,mean,na.rm=T) * 0.2
  
}

remove_outlier <- function(count){
  
  count = count[!is.na(count)]
  
  if(length(count) <= 1){return(count)}
  
  bound = 1.96
  index = 1:length(count)
  out_index = index[count > mean(count) + sd(count) * bound | count < mean(count) - sd(count) * bound]
  if(length(out_index) == 0){return(count)}
  else{
    count = count[-out_index]
    count = remove_outlier(count)
  }
  return(count)
}

generate_reference <- function(data,parameter,deter,sex){
  
  reference_all = array(0,dim=c(26*2,12500/parameter$nStep))
  
  for(i in 1:26){
    for(j in 1:12500/parameter$nStep){
      if(i == 23 | i == 24) {
        temp = data[i,j,sex == 'M' & deter[,i] == TRUE]
      }else if(i == 25 | i == 26){
        temp = data[i-2,j,sex == 'F' & deter[,i-2] == TRUE]
      }else{
        temp = data[i,j, deter[,i] == TRUE]
      }
      temp = remove_outlier(temp)
      if(i == 23 | i == 24) {
        reference_all[i,j]=ifelse(length(temp)>parameter$no_na_ratio*length(sex[sex == 'M']),mean(temp),0)
        reference_all[i+26,j]=ifelse(reference_all[i,j]==0,NA,sd(temp)/reference_all[i,j])
      }else if(i == 25 | i == 26) {
        reference_all[i,j]=ifelse(length(temp)>parameter$no_na_ratio*length(sex[sex == 'F']),mean(temp),0)
        reference_all[i+26,j]=ifelse(reference_all[i,j]==0,NA,sd(temp)/reference_all[i,j])
      }else{
        reference_all[i,j]=ifelse(length(temp)>parameter$no_na_ratio*length(sex),mean(temp),0)
        reference_all[i+26,j]=ifelse(reference_all[i,j]==0,NA,sd(temp)/reference_all[i,j])
      }
    }
  }
  
  bili = 5000000 / sum(reference_all[1:22,],na.rm=T)
  reference_all = reference_all * bili
  
  reference_all = round(reference_all,5)
  
  write.table(reference_all,parameter$output_reference,sep='\t',quote=FALSE,col.names = F,row.names = T)
  
  return(reference_all)
}

filepath_20K <- as.list(paste(input_dir,parameter$file_20K,sep="/"))
if(parameter$library){
  parameter$da_li <- t(as.matrix(read.table(CHRs_block, stringsAsFactors=F, header=F)))
}


print('[read data ...]')

cluster <- makeCluster(parameter$Ncores,"PSOCK")
alldata = clusterApplyLB(cluster,filepath_20K,read_data,parameter)
stopCluster(cluster)

#for(i in 1:length(filepath_20K)){
#  print(i)
#  read_data(filepath_20K[[i]],parameter,da_li)
#}

data = array(0,dim = c(24,12500/parameter$nStep,length(alldata)))
for(i in 1:length(alldata)){
  data[,,i] = alldata[[i]]$data
}
sex = sapply(alldata,function(x) x$sex)
deter = t(sapply(alldata,function(x) x$deter))
ratio = t(sapply(alldata,function(x) x$ratio))

print('[Calculate reference ...]')

a = generate_reference(data,parameter,deter,sex)



