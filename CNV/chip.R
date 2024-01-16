
GD_score_cutoff = 15
GD_length_cutoff = 200
SW_cutoff = 5
SW_break_num = 20
merge_cutoff = 0.2
bound = 0.25


Length_chr = c(
  249250621, 243199373, 198022430, 191154276, 180915260,
  171115067, 159138663, 146364022, 141213431, 135534747,
  135006516, 133851895, 115169878, 107349540, 102531392,
  90354753, 81195210, 78077248, 59128983, 63025520,
  48129895, 51304566, 155270560, 59373566)
centromereS = c(
  121535434,92326171,90504854,49660117,46405641,
  58830166,58054331,43838887,47367679,39254935,
  51644205,34856694,16000000,16000000,17000000,
  35335801,22263006,15460898,24681782,26369569,
  11288129,13000000,58632012,10104553)
centromereE = c(
  124535434,95326171,93504854,52660117,49405641,
  61830166,61054331,46838887,50367679,42254935,
  54644205,37856694,19000000,19000000,20000000,
  38335801,25263006,18460898,27681782,29369569,
  14288129,16000000,61632012,13104553)

absolute_error <- function(breakpoint,fd){
  # calculate absolute error in GD
  
  if(length(fd) == 1){return(0)}
  
  if(length(breakpoint) == 0){return(sum(abs(fd - median(fd))))
  }else{
    
    dimension = length(fd)
    mean_sequence = rep(0,dimension)
    index_s = c(1,breakpoint+1)
    index_e = c(breakpoint,dimension)
    
    for(i in 1:length(index_s)){
      mean_sequence[index_s[i]:index_e[i]] = median(fd[index_s[i]:index_e[i]])
    }
    
    return(sum(abs(fd - mean_sequence)))
  }
}

gradient_descent_recursive <- function(fd){
  # GD main flow
  
  fast_gradient <- function(fd,index_s,index_e){
    
    temp = sort(fd,index.return = TRUE)
    fd_sort = temp$x
    fd_index = temp$ix
    fd_right = fd_sort
    fd_left = rep(NA,length(fd_right))
    ori_index = sort(fd_index,index.return = TRUE)$ix
    
    left_median = rep(0,length(fd)-1)
    right_median = rep(0,length(fd)-1)
    
    ori_error = sum(abs(fd - median(fd)))
    gradient_error = rep(0,length(fd)-1)
    
    for(i in 1:(length(fd)-1)){
      
      fd_left[ori_index[i]] = fd[i]
      fd_right[ori_index[i]] = NA
      
      temp = fd_left[!is.na(fd_left)]
      if(length(temp) %% 2 == 1){
        left_median[i] = temp[(length(temp)+1)/2]
      }else{
        left_median[i] = (temp[length(temp)/2]+temp[length(temp)/2+1])/2
      }
      
      temp = fd_right[!is.na(fd_right)]
      if(length(temp) %% 2 == 1){
        right_median[i] = temp[(length(temp)+1)/2]
      }else{
        right_median[i] = (temp[length(temp)/2]+temp[length(temp)/2+1])/2
      }
      
      gradient_error[i] = sum(abs(fd[1:i]-left_median[i])) + sum(abs(fd[(i+1):length(fd)]-right_median[i]))
      
    }
    
    diff_error = ori_error - gradient_error
    grad = max(diff_error)
    breakpoint = which.max(diff_error)
    
    return(list(seg_error = ori_error, gradient = grad,breakpoint = breakpoint+index_s-1))
  }
  
  
  dimension = length(fd)
  breakpoint = breakpoint_old = c()
  
  index_s = 1
  index_e = length(fd)
  
  iter_times = 0
  seg_active = 1      # if seg_active == 0, then do not apply GD in this segment
  
  while(TRUE){
    
    iter_times = iter_times + 1;
    seg_gradient = seg_error = rep(0,length(index_s))
    breakpoint_old = breakpoint
    
    for(i in 1:(length(breakpoint_old)+1)){
      
      if(seg_active[i] == 0){next}
      
      #seg_info = fast_gradient(fd[index_s[i]:index_e[i]],index_s[i],index_e[i])
      seg_info = .Call("fast_gradient_so",fd[index_s[i]:index_e[i]],index_s[i],index_e[i]);
      chang = index_e[i] - index_s[i] + 1
      seg_gradient[i] = seg_info$gradient
      
      if(seg_gradient[i] < max(GD_score_cutoff * 0.7^(iter_times-1),0.3*GD_score_cutoff)){
        seg_active[i] = 0
      }else{
        position = seg_info$breakpoint
        breakpoint = c(breakpoint,position)
        index_s = c(index_s,position+1)
        index_e = c(index_e,position)
        seg_active = c(seg_active,1)
      }
    }
    
    if(!is.null(breakpoint)){breakpoint = sort(breakpoint)}
    this_order = order(index_s)
    index_s = sort(index_s)
    index_e = sort(index_e)
    seg_active = seg_active[this_order]
    seg_active[(index_e-index_s+1)<=GD_length_cutoff] = 0
    
    if(sum(seg_active) < 0.5){break}
    if(iter_times > 100){break}
  }
  
  segments = cal_seginfo(fd,index_s,index_e) 
  
  return(segments)
}

find_contig <- function(temp_fd,s,e){
  # find segment in SW
  
  if(s == e){return(list(index_s = s,index_e = e))
  }else{
    
    ave_fd = median(temp_fd)
    dimension = length(temp_fd)
    
    wucha1 = round(abs(temp_fd-ave_fd)-abs(temp_fd-ave_fd+bound),5) / bound
    wucha3 = round(abs(temp_fd-ave_fd)-abs(temp_fd-ave_fd-bound),5) / bound
    
    breakpoint1 = find_position(wucha1)
    breakpoint3 = find_position(wucha3)
    breakpoint = c(breakpoint1,breakpoint3)
    
    if(is.null(breakpoint)){
      return(list(index_s = s,index_e = e))
    }else{
      breakpoint = sort(breakpoint)
      return(list(index_s = c(s,breakpoint + s),index_e = c(breakpoint + s - 1,e)))
    }
  }
}

find_position <- function(wucha){
  # calculate region score in SW
  
  dimension = length(wucha)
  
  wucha[wucha > 0] = abs(wucha[wucha > 0])
  wucha[wucha < 0] = -abs(wucha[wucha < 0])
  
  wucha = wucha
  breakpoint = c()
  
  i = 1
  while(i <= dimension-1){
    match_array = rep(0,dimension)
    match_array[i] = wucha[i]
    for(j in (i+1):length(wucha)){
      match_array[j] = ifelse(match_array[j-1] < 0,break,match_array[j-1] + wucha[j])
      if(j > SW_break_num){
        deter = sapply(1:(SW_break_num),function(x) match_array[j-x+1]<match_array[j-x])
        if(length(deter[deter==FALSE])==0){break}
      }
    }
    if(max(match_array) < SW_cutoff){
      i = i + 1
    }else{
      position = which.max(match_array)
      breakpoint = c(breakpoint,i-1,position)
      i = position + 1
    }
  }
  breakpoint = breakpoint[breakpoint > 0 & breakpoint < dimension]
  return(breakpoint)
}

match_sequence <- function(gradient_list,fd){
  # main process in SW
  
  index_s_new = c()
  index_e_new = c()
  
  index_s = gradient_list$index_s
  index_e = gradient_list$index_e
  
  for(i in 1:length(index_s)){
    temp_fd = fd[index_s[i]:index_e[i]]
    index_temp = find_contig(temp_fd,index_s[i],index_e[i])
    index_s_new = c(index_s_new, index_temp$index_s)
    index_e_new = c(index_e_new, index_temp$index_e)
  }
  
  segments = cal_seginfo(fd,index_s_new,index_e_new)
  
  return(segments)
}

cal_seginfo <- function(fd,index_s,index_e){
  # calculate segment information
  # include start index (index_s), end index (index_e), copy number (copy), segment lenght(segLen) 
  
  copy = sapply(1:length(index_s),function(x) median(fd[index_s[x]:index_e[x]]))
  segLen = sapply(1:length(index_s),function(x) (index_e[x] - index_s[x] + 1))
  
  change = c()
  if(length(copy) == 1){change = NULL
  }else{change = sapply(1:(length(index_s)-1),function(x) abs(copy[x+1]-copy[x]))}
  
  return(list(index_s = index_s, index_e = index_e, copy = copy, segLen = segLen, change = change))
}

combine_segments <- function(segments,fd){
  # merge two segment when copy number different < merge cutoff
  
  index_s = segments$index_s
  index_e = segments$index_e
  copy = segments$copy
  change = segments$change
  
  if(length(copy) == 1){return(segments)}
  if(min(change) > merge_cutoff){return(segments)}
  
  com_point = which.min(change)
  
  return(combine_segments(cal_seginfo(fd,index_s[-(com_point+1)],index_e[-com_point]),fd))
  
}

Call_cnv_GD = function(ldata,i){
  # main workflow for GD algorithm
  
  print(paste('Chromosome:',i))
  
  gradient_list = gradient_descent_recursive(ldata)    # GD
  sequence_list = match_sequence(gradient_list,ldata)   # sw
  gradient_merge = combine_segments(sequence_list,ldata)    # 
  
  index_s = gradient_merge$index_s
  index_e = gradient_merge$index_e
  copy = gradient_merge$copy
  EstiCopy = 1:length(ldata)
  for(j in 1:length(copy)){EstiCopy[index_s[j]:index_e[j]] = copy[j]}
  
  return(EstiCopy)
}

library('zoo')
library('Rcpp')
library('rmarkdown')

dyn.load("F:/kenuoan/MyRscript/v4.2_tumor/HeapMedian.dll");

input_dir = "F:/kenuoan/tumor_data/chip_test"
cytoband_file = "F:/kenuoan/MyRscript/v4.2/cytoBand_hg19_ucsc.txt"
html_plot = "F:/kenuoan/MyRscript/v4.2_tumor/tumor_cnv_html_plot.Rmd"
files = list.files(input_dir,'.lrr.txt')
filepath <- as.list(paste(input_dir,files,sep="/"))

ideoCyto = read.table(cytoband_file, header=F, stringsAsFactors=F)

out_all_path = paste(input_dir,'/','out_all.txt',sep='')
cat(file=out_all_path,paste('sampleName','chromosome','spos','epos','copy','cnv_length',sep='\t'),end='\n')

chromosome = c(1:22,'X','Y')

for(n in 1:length(filepath)){
  
  print(paste('Begin Sample: ',n))
  
  sampleName = strsplit(basename(filepath[[n]]),split='\\.')[[1]][1]
  out_path = paste(input_dir,'/',sampleName,'.cnv.txt',sep='')
  cat(file=out_path,paste(sampleName,'M','NA',sep='\t'),end='\n')
  
  data = read.table(filepath[[n]],header=TRUE)
  
  copy_all = c()
  EstiCopy_all = c()
  chr_all = c()
  position_all = c()
  
  for(i in 1:24){
    
    ldata = data[data[,1] == chromosome[i] ,3]
    X = data[data[,1] == chromosome[i],2]
    if(i == 23 | i == 24){
      copy = 10^(ldata)
      param = 1
    }else{
      copy = 10^(ldata) * 2
      param = 2
    }
    
    EstiCopy = Call_cnv_GD(ldata,i)
    plot(X,copy,cex=0.25)
    
    if(i == 23 | i == 24){
      EstiCopy = 10^(EstiCopy)
    }else{
      EstiCopy = 10^(EstiCopy) * 2
    }
    
    copy_all = c(copy_all,copy)
    EstiCopy_all = c(EstiCopy_all,EstiCopy)
    chr_all = c(chr_all,rep(chromosome[i],length(copy)))
    position_all = c(position_all,X)
    
    index_s = 1
    index_e = c()
    spos = X[1]
    epos = c()
    copy = EstiCopy[1]
    
    num = 1
    for(j in 1:(length(EstiCopy)-1)){
      if(EstiCopy[j] != EstiCopy[j+1]){
        spos = c(spos,X[j])
        epos = c(epos,X[j+1])
        index_e = c(index_e,j)
        index_s = c(index_s,j+1)
        copy = c(copy,EstiCopy[j+1])
        num = j + 1
      }
    }
    index_e = c(index_e,length(X))
    epos = c(epos,X[length(X)])
    
    sentence = c()
    for(j in 1:length(copy)){
      if(abs(copy[j]-param) > 0.05){
        sentence = c(sentence,paste('>GDSW','\t','chr',chromosome[i],':',spos[j],'-',epos[j],'\t',
                                    round(copy[j],4),'\t', index_e[j]-index_s[j]+1,
                                    '[', index_s[j], '-', index_e[j], '(', 0, ' NA', ']',sep=''))
        cat(file=out_path,paste(sentence[length(sentence)],sep='\t'),append=T,end='\n')
        chang = round((index_e[j]-index_s[j])/1e6,3)
        if(chang > 0.1){
          cat(file=out_all_path,paste(sampleName,sentence[length(sentence)],sep='\t'),append=T,end='\n')
        }
      }
    }
  }
  
  info = list()
  info$position = position_all
  info$EstiCopy = EstiCopy_all
  info$obsCopy = copy_all
  info$chr = chr_all
  info$chromosome = chromosome
  info$sampleName = sampleName
  info$sex = 'M'
  info$ideoCyto = ideoCyto
  info$Length_chr = Length_chr
  info$centromereS = centromereS
  info$centromereE = centromereE
  rmarkdown::render(html_plot,output_file = paste(input_dir,'/',sampleName,'.cnv.html',sep=''),params = info,intermediates_dir = tempdir(),
                    html_document(fig_width = 18,fig_height = 4.5))
  
  
}




#X = data[data[,1]=="5",2]
#Y = data[data[,1]=="5",3]
#fd = 2^Y
#copy = 10^(Y)*2

#plot(X,copy,ylim =c(0,6),cex=0.25)

#plot(X,copy,cex=0.25)


#plot(X,EstiCopy,cex=0.25)
