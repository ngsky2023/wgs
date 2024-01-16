.libPaths("/share/work1/wangrr/local/Rlib/")

Args <- commandArgs(TRUE)
if (length(Args) != 5) {
        cat("usage: <infile> <outdir> <nBlock> <nStep> <reference> <IT_version>\n", file=stderr())
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

parameter = list()

## set the parameters
input_dir = Args[1]
output_dir = Args[2]
parameter$nBlock = as.numeric(Args[3])
parameter$nStep = as.numeric(Args[4])
parameter$reference_file = Args[5]
program_dir = getScriptPath()


parameter$IT_version = FALSE



Rscript_lsb = ""
cytoband_file = "cytoBand_hg19_ucsc.txt";
preProcess = "preProcess"
Callcnv_dir = "callcnv"
result_dir = "result"
ratio = "chromosome_ratio"

### parameter

parameter$cytoband_file = paste(program_dir,cytoband_file,sep="/")
parameter$HeapMedian = paste(program_dir,'HeapMedian.so',sep="/")
parameter$html_plot = paste(program_dir,'cnv_html_plot.Rmd',sep="/")
parameter$output_copy = paste(output_dir,"chromosome_copy_number",sep='/')
parameter$output_DeAR = paste(output_dir,"adjust_distance",sep='/')
parameter$output_QC = paste(output_dir,"QC",sep='/')
parameter$preProcess = paste(output_dir,preProcess,sep="/")
parameter$Callcnv_dir = paste(output_dir,Callcnv_dir,sep="/")
parameter$result_dir = paste(output_dir,result_dir,sep="/")
parameter$result_graph = paste(parameter$result_dir,'graph',sep="/")
if(parameter$IT_version){
	parameter$output_sex = paste(output_dir,"uniq_mapped.sex",sep='/')
}else{
	parameter$output_sex = paste(output_dir,"sex",sep='/')
}


parameter$input_dir = input_dir
parameter$output_dir = output_dir
parameter$program_dir = program_dir
parameter$Ncores = 4
parameter$Totalbin = 1250
parameter$Nbin = parameter$Totalbin / parameter$nStep
parameter$window_len = 200000
parameter$GD_score_cutoff = 5       # GD cutoff
parameter$GD_length_cutoff = 100    # GD segment length cutoff, if segment length < this cutoff, then stop GD
parameter$SW_cutoff = 5             # SW cutoff
parameter$SW_score_penalty = 3      # SW penalty, for my suggestion, set 1,2 or 3
parameter$SW_break_num = 5          # SW break CNV number, suggest 3~10
parameter$merge_cutoff = 0.0        # merge cutoff, if copy number difference in two adjacent cnv, then merge it.
parameter$ra_off = 0.6
parameter$cn_off = 0.01
parameter$no_na_ratio = 0.2         # used in reference generation
parameter$DeAR = FALSE              # whether do DeAR correction
parameter$Na_penalty = 0.01
parameter$bound = 0.3

parameter$file_20K = input_dir

parameter$chrs = paste('chr', c(1:24), sep='')

parameter$Length_chr = c(
  249250621, 243199373, 198022430, 191154276, 180915260,
  171115067, 159138663, 146364022, 141213431, 135534747,
  135006516, 133851895, 115169878, 107349540, 102531392,
  90354753, 81195210, 78077248, 59128983, 63025520,
  48129895, 51304566, 155270560, 59373566)
parameter$centromereS = c(
  121535434,92326171,90504854,49660117,46405641,
  58830166,58054331,43838887,47367679,39254935,
  51644205,34856694,16000000,16000000,17000000,
  35335801,22263006,15460898,24681782,26369569,
  11288129,13000000,58632012,10104553)
parameter$centromereE = c(
  124535434,95326171,93504854,52660117,49405641,
  61830166,61054331,46838887,50367679,42254935,
  54644205,37856694,19000000,19000000,20000000,
  38335801,25263006,18460898,27681782,29369569,
  14288129,16000000,61632012,13104553)

parameter$DeAR_direction = c(0.056007117093369,-0.059805936008050,-0.088717419624836,-0.175538875621678,-0.095587871828156,
-0.088444219480323,-0.040526153217012,-0.057790761296531,0.031051592895438,0.041412645685552,0.036496815962586,0.002061002747194,
-0.148555952650169,0.010592009355499,0.088636283980166,0.256592289322018,0.338112146434760,-0.085139495020531,0.657450223485251,
0.197581344913872,0.041171365398495,0.488082842232818)
# the direction vector in DeAR correction (before GC)

parameter$DeAR_direction = c(0.031251326562927,-0.086005650384720,-0.088563219949689,-0.100296873582170,-0.087498452101635,
-0.072330536841494,-0.036403402903086,-0.038124636542417,0.014691400549564,0.006504200841056,-0.023355474054425,0.040938841323958,
-0.084371673051938,0.049967216270285,0.071372859189773,0.256910405288226,0.333947842492641,-0.071994789411517,0.789508274540154,
0.097878602174229,0.091721782573936,0.341824969208459)
# the direction vector in DeAR correction (after GC)

parameter$DeAR_triploid = c(2.879451752386346,2.869896128892775,2.891286193481981,2.897144140163833,2.902073870355654,
2.907359884635960,2.919276175203279,2.919752074118879,2.938649807959235,2.928813718172947,2.927217660099291,2.928509232613440,
2.945497978161096,2.950894170115114,2.956821535122604,2.959014687562571,2.960073054602987,2.956704969136927,2.974533531172864,
2.965437822182463,2.980697006187927,2.981824290040060)
# the observe triploid value
parameter$DeAR_haploid = c(1.043692081612546,1.047484570552276,1.039067374687176,1.036807173046746,1.034919744561109,
1.032910597456148,1.028436262662622,1.028259131583750,1.021320056172603,1.024908872789848,1.025495878532064,1.025020756670322,
1.018850258049299,1.016920540879697,1.014817348099602,1.014043488269038,1.013670872218170,1.014858543941853,1.008633397566470,
1.011790393454930,1.006516185786784,1.006130858385573)
# the observe haploid value


library(parallel)

### check files

message("[check data ...]")



if(!file.exists(parameter$output_dir)) dir.create(parameter$output_dir)
if(file.exists(parameter$output_sex)) unlink(parameter$output_sex)

if(!parameter$IT_version){

  if(!file.exists(parameter$preProcess)) dir.create(parameter$preProcess)
  if(!file.exists(parameter$Callcnv_dir)) dir.create(parameter$Callcnv_dir)
  if(!file.exists(parameter$result_dir)) dir.create(parameter$result_dir)
  if(!file.exists(parameter$result_graph)) dir.create(parameter$result_graph)
  if(file.exists(parameter$output_copy)) unlink(parameter$output_copy)
  if(file.exists(parameter$output_DeAR)) unlink(parameter$output_DeAR)
  if(file.exists(parameter$output_QC)) unlink(parameter$output_QC)
  
}


detectCNV_core <- function(alldata,parameter,ideoCyto){
  # core function in CNV detection
  
  if(alldata$sumCount == 'None'){
    return(list(sd = 'None', sampleName = alldata$sampleName, sumCount = 'None'))
  }
  
  .libPaths("/share/work1/wangrr/local/Rlib/")
  library(zoo)
  library(Rcpp)
  library(rmarkdown)
  
  nBlock = parameter$nBlock
  nStep = parameter$nStep
  Nbin = parameter$Nbin
  window_len = parameter$window_len
  Totalbin = parameter$Totalbin
  GD_score_cutoff = parameter$GD_score_cutoff
  GD_length_cutoff = parameter$GD_length_cutoff
  SW_cutoff = parameter$SW_cutoff
  SW_score_penalty = parameter$SW_score_penalty
  SW_break_num = parameter$SW_break_num
  merge_cutoff = parameter$merge_cutoff
  ra_off = parameter$ra_off
  chrs = parameter$chrs
  cn_off = parameter$cn_off
  bound = parameter$bound
  Na_penalty = parameter$Na_penalty
  Length_chr = parameter$Length_chr
  centromereS = parameter$centromereS
  centromereE = parameter$centromereE
  input_dir = parameter$input_dir
  output_dir = parameter$output_dir
  program_dir = parameter$program_dir
  
  cytoband_file = parameter$cytoband_file
  html_plot = parameter$html_plot
  output_sex = parameter$output_sex
  preProcess = parameter$preProcess
  Callcnv_dir = parameter$Callcnv_dir
  result_dir = parameter$result_dir
  result_graph = parameter$result_graph
  CV = parameter$CV
  reference = parameter$reference
  reference_file = parameter$reference_file
  file_20K = parameter$file_20K
  file_GC = parameter$file_GC
  
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
  
  gradient_descent_recursive <- function(fd,rdata,ref_data,score,xcoor,distance,distance_ref,yh){
    # GD main flow
    
    distance[distance > 2] = 2
    
    fd_new = c()
    index_new = c()
    num = 1
    not_na = c()
    
    zhong = median(fd)
    for(i in 1:length(fd)){
      fd_new = c(fd_new,fd[i])
      index_new = c(index_new,i)
      not_na = c(not_na,1)
      num = num + 1
      if(i == length(fd)) break
      if(distance[i] == 1) next
      for(j in 1:(distance[i]-1)){
        fd_new = c(fd_new,zhong)
        not_na = c(not_na,NA)
        index_new = c(index_new,NA)
        num = num + 1
      }
    }
    
    index_new_s = index_new
    for(i in length(index_new):1){
      if(is.na(index_new[i])){
        index_new_s[i] = index_new_s[i+1]
      }
    }
    
    index_new_e = index_new
    for(i in 1:length(index_new)){
      if(is.na(index_new[i])){
        index_new_e[i] = index_new_e[i-1]
      }
    }
    
    
    dimension = length(fd_new)
    #dimension = length(fd)
    breakpoint = breakpoint_old = c()
    
    index_s = 1
    index_e = dimension
    
    iter_times = 0
    seg_active = 1      # if seg_active == 0, then do not apply GD in this segment
    
    
    while(TRUE){
      
      iter_times = iter_times + 1;
      seg_gradient = seg_error = rep(0,length(index_s))
      breakpoint_old = breakpoint
      
      for(i in 1:(length(breakpoint_old)+1)){
        
        if(seg_active[i] == 0){next}
        
        #seg_info = fast_gradient(fd[index_s[i]:index_e[i]],index_s[i],index_e[i])
        seg_info = .Call("fast_gradient_so",fd_new[index_s[i]:index_e[i]],index_s[i],index_e[i]);
        #seg_info = .Call("fast_gradient_so",fd[index_s[i]:index_e[i]],index_s[i],index_e[i]);
        chang = index_e[i] - index_s[i] + 1
        seg_gradient[i] = seg_info$gradient
        
        if(seg_gradient[i] < max(GD_score_cutoff * log(chang) / log(10000))){
          seg_active[i] = 0
        }else{
          position = seg_info$breakpoint
          breakpoint = c(breakpoint,position)
          index_s = c(index_s,position+1)
          index_e = c(index_e,position)
          seg_active = c(seg_active,1)
          temp_fd = fd_new * not_na
          for(j in 1:length(index_s)){
            temp_fd2 = temp_fd[index_s[j]:index_e[j]]
            zhong = median(temp_fd2,na.rm=T)
            temp_fd2[is.na(temp_fd2)] = zhong
            fd_new[index_s[j]:index_e[j]] = temp_fd2
          }
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
    
    index_s = index_new_s[index_s]
    index_e = index_new_e[index_e]
    
    segments = cal_seginfo(fd,rdata,ref_data,index_s,index_e,score,xcoor,distance_ref,yh) 
    
    return(segments)
  }
  
  gradient_descent_recursive2 <- function(fd,rdata,ref_data,score,xcoor,distance,distance_ref,yh){
    # GD main flow
    
    #fd_new = c()
    #index_new = c()
    #num = 1
    #not_na = c()
    
    #zhong = median(fd)
    #for(i in 1:length(fd)){
    #  fd_new = c(fd_new,fd[i])
    #  index_new = c(index_new,i)
    #  not_na = c(not_na,1)
    #  num = num + 1
    #  if(i == length(fd)) break
    #  if(distance[i] == 1) next
    #  for(j in 1:(distance[i]-1)){
    #    fd_new = c(fd_new,zhong)
    #    not_na = c(not_na,NA)
    #    index_new = c(index_new,NA)
    #    num = num + 1
    #  }
    #}
    
    #for(i in length(index_new):1){
    #  if(is.na(index_new[i])){
    #    index_new[i] = index_new[i+1]
    #  }
    #}
    
    #dimension = length(fd_new)
    dimension = length(fd)
    breakpoint = breakpoint_old = c()
    
    index_s = 1
    index_e = dimension
    
    iter_times = 0
    seg_active = 1      # if seg_active == 0, then do not apply GD in this segment
    
    
    while(TRUE){
      
      iter_times = iter_times + 1;
      seg_gradient = seg_error = rep(0,length(index_s))
      breakpoint_old = breakpoint
      
      for(i in 1:(length(breakpoint_old)+1)){
        
        if(seg_active[i] == 0){next}
        
        #seg_info = fast_gradient(fd[index_s[i]:index_e[i]],index_s[i],index_e[i])
        #seg_info = .Call("fast_gradient_so",fd_new[index_s[i]:index_e[i]],index_s[i],index_e[i]);
        seg_info = .Call("fast_gradient_so",fd[index_s[i]:index_e[i]],index_s[i],index_e[i]);
        chang = index_e[i] - index_s[i] + 1
        seg_gradient[i] = seg_info$gradient
        
        if(seg_gradient[i] < max(GD_score_cutoff * log(chang) / log(10000))){
          seg_active[i] = 0
        }else{
          position = seg_info$breakpoint
          breakpoint = c(breakpoint,position)
          index_s = c(index_s,position+1)
          index_e = c(index_e,position)
          seg_active = c(seg_active,1)
          #temp_fd = fd_new * not_na
          #for(j in 1:length(index_s)){
          #  temp_fd2 = temp_fd[index_s[j]:index_e[j]]
          #  zhong = median(temp_fd2,na.rm=T)
          #  temp_fd2[is.na(temp_fd2)] = zhong
          #  fd_new[index_s[j]:index_e[j]] = temp_fd2
          #}
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
    
    #index_s = index_new[index_s]
    #index_e = index_new[index_e]
    
    segments = cal_seginfo(fd,rdata,ref_data,index_s,index_e,score,xcoor,distance_ref,yh) 
    
    return(segments)
  }
  
  find_contig <- function(temp_fd,temp_rd,temp_ref,temp_score,temp_distance,s,e,yh){
    # find segment in SW
    
    if(s == e){return(list(index_s = s,index_e = e))
    }else{
      
      ave_fd = mean(temp_fd)
      dimension = length(temp_fd)
      
      
      if(ave_fd < homo_ratio) return(list(index_s = s,index_e = e))
      
      ratio3 = (ave_fd+bound) / ave_fd
      ratio1 = (ave_fd-bound) / ave_fd
      
      if(ratio1<homo_ratio) ratio1 = homo_ratio
      
      lambda = temp_ref * (ave_fd / yh)
      
      lambda3 = lambda * ratio3
      lambda1 = lambda * ratio1
      #lambda1[lambda1 < 1] = 1
      
      prob3 = sapply(1:dimension,function(x) dnorm(temp_rd[x],mean=lambda3[x],sd=sqrt(lambda3[x])))
      prob2 = sapply(1:dimension,function(x) dnorm(temp_rd[x],mean=lambda[x],sd=sqrt(lambda[x])))
      prob1 = sapply(1:dimension,function(x) dnorm(temp_rd[x],mean=lambda1[x],sd=sqrt(lambda1[x])))
      
      prob3[prob3 < 1e-300] = 1e-300
      prob2[prob2 < 1e-300] = 1e-300
      prob1[prob1 < 1e-300] = 1e-300
      
      
      wucha3 = (prob3 / (prob2 + prob3) - 0.5) * 2
      wucha1 = (prob1 / (prob2 + prob1) - 0.5) * 2
      
      breakpoint1 = find_position(wucha1,temp_distance)
      breakpoint3 = find_position(wucha3,temp_distance)
      breakpoint = union(breakpoint1,breakpoint3)
      
      if(is.null(breakpoint)){
        return(list(index_s = s,index_e = e))
      }else{
        breakpoint = sort(breakpoint)
        return(list(index_s = c(s,breakpoint + s),index_e = c(breakpoint + s - 1,e)))
      }
    }
  }
  
  find_position <- function(wucha,temp_distance){
    # calculate region score in SW
    
    dimension = length(wucha)
  
    breakpoint = c()
    
    i = 1
    while(i <= dimension-1){
      match_array = rep(0,dimension)
      match_array[i] = wucha[i]
      for(j in (i+1):length(wucha)){
        match_array[j] = ifelse(match_array[j-1] < 0,break,
                                match_array[j-1] + wucha[j] - (temp_distance[j-1]-1) * Na_penalty)
        if(j > SW_break_num){
          deter = sapply(1:(SW_break_num),function(x) match_array[j-x+1]<match_array[j-x]-(temp_distance[j-x]-1)*Na_penalty)
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
  
  match_sequence <- function(gradient_list,fd,rdata,ref_data,score,xcoor,distance,distance_ref,yh){
    # main process in SW
    
    index_s_new = c()
    index_e_new = c()
    
    index_s = gradient_list$index_s
    index_e = gradient_list$index_e
    
    for(i in 1:length(index_s)){
      temp_fd = fd[index_s[i]:index_e[i]]
      temp_score = score[index_s[i]:index_e[i]]
      temp_rd = rdata[index_s[i]:index_e[i]]
      temp_ref = ref_data[index_s[i]:index_e[i]]
      temp_distance = distance[index_s[i]:(index_e[i]-1)]
      index_temp = find_contig(temp_fd,temp_rd,temp_ref,temp_score,temp_distance,index_s[i],index_e[i],yh)
      index_s_new = c(index_s_new, index_temp$index_s)
      index_e_new = c(index_e_new, index_temp$index_e)
    }
    
    segments = cal_seginfo(fd,rdata,ref_data,index_s_new,index_e_new,score,xcoor,distance_ref,yh)
    
    return(segments)
  }
  
  cal_seginfo <- function(fd,rdata,ref_data,index_s,index_e,score,xcoor,distance_ref,yh){
    # calculate segment information
    # include start index (index_s), end index (index_e), copy number (copy), segment lenght(segLen) 
    
    #copy = sapply(1:length(index_s),function(x) sum(fd[index_s[x]:index_e[x]]*score[index_s[x]:index_e[x]])/sum(score[index_s[x]:index_e[x]]))
    #copy = sapply(1:length(index_s),function(x) median(fd[index_s[x]:index_e[x]]))
    copy = rep(0,length(index_s))
    score_region = rep(0,length(index_s))
    obs = rep(0,length(index_s))
    ref = rep(0,length(index_s))
    for(i in 1:length(index_s)){
      obs[i] = sum(rdata[index_s[i]:index_e[i]])
      if(index_s[i] == index_e[i]){
        ref[i] = ref_data[index_s[i]]
      }else{
        ref[i] = sum(ref_data[index_s[i]:index_e[i]]) + sum(distance_ref[index_s[i]:(index_e[i]-1)],na.rm=T)
      }
      copy[i] = obs[i] / ref[i] * yh
      score_region[i] = -log(dnorm(obs[i],mean=ref[i],sd=sqrt(ref[i])),10)
    }
    segLen = sapply(1:length(index_s),function(x) (index_e[x] - index_s[x] + 1))
    
    change = c()
    if(length(copy) == 1){change = NULL
    }else{change = sapply(1:(length(index_s)-1),function(x) abs(copy[x+1]-copy[x]))}
    
    return(list(index_s = index_s, index_e = index_e, copy = copy, segLen = segLen, score = score_region,obs=obs,ref=ref,
                change = change, xcoor_s = xcoor[index_s]-xcoor.offset, xcoor_e = xcoor[index_e]+xcoor.offset))
  }
  
  combine_segments <- function(segments,fd,rdata,ref_data,score,xcoor,distance_ref,yh){
    # merge two segment when copy number different < merge cutoff
    
    index_s = segments$index_s
    index_e = segments$index_e
    copy = segments$copy
    change = segments$change
    score_region = segments$score
    
    if(length(copy) == 1){return(segments)}
    
    com_point = c()
    change_max = 999
    
    for(i in 1:length(change)){
      if(change[i] > merge_cutoff){
        next
      }else{
        if(abs(copy[i]-yh)>cn_off & abs(copy[i+1]-yh)>cn_off & change[i] < change_max){
          com_point = i
          change_max = change[i]
        }
      }
    }
    
    if(length(com_point) == 0){return(segments)}
    
    seg_info = cal_seginfo(fd,rdata,ref_data,index_s[-(com_point+1)],index_e[-com_point],score,xcoor,distance_ref,yh)
    
    return(combine_segments(seg_info,fd,rdata,ref_data,score,xcoor,distance_ref,yh))
  }
  
  mhtplot <- function(out, da, name, hasxy) {
    
    x.coor = da[1,]
    name = strsplit(name, '_', perl=T)[[1]][1]
    if (hasxy == "noXY") {
      chrNum = 22
    }else if (hasxy == "XY") {
      chrNum = 24
    }else{
      cat("the type is not XY or noXY", stderr())
      quit(status=1)
    }
    x.all = c()
    y.all = c()
    breaks.pos = cumsum(c(0, Length_chr[1:(chrNum-1)]+20000000)) + Length_chr[1:chrNum]/2
    chrs.pos = cumsum(Length_chr[1:chrNum] + 20000000) - 10000000
    for (i in 1:chrNum){
      x.index = which(x.coor < Length_chr[i])
      x.cur = x.coor[x.index]
      y.cur = da[2+(i-1)*2, ][x.index]
      lstMrk = max(x.cur, na.rm=T)
      x.cur = rollapply(x.cur, width=5, by=5, FUN=mean, na.rm=T, partial=T, align='left')
      y.cur = rollapply(y.cur, width=5, by=5, FUN=mean, na.rm=T, partial=T, align='left')
      if (i == 1) {
        x.all = x.cur
      }else{
        x.all = c(x.all, x.cur + max(x.all, na.rm=T) + 20000000)
      }
      y.all = c(y.all, y.cur)
    }
    y.all[!is.na(y.all)&y.all>4] = 4
    png(out, width=2000, height=400, units='px', type='cairo', pointsize=12)
    par(mai=c(0.82, 1.02, 0.82,0.42))
    x.max = max(x.all, na.rm=T)
    Name = 'Whole Genome Detection Results'
    Encoding(Name) = 'UTF-8'
    Ylab = 'Copy Number'
    Encoding(Ylab) = 'UTF-8'
    plot(x.all, y.all, type='n', main=Name, cex.main=2, cex.lab=2, ylim=c(0, 4), axe=F, 
         ylab=Ylab, xlab='', xlim=c(100000000, x.max-100000000), cex.lab=1.75)
    chr_all_name = c(1:22, 'X', 'Y')
    axis(1, at=breaks.pos, labels=chr_all_name[1:chrNum], cex.axis=1.5)
    lbl = c(0:4)
    axis(2, at=lbl, labels=lbl, cex.axis=1.5)
    abline(h=lbl, lty=2, col='gray')
    points(x.all[!is.na(y.all)&y.all>=2], y.all[!is.na(y.all)&y.all>=2], pch=20, col='blue', cex=0.5)
    points(x.all[!is.na(y.all)&y.all<2], y.all[!is.na(y.all)&y.all<2], pch=20, col='red', cex=0.5)
    
    x.moving.average = rollapply(x.all, width=100, by=50, FUN=mean, na.rm=T, partial=T, align="left")
    y.moving.average = rollapply(y.all, width=100, by=50, FUN=mean, na.rm=T, partial=T, align="left")
    
    lines(x.moving.average, y.moving.average, lwd=2, col=1)
    abline(v=c(0,chrs.pos), lwd=2, col='gray')
    box()
    dev.off()
  }
  
  scatterplot <- function(outdir, da, name, sex, ideoCyto) {
    
    x.coor = da[1, ]
    chrNum = 24
    for (i in 1:chrNum) {
      x.index = which(x.coor<Length_chr[i])
      x.cur = x.coor[x.index]
      y.cur = da[2+(i-1)*2, ][x.index]
      l.cur = da[3+(i-1)*2, ][x.index]
      Nchr = paste('chr', i, sep='')
      if (i == 23) { Nchr = 'chrX'}
      if (i == 24) { Nchr = 'chrY'}
      out.pdf = paste(outdir, '/', name, '.', Nchr, '.png', sep='')
      png(out.pdf, width=2000, height=400, units='px', type='cairo', pointsize=12)
      x.lab.offset = max(x.cur)/2000
      par(mai=c(0.82, 1.02, 0.82,0.42))
      if (i==23) {
        Name = 'Chromosome X'
      }else if (i==24) {
        Name = 'Chromosome Y'
      }else{
        Name = paste('Chromosome ',i, sep='')
      }
      Encoding(Name) = 'UTF-8'
      Ylab = 'Copy Number'
      Encoding(Ylab) = 'UTF-8'
      Xlab = 'Position (Mb)'
      Encoding(Xlab) = 'UTF-8'
      plot(x.cur, y.cur, ylim=c(-1,4), type='n', main=Name, cex.main=2, cex.lab=2, xlab=Xlab, 
           ylab=Ylab, axe=F, xlim=c(50*x.lab.offset, max(x.cur)-50*x.lab.offset), cex.lab=1.75)
      if (max(x.cur) >= 100000000) {
        x.at = seq(from=0, by=20000000, to=max(x.cur)+20000000)
        x.labels = seq(from=0, by=20, to=(max(x.cur)+20000000)/1000000)
      }else{
        x.at = seq(from=0, by=10000000, to=max(x.cur)+20000000)
        x.labels = seq(from=0, by=10, to=(max(x.cur)+20000000)/1000000)
      }
      
      lbl = c(0:4)
      axis(2, at=lbl, labels=lbl, cex.axis=1.5)
      axis(1, at=x.at, labels=x.labels, cex.axis=1.5)
      abline(h=lbl, lwd=2, lty=2, col='gray')
      points(x.cur, y.cur, pch=20, cex=0.75, col=rgb(0.5,0.5,0.5,1))
      lines(x.cur[!is.na(l.cur)&x.cur<centromereS[i]], l.cur[!is.na(l.cur)&x.cur<centromereS[i]],
            lwd=2, col='blue')
      lines(x.cur[!is.na(l.cur)&x.cur>centromereE[i]], l.cur[!is.na(l.cur)&x.cur>centromereE[i]],
            lwd=2, col='blue')
      
      shifting = 0.02
      yh = 0
      if ((i==23&sex=='M')|i==24) { yh = -1 }
      for (j in x.cur[is.na(y.cur)]) {
        segments(j-shifting, 2+yh, j+shifting, 2+yh, lwd=4, col='red')
      }
      segments(centromereS[i], 2+yh, centromereE[i], 2+yh, lwd=8)
      #### add cytoband plot
      colnames(ideoCyto) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
      chrCyto <- subset(ideoCyto, chrom==Nchr)
      gieStainCol <- data.frame(type=c("stalk", "gvar", "acen", "gneg", "gpos25", "gpos50", "gpos75", "gpos100"),
                                Col=c("gray50", "gray50", "gray50", "gray100", "gray75", "gray50", "gray25", "gray0"),
                                stringsAsFactors=F)
      chrCyto$Col <- unlist(lapply(chrCyto$gieStain, FUN=function(x) {gieStainCol$Col[gieStainCol$type == x]}))
      
      for(i in seq(1, nrow(chrCyto))) {
        tmpStart <- chrCyto[i,2]
        tmpEnd <- chrCyto[i,3]
        tmpCol <- chrCyto[i,6]
        if (chrCyto[i,5]=="acen" & substring(chrCyto[i,4], 1, 1)=="p") {
          polygon( x=c(tmpStart, tmpEnd, tmpStart), y=c(-0.3, -0.6 , -0.9), col="orangered4" )
        }else if (chrCyto[i,5]=="acen" & substring(chrCyto[i,4], 1, 1)=="q") {
          polygon( x=c(tmpEnd, tmpStart, tmpEnd), y=c(-0.3, -0.6 , -0.9), col="orangered4" )
        }else{
          rect(chrCyto[i,2], -0.3, chrCyto[i,3], -0.9, col=tmpCol, lwd=1)
        }
        if (chrCyto[i,5]!='acen' & tmpEnd-tmpStart > 60*max(x.cur)/2000) {
          if (chrCyto[i,5] %in% c("stalk", "gvar", "acen", "gpos50", "gpos75", "gpos100")) {
            text(x=(tmpEnd+tmpStart)/2, y=-0.6, labels=chrCyto[i, 4], cex=1.5, col='white')
          }else{
            text(x=(tmpEnd+tmpStart)/2, y=-0.6, labels=chrCyto[i, 4], cex=1.5, col='black')
          }
        }
      }
      
      box()
      dev.off() 
    }
  }
  
  homo_supplement <- function(Param,data){
    
    supplement = c()
    fd = data/Param
    
    for(i in 1:24){
      deter = rep(0,Nbin)
      sum_deter = rep(0,Nbin)
      deter[!is.na(Param[i,]) & is.na(fd[i,])] = 1
      sum_deter[1] = ifelse(deter[1] == 0,0,1)
      num_notNA = 0;
      for(j in 2:Nbin){
        if(deter[j] == 1){
          sum_deter[j] = sum_deter[j-1] + 1
          num_notNA = 0
        }else{
          if(is.na(fd[i,j])){
            temp = sum_deter[j-1] - 0.5
          }else{
            temp = sum_deter[j-1] + 1 - fd[i,j]
          }
          if(num_notNA >= 5){
            sum_deter[j] = 0
            num_notNA = 0
            next
          }
          if(temp > 0){
            sum_deter[j] = temp
            num_notNA = num_notNA + 1
          }else{
            sum_deter[j] = 0
            num_notNA = 0
          }
        }
      }
      index = c()
      max_value = 0
      for(j in 1:Nbin){
        if(sum_deter[j] != 0){
          index = c(index,j)
          if(sum_deter[j] > max_value){
            max_value = sum_deter[j]
          }
          next
        }
        if(max_value < 5){
          index = c();
          max_value = 0;
          next;
        }
        start = index[1]
        for(k in index[1]:index[length(index)]){
          if(deter[k] == 1){
            data[i,k] = 0.01
          }
          if(sum_deter[k] >= max_value){
            max_value = 0;
            index = c();
            stop = k
            supplement = c(supplement,paste(i,start,stop))
            break
          }
        }
      }
    }
    
    return(list(data = data, supplement = supplement))
  }
  
  ContinueRegion_GD <- function(gradient_list,yh,dimension){
    
    index_s = gradient_list$index_s
    index_e = gradient_list$index_e
    copy = gradient_list$copy
    
    EstiCopy = 1:dimension
    for(j in 1:length(copy)){EstiCopy[index_s[j]:index_e[j]] = copy[j]}
    EstiCopy[abs(EstiCopy - yh) < cn_off] = yh
    
    return(EstiCopy)
  }
  
  Call_cnv_GD <- function(xcoor, ldata,rdata,ref_data, i, sex){
    
    # main workflow for GD algorithm
    
    #score = ref_data / mean(ref_data,na.rm=T)
    #score[score > 1] = 1
    score = rep(1,length(ref_data))
    
    index = which(xcoor <= Length_chr[i])
    if (i == 24) {ldata[xcoor > 25000000] = NA}
    
    index.p = which((xcoor<=centromereS[i])&!is.na(ldata))
    index.q = which(xcoor>=centromereE[i]&xcoor<=Length_chr[i]&!is.na(ldata))
    xcoor.p = xcoor[index.p]
    xcoor.q = xcoor[index.q]
    if(length(index.p)>2) distance.p = index.p[2:length(index.p)] - index.p[1:(length(index.p)-1)]
    if(length(index.q)>2) distance.q = index.q[2:length(index.q)] - index.q[1:(length(index.q)-1)]
    if(length(index.p)>2) distance_ref.p = round(sapply(1:(length(index.p)-1),function(x) sum(ref_data[(index.p[x]):(index.p[x+1])],na.rm=T)-ref_data[(index.p[x])]-ref_data[(index.p[x+1])]),3)
    if(length(index.q)>2) distance_ref.q = round(sapply(1:(length(index.q)-1),function(x) sum(ref_data[(index.q[x]):(index.q[x+1])],na.rm=T)-ref_data[(index.q[x])]-ref_data[(index.q[x+1])]),3)
    
    ldata.p = ldata[index.p]
    ldata.q = ldata[index.q]
    rdata.p = rdata[index.p]
    rdata.q = rdata[index.q]
    ref_data.p = ref_data[index.p]
    ref_data.q = ref_data[index.q]
    score.p = score[index.p]
    score.q = score[index.q]
    
    chrN = paste('chr', i, sep='')
    if (i == 23) { chrN = 'chrX' }
    if (i == 24) { chrN = 'chrY' }
    
    ## sex type param
    yh = ifelse((i==23&sex=='M')|(i==24&sex=='M'), 1, 2)
    res.fin = rep(NA, length(xcoor))
    if (i==24&sex=='F'){
      return(res.fin)
    }
    res.gd = rep(NA, length(xcoor))
    res.sw = rep(NA, length(xcoor))

    if(length(ldata[!is.na(ldata)]) < 15){
      CN = rep(NA,length(xcoor))
      CN[!is.na(ldata)] = 0
      if(!parameter$IT_version){
        #scatterplot_pdf(CN,ldata,xcoor,xcoor.p,xcoor.q,i,index,index.p,index.q,yh,'GD')
        #scatterplot_pdf(CN,ldata,xcoor,xcoor.p,xcoor.q,i,index,index.p,index.q,yh,'SW')
        #scatterplot_pdf(CN,ldata,xcoor,xcoor.p,xcoor.q,i,index,index.p,index.q,yh,'Merge')
      }
      return(CN)
    }
    
    if (i %in% c(13,14,15,21,22)) {
      
      gradient_list = gradient_descent_recursive(ldata.q,rdata.q,ref_data.q,score.q,xcoor.q,distance.q,distance_ref.q,yh)    # GD
      sequence_list = match_sequence(gradient_list,ldata.q,rdata.q,ref_data.q,score.q,xcoor.q,distance.q,distance_ref.q,yh)   # sw
      gradient_merge = combine_segments(sequence_list,ldata.q,rdata.q,ref_data.q,score.q,xcoor.q,distance_ref.q,yh)    # Merge result
      
      res.gd[index.q] = ContinueRegion_GD(gradient_list,yh,length(ldata.q))
      res.sw[index.q] = ContinueRegion_GD(sequence_list,yh,length(ldata.q))
      res.fin[index.q] = ContinueRegion_GD(gradient_merge,yh,length(ldata.q))
      
      combine_result = list(xcoor_s = gradient_merge$xcoor_s,xcoor_e = gradient_merge$xcoor_e,
                            copy = (gradient_merge$copy), score = (gradient_merge$score),
                            obs = (gradient_merge$obs),ref = (gradient_merge$ref))
      
    }else{
      
      gradient_list_p =  gradient_descent_recursive(ldata.p,rdata.p,ref_data.p,score.p,xcoor.p,distance.p,distance_ref.p,yh)
      sequence_list_p = match_sequence(gradient_list_p,ldata.p,rdata.p,ref_data.p,score.p,xcoor.p,distance.p,distance_ref.p,yh)
      gradient_merge_p = combine_segments(sequence_list_p,ldata.p,rdata.p,ref_data.p,score.p,xcoor.p,distance_ref.p,yh)
      
      res.gd[index.p] = ContinueRegion_GD(gradient_list_p,yh,length(ldata.p))
      res.sw[index.p] = ContinueRegion_GD(sequence_list_p,yh,length(ldata.p))
      res.fin[index.p] = ContinueRegion_GD(gradient_merge_p,yh,length(ldata.p))
      
      gradient_list_q = gradient_descent_recursive(ldata.q,rdata.q,ref_data.q,score.q,xcoor.q,distance.q,distance_ref.q,yh)
      sequence_list_q = match_sequence(gradient_list_q,ldata.q,rdata.q,ref_data.q,score.q,xcoor.q,distance.q,distance_ref.q,yh)
      gradient_merge_q = combine_segments(sequence_list_q,ldata.q,rdata.q,ref_data.q,score.q,xcoor.q,distance_ref.q,yh)
      
      res.gd[index.q] = ContinueRegion_GD(gradient_list_q,yh,length(ldata.q))
      res.sw[index.q] = ContinueRegion_GD(sequence_list_q,yh,length(ldata.q))
      res.fin[index.q] = ContinueRegion_GD(gradient_merge_q,yh,length(ldata.q))
      
      combine_result = list(xcoor_s = c(gradient_merge_p$xcoor_s,gradient_merge_q$xcoor_s),
                            xcoor_e = c(gradient_merge_p$xcoor_e,gradient_merge_q$xcoor_e),
                            copy = (c(gradient_merge_p$copy,gradient_merge_q$copy)),
                            score = (c(gradient_merge_p$score,gradient_merge_q$score)),
                            obs = (c(gradient_merge_p$obs,gradient_merge_q$obs)),
                            ref = (c(gradient_merge_p$ref,gradient_merge_q$ref)))
      
    }
      ##### output cnv region and filter some region

      CNV_region(combine_result,yh,ldata,i)

      if(!parameter$IT_version){
        #scatterplot_pdf(res.gd,ldata,xcoor,xcoor.p,xcoor.q,i,index,index.p,index.q,yh,'GD')
        #scatterplot_pdf(res.sw,ldata,xcoor,xcoor.p,xcoor.q,i,index,index.p,index.q,yh,'SW')
        #scatterplot_pdf(res.fin,ldata,xcoor,xcoor.p,xcoor.q,i,index,index.p,index.q,yh,'Merge')
      }
      
      return(res.fin)
  }
  
  scatterplot_pdf <- function(res,ldata,xcoor,xcoor.p,xcoor.q,i,index,index.p,index.q,yh,name_this){
    
    chrN = paste('chr', i, sep='')
    if (i == 23) { chrN = 'chrX' }
    if (i == 24) { chrN = 'chrY' }
    
    shifting = 0.01*nStep
    lbl = c(0:4)
    plot(xcoor[index]/1000000, ldata[index], pch=20, type='p',col=rgb(0.5,0.5,0.5,1),xlab='chromosome pos(Mb)', 
         ylab='copy number', main=paste(chrN,': ',name_this,sep=''), cex=0.25, ylim=c(0,4))
    abline(h=lbl, lty=c(2,2,1,2,2), col=c('gray', 'gray', 'black', 'gray', 'gray'))
    lines(xcoor.q/1000000, res[index.q], lwd=2, col='blue')
    if (i %in% c(13,14,15,21,22)) {
    }else{
      lines(xcoor.p/1000000, res[index.p], lwd=2, col='blue')
    }
    for (j in xcoor[index][is.na(ldata[index])]/1000000) {
      segments(j-shifting, yh, j+shifting, yh, lwd=4, col='red')
    }
    segments(centromereS[i]/1000000, yh, centromereE[i]/1000000, yh, lwd=4)
  }
  
  CNV_region <- function(combine_result,yh,ldata,i){
    # generate CNV region file
    
    chrN = paste('chr', i, sep='')
    if (i == 23) { chrN = 'chrX' }
    if (i == 24) { chrN = 'chrY' }
    
    stp_x = combine_result$xcoor_s
    edp_x = combine_result$xcoor_e
    copy = combine_result$copy
    score = combine_result$score
    obs = combine_result$obs
    ref = combine_result$ref
    
    for (k in 1:length(stp_x)) {
      
      stp_y = min(which(x.coor>stp_x[k]&x.coor<edp_x[k]))
      edp_y = max(which(x.coor>stp_x[k]&x.coor<edp_x[k]))
      
      stp_y = stp_x[k] / window_len / nStep + 1
      edp_y = (edp_x[k]-window_len*nBlock) / window_len / nStep + 1
      
      Length_region = edp_x[k] - stp_x[k]
      L_y = edp_y - stp_y + 1
      
      lf.l = ldata[seq(stp_y, edp_y)]
      lf.na = sum(is.na(lf.l))

      if (i == 24) {
        ra_off = 1
      }
      
      if (copy[k]<yh+cn_off&copy[k]>yh-cn_off) { 
        next
      } 
  
      #if (lf.na/L_y < ra_off) {
      if(TRUE){
        zvalid = ifelse((copy[k]>2)|(copy[k]<2), 'zValid', 'z=0')
        edp_x[k] = ifelse(edp_x[k]>Length_chr[i], Length_chr[i], edp_x[k])
        cat(file=out_cnv_result, append=T, paste(paste('>','GDSW',sep=''), '\t', chrN, ':', stp_x[k], '-',edp_x[k], '\t', 
                      round(copy[k],3), '\t', L_y,'[', stp_y, '-', edp_y, '(', lf.na, ' NA', ']\t', 
                      round(mean(CV[i,stp_y:edp_y],na.rm=T),3), '(CV)', '\t', round(score[k],3), '\t',
                      round(obs[k]),'\t',round(ref[k]),'\n', sep=''))
      }
    }
  }
  
  calculate_fd <- function(data,reference,sex){
    
    data[is.na(reference)] = NA
    
    fd = t(sapply(1:24,function(x) data[x,]/reference[x,]))

    write.table(fd,output_fd,col.names = FALSE,row.names = T,sep = "\t",quote=FALSE)
    
 
    fd.data = fd * 2
    if(sex == 'M'){fd.data[23:24,] = fd[23:24,]}
    
    return(list(data = data, fd.data = fd.data, reference = reference))
  }

  #----------------------------------------   get all necessary data

  data = as.matrix(alldata$data)
  sampleName = alldata$sampleName
  sumCount = alldata$sumCount
  sex = alldata$sex
  dyn.load(parameter$HeapMedian);

  #-----------------------------------------  manage output path

  if(parameter$IT_version){
    outdir_cnv = output_dir
    outdir_view = output_dir
    out_for_soft_data = paste(output_dir, '/', sampleName, '.view', sep='')
    output_fd = paste(output_dir,paste(sampleName,'.fd',sep = ""),sep = "/")
    out_cnv_result = paste(outdir_cnv, '/', sampleName, '.cnv', sep='')
    view.file = output_dir
  }else{
    outdir_cnv = paste(output_dir,'callcnv',sep='/')
    outdir_view = paste(output_dir,'result',sep='/')
    out_for_soft_data = paste(outdir_view, '/', sampleName, '.view.txt', sep='')
    output_fd = paste(output_dir,'preProcess',paste(sampleName,'.fd.txt',sep = ""),sep = "/")
    out_cnv_result = paste(outdir_cnv, '/', sampleName, '.cnv.txt', sep='')
    graph_out = paste(outdir_cnv, '/', sampleName, '.cnv.pdf', sep='')
    html_out = paste(outdir_cnv, '/', sampleName, '.cnv.html', sep='')
    view.file = paste(output_dir,'result','graph',sep="/")
  }

  
  #------------------------------------------ calculate fd value and Score
  
  reference[reference < 0.001] = NA
  
  homo_ratio = 0.1
  
  if(sex == 'M'){reference = reference[1:24,]
  }else{reference = rbind(reference[1:22,],reference[25:26,])}
  
  data_info = calculate_fd(data,reference,sex)
  
  data = data_info$data
  fd.data = data_info$fd.data
  reference = data_info$reference

  biaozhuncha = sd(fd.data[1:22,],na.rm=T)

  #----------------------------------- huaiyu code
  
  ###### set middle of windows as x coordinater
  
  x.coor = (seq(1, length(fd.data[1,]))-1)*nStep*window_len + nBlock*window_len/2 
  xcoor.offset = nBlock*window_len/2
  
  ###### set a list for save final result line
  result.lines = list()
  
  options(scipen=20)
  
  if (file.exists(out_cnv_result)) unlink(out_cnv_result)
  wave.degree = sum(apply(fd.data, 1, sd, na.rm=T), na.rm=T)
  cat(file=out_cnv_result, paste(strsplit(basename(out_cnv_result), '[.]', perl=T)[[1]][1],
                                 sex, wave.degree, sep='\t'), end='\n')
  
  #------------------------------------ GDSW algorithm
  
  if(parameter$IT_version){
    result.lines = lapply(1:24,function(x) Call_cnv_GD(x.coor, fd.data[x, ],data[x,],reference[x,], x, sex))
  }else{
    #pdf(graph_out, width=12, height=9,useDingbats=FALSE)
    #par(mfrow=c(3,1), lend=1)
    result.lines = lapply(1:24,function(x) Call_cnv_GD(x.coor, fd.data[x, ],data[x,],reference[x,], x, sex))
    #dev.off()
  }
  
  
  
  #-------------------------------------  generate view files
  
  if (file.exists(out_for_soft_data)) {unlink(out_for_soft_data)}
  param.cal = rep(2, 24)
  param.cal[24] = ifelse(sex == 'M', 1, 2)
  param.cal[23] = ifelse(sex == 'M', 1, 2)
  
  copyNumber = fd.data
  
  # print(result.lines)
  finalLine = as.matrix(do.call(rbind, result.lines))
  
  view = array(0,dim=c(49,Nbin))
  view[1,] = x.coor
  
  #----------------------------------------- plot png
  
  Out.put = paste(view.file, '/', sampleName, '.mht.png', sep='')
  Out.put.noxy = paste(view.file, '/', sampleName, '.mht.noXY.png', sep='')

  for (j in c(1:24)) {
    view[2*j,] = copyNumber[j, ]
    view[2*j+1,] = finalLine[j, ]
  }
  
  mhtplot(Out.put, view, sampleName, 'XY')
  mhtplot(Out.put.noxy, view, sampleName, 'noXY')
  scatterplot(view.file, view, sampleName, sex, ideoCyto)
  
  if(!parameter$IT_version){
    info = list()
    info$view = view
    info$sampleName = sampleName
    info$sex = sex
    info$ideoCyto = ideoCyto
#    rmarkdown::render(html_plot,output_file = html_out,params = info,intermediates_dir = tempdir(),
#                    html_document(fig_width = 18,fig_height = 4.5))
  }
  
  #-----------------------------------------  view files
  
  cat(paste(c('xcoor', x.coor), sep='\t'), file=out_for_soft_data, end='\n')
  
  for (j in c(1:24)) {
    cat(file=out_for_soft_data, paste(c(paste('chr', j, '.CN', sep=''), copyNumber[j, ]), sep=' '),append=T, end='\n')
    cat(file=out_for_soft_data, paste(c(paste('chr', j, '.L', sep=''), finalLine[j, ]), sep=' '),append=T, end='\n')
  }
  
  return(list(sd = biaozhuncha, sampleName = sampleName, sumCount = sumCount))
}

do_preProcess <- function(filepath_20K,parameter){

  correct_sex_type <- function(data,reference,nStep){
    temp = data[24,] / reference[24,]
    temp = temp[!is.na(temp) & temp > 0.1]
    return(ifelse(length(temp) > (20 + 80 / nStep),"M","F"))
  }
  
  window_len <- 200000
  limit_cutoff <- 0.2
  
  library(stringr)
  output_dir = parameter$output_dir
  chrs = parameter$chrs
  nBlock = parameter$nBlock
  nStep = parameter$nStep
  reference = parameter$reference
  reference[reference < 1] = NA
  
  if(parameter$IT_version){
    outpath = output_dir
  }else{
    outpath = paste(output_dir,'preProcess',sep='/')
  }
  
  # read file
  data = as.matrix(read.table(filepath_20K, row.names=1, stringsAsFactors=F, header=F))
  sampleName = strsplit(basename(filepath_20K),split='\\.')[[1]][1]
  sumCount = sum(data, na.rm=T)
  
  sex = correct_sex_type(data,reference,nStep)
  sex_info = paste(sampleName,'\t',sex,'\t',ratio,'%=',sum(data[24,],na.rm=T),'/',sum(data[1:22,],na.rm=T),sep='')

  if(sex == 'F'){
    reference = rbind(reference[1:22,],reference[25:26,])
  }else{
    reference = reference[1:24,]
  }
  
  index = 1:22
  alldata = data
  
  return(list(data = alldata, sex = sex, sampleName = sampleName, sumCount = sumCount))
}

output_copy <-function(alldata){
  for(i in 1:length(alldata)){
    #cat(alldata[[i]]$sampleName,append = T,end = '\n',file = parameter$output_copy)
    write.table(alldata[[i]]$ob_chr_dosage,col.names = F,quote=F,append = T,file = parameter$output_copy)
    #cat('\n',append = T,file = parameter$output_copy)
  }
}

output_distance <- function(alldata){
  
  SampleName = sapply(alldata,function(x) x$sampleName)
  DeAR_distance = sapply(alldata,function(x) x$DeAR_distance)
  GC_distance = sapply(alldata,function(x) x$GC_distance)
  aneuploid_distance = sapply(alldata,function(x) x$aneuploid_distance)
  out_data = cbind(SampleName,GC_distance,aneuploid_distance,DeAR_distance)
  write.table(out_data,col.names = T,quote=F,row.names = F,sep = "\t",file = parameter$output_DeAR)
  
}

output_QC <- function(x){
  
  SampleName = x$sampleName
  sumCount = x$sumCount
  biaozhuncha = x$sd
  out_data = cbind(SampleName,sumCount,biaozhuncha)
  colnames(out_data) = c('SampleName','sumCount','sd')
  write.table(out_data,col.names = T,quote=F,row.names = F,sep = "\t",file = parameter$output_QC)
  
}

message("[read data and do preProcess ...]")
reference_all = as.matrix(read.table(parameter$reference_file, row.names=1, stringsAsFactors=F, header=F))
parameter$reference = reference_all[1:26,]
parameter$CV = reference_all[27:52,]

ideoCyto = read.table(parameter$cytoband_file, header=F, stringsAsFactors=F)

filepath_20K <- parameter$file_20K

alldata = do_preProcess(filepath_20K, parameter)


message('[Call CNV ...]')

info = detectCNV_core(alldata,parameter,ideoCyto)

output_QC(info)



