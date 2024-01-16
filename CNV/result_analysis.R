
Args <- commandArgs(TRUE)
if (length(Args) != 3) {
  cat("usage: <infile> <outdir> <reference_file>\n", file=stderr())
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

#program_dir = getScriptPath()

library('stringr')

program_dir = getScriptPath()
cyto_file = "cytoBand_hg19_ucsc.txt"
IT_version = FALSE

inpath = Args[1]
outpath = Args[2]
reference_file = Args[3]
cyto_file = paste(program_dir,'cytoBand_hg19_ucsc.txt',sep='/')

merge_output = paste(outpath,'callcnv.merge_out.txt',sep = '/')
detail_output = paste(outpath,'cnv.detail.txt',sep = '/')
low_ratio_output = paste(outpath,'low_ratio_aneuploid.txt',sep = '/')

if(file.exists(merge_output)) unlink(merge_output)
if(file.exists(detail_output)) unlink(detail_output)
if(file.exists(low_ratio_output)) unlink(low_ratio_output)

parameter = list()

#parameter$Statistics_available_fd_min = 0.05;         
parameter$cnv_fd_cutoff = 0.05; 
parameter$minimum_cnv_length = 2e6                  ## minimum cnv length cutoff
parameter$cnv_mosaic_cn_cutoff = 0.03;               ## cnv mosaic copy number cutoff
parameter$cnv_mosaic_length_cutoff = 5e6           ## cnv mosaic length cutoff
parameter$mosaic_fd_cutoff = 0.02;                   ## mosaic fd cutoff, used for detect mosaic chromosome
parameter$mosaic_region_per_chr_cutoff = 0.8;       ## when mutation region  > 90% chromosome, is a ploid mutation
## if a cnv region content of this ratio NA bin, will remove
parameter$NA_max_ratio_1 = 20;              ## when size 
parameter$NA_max_ratio_2 = 30;              ## when size 
parameter$NA_max_ratio_3 = 50;              ## when size 

parameter$chr_length = c(
  249250621, 243199373, 198022430, 191154276, 180915260, 
  171115067, 159138663, 146364022, 141213431, 135534747, 
  135006516, 133851895, 97269879, 89749541, 83531393, 
  90354753, 81195210, 78077248, 59128983, 63025520, 
  34929896, 36604567, 155270560, 59373566)
parameter$chr_centromere_start = c(
  121535434,92326171,90504854,49660117,46405641,
  58830166,58054331,43838887,47367679,39254935,
  51644205,34856694,16000000,16000000,17000000,
  35335801,22263006,15460898,24681782,26369569,
  11288129,13000000,58632012,10104553)
parameter$chr_centromere_end = c(
  124535434,95326171,93504854,52660117,49405641,
  61830166,61054331,46838887,50367679,42254935,
  54644205,37856694,19000000,19000000,20000000,
  38335801,25263006,18460898,27681782,29369569,
  14288129,16000000,61632012,13104553)
parameter$yangshui_mosaic_cutoff = c(
  0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
  0.05,0.1,0.1,0.05,0.15,0.05,0.05,0.1,0.1,0.1)
parameter$liuchanwu_mosaic_cutoff = c(
  0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
  0.1,0.15,0.2,0.1,0.3,0.1,0.1,0.2,0.15,0.1)

chrs = paste('chr',1:22,sep='')
parameter$chrs = c(chrs,'chrX','chrY')

read_cytoBand <- function(){
  
  cyto_info = list()
  
  chr = c()
  spos = c()
  epos = c()
  cyto = c()
  
  fr <- file(cyto_file, "r")
  line = readLines(fr,n=1)
  while(TRUE){
    line = readLines(fr,n=1)
    if(length(line) == 0) break
    info = strsplit(line,'\t')[[1]]
    chr = c(chr,info[1])
    spos = c(spos,as.numeric(info[2]))
    epos = c(epos,as.numeric(info[3]))
    cyto = c(cyto,info[4])
  }
  close(fr)
  
  cyto_info$chr = chr
  cyto_info$spos = spos
  cyto_info$epos = epos
  cyto_info$cyto = cyto

  return(cyto_info)
}

read_cnv_result <- function(file_path,cytoBand){
  
  detectBand <- function(chr,spos,epos,cytoBand){
    
    spos_add = TRUE
    epos_add = TRUE
    
    for(i in 1:length(cytoBand$chr)){
      if(cytoBand$chr[i] != chr) next
      if(spos >= cytoBand$spos[i] & spos < cytoBand$epos[i] & spos_add){
        cyto_chr = stringr::str_extract_all(cytoBand$chr[i], "[XY\\d]+")[[1]]
        cyto_spos = cytoBand$cyto[i]
        spos_add = FALSE
      }
      if(epos > cytoBand$spos[i] & epos <= cytoBand$epos[i] & epos_add){
        cyto_epos = cytoBand$cyto[i]
        epos_add = FALSE
      }
    }
    
    if(cyto_spos == cyto_epos){
      sentence = paste(cyto_chr,cyto_spos,sep='')
    }else{
      sentence = paste(cyto_chr,cyto_spos,'-',cyto_epos,sep='')
    }
    
    return(sentence)
    
  }
  
  file_info = list()
  
  file_info$sampleName = strsplit(basename(file_path),split='\\.')[[1]][1]
  
  fr <- file(file_path, "r")
  line = readLines(fr,n=1)
  file_info$sex = strsplit(line,split='\t')[[1]][2]
  
  NA_bins = c()
  all_bins = c()
  cnv_length = c()
  spos = c()
  epos = c()
  copy_number = c()
  chr = c()
  cyto = c()
  
  while(TRUE){
    line = readLines(fr,n=1)
    if(length(line) == 0) break
    info = stringr::str_extract_all(line, "[chrXchrYchr\\d\\.]+")[[1]]
    chr = c(chr,info[1])
    spos = c(spos,as.numeric(info[2]))
    epos = c(epos,as.numeric(info[3]))
    copy_number = c(copy_number,as.numeric(info[4]))
    #NA_bins = c(NA_bins,as.numeric(info[8]))
    NA_bins = c(NA_bins,0)
    all_bins = c(all_bins,as.numeric(info[5]))
    cnv_length = c(cnv_length,as.numeric(info[3])-as.numeric(info[2]))
    temp = detectBand(info[1],as.numeric(info[2]),as.numeric(info[3]),cytoBand)
    cyto = c(cyto,temp) 
  }
  
  close(fr)
  
  file_info$NA_bins = NA_bins
  file_info$all_bins = all_bins
  file_info$cnv_length = cnv_length
  file_info$spos = spos
  file_info$epos = epos
  file_info$copy_number = copy_number
  file_info$chr = chr
  file_info$cyto = cyto
  
  return(file_info)
  
}

output_cnv_merge_out <- function(all_cnv_info){
  
  options(scipen = 200)
  for(i in 1:length(all_cnv_info)){
    
    info = all_cnv_info[[i]]
    cat(file=merge_output, append=T,paste(info$sampleName,info$sex,'\n',sep = '\t'))
    chr = info$chr
    spos = info$spos
    epos = info$epos
    all_bins = info$all_bins
    NA_bins = info$NA_bins
    cnv_length = info$cnv_length
    copy_number = info$copy_number
    for(j in 1:length(chr)){
      cat(file=merge_output, append=T,paste(chr[j],':',spos[j],'-',epos[j],'_',copy_number[j],
                '_',all_bins[j],'_',NA_bins[j],'_',cnv_length[j],'\n',sep = ''))
    }
    
  }
  
}

output_files <- function(cnv_info_output){
  
  options(scipen = 200)
  cat(file=detail_output,paste('Sample','Final_result','Detail_result',sep = '\t'),end='\n')
  cat(file=low_ratio_output,paste('Sample','Detail_result','Mosaic_type','Mosaic_ratio',sep = '\t'),end='\n')
  for(i in 1:length(cnv_info_output)){
    sentence_detail = cnv_info_output[[i]][1]
    sentence_low_ratio = cnv_info_output[[i]][2]
    cat(file=detail_output, append=T,sentence_detail,end='\n')
    if(sentence_low_ratio != 'none') cat(file=low_ratio_output, append=T,sentence_low_ratio,end='\n')
  }
  
}

manage_output <- function(cnv_info,parameter){
  
  chromosome_status <- function(cnv_info,uselen,min_mosaic_cutoff,chr)
  {
    
    chromosome_status <- new.env()
    chromosome_status$sex = cnv_info$sex
    chromosome_status$uselen = uselen
    chromosome_status$chr = chr
    chromosome_status$spos = c()
    chromosome_status$epos = c()
    chromosome_status$copy_number = c()
    chromosome_status$NA_bins = c()
    chromosome_status$all_bins = c()
    chromosome_status$cnv_length = c()
    chromosome_status$cyto = c()
    chromosome_status$ctype = c()   # cnv type
    chromosome_status$num_cnv = 0
    chromosome_status$type = 'cnv'
    
    if(chromosome_status$sex == 'M' && 
       (chromosome_status$chr == 'chrX' | chromosome_status$chr == 'chrY')){
      chromosome_status$yh = 1
    }else{
      chromosome_status$yh = 2
    }
    
    chromosome_status$add_cnv <- function(cnv_info,j){
      chromosome_status$spos = c(chromosome_status$spos,cnv_info$spos[j])
      chromosome_status$epos = c(chromosome_status$epos,cnv_info$epos[j])
      chromosome_status$copy_number = c(chromosome_status$copy_number,cnv_info$copy_number[j])
      chromosome_status$NA_bins = c(chromosome_status$NA_bins,cnv_info$NA_bins[j])
      chromosome_status$all_bins = c(chromosome_status$all_bins,cnv_info$all_bins[j])
      chromosome_status$cnv_length = c(chromosome_status$cnv_length,cnv_info$cnv_length[j])
      chromosome_status$cyto = c(chromosome_status$cyto,cnv_info$cyto[j])
      chromosome_status$num_cnv = length(chromosome_status$copy_number)
      temp = ifelse(cnv_info$copy_number[j]>chromosome_status$yh,'dup','del')
      chromosome_status$ctype = c(chromosome_status$ctype,temp)
    }
    
    chromosome_status$deter_chromosome_mosaic <- function(){
      # determine whether mosaic, aneuploid or cnv
      
      temp = sort(chromosome_status$copy_number,index = TRUE)
      copy = temp$x
      index = temp$ix
      NA_bins = chromosome_status$NA_bins[index]
      all_bins = chromosome_status$all_bins[index]
      
      for(i in 1:chromosome_status$num_cnv){
        if(i == 1){
          mosaic_class = 1
        }else{
          if(copy[i] - copy[i-1] < parameter$mosaic_fd_cutoff & 
             ((copy[i]>chromosome_status$yh & copy[i-1]>chromosome_status$yh)|
              (copy[i]<chromosome_status$yh & copy[i-1]<chromosome_status$yh))){
            mosaic_class = c(mosaic_class,mosaic_class[i-1])
          }else{
            mosaic_class = c(mosaic_class,mosaic_class[i-1]+1)
          }
        }
      }
      for(i in 1:max(mosaic_class)){
        seg_index = mosaic_class == i
        useratio = (sum(all_bins[seg_index]) - sum(NA_bins[seg_index]))/chromosome_status$uselen
        useratio_cutoff = parameter$mosaic_region_per_chr_cutoff
        if(chromosome_status$chr == 'chrY') useratio_cutoff = useratio_cutoff - 0.1
        if(useratio > useratio_cutoff){
          chromosome_status$type = 'mosaic'
          chromosome_status$mosaic_copy = round(sum(copy[seg_index] * 
              (all_bins[seg_index] - NA_bins[seg_index])) / sum(all_bins[seg_index] - NA_bins[seg_index]),3)
        }
      }
    }
    
    chromosome_status$deter_output <- function(){
      # manage output information for single chromosome
      
      if(chromosome_status$type == 'mosaic'){
        shift_value = chromosome_status$mosaic_copy - chromosome_status$yh
        if(abs(shift_value) <= parameter$mosaic_fd_cutoff){
          chromosome_status$type = 'cnv'
          if(abs(shift_value) >= min_mosaic_cutoff){
            chromosome_status$type = 'cnv and low_ratio_mosaic'
            chromosome_status$mosaic_ctype = ifelse(shift_value > 0,'dup','del')
          }
        }else{
          cutoff = parameter$mosaic_fd_cutoff
          region = c(c(-2,-1,-0)-cutoff,c(0:100)+cutoff)
          for(i in 1:(length(region)-1)){
            if(shift_value>region[i] & shift_value<=region[i+1]){
              chromosome_status$ctype = i-3
            }
          }
          chromosome_status$type = 'high_ratio_mosaic'
          if(chromosome_status$yh == 1){
            aneuploid_region = c(2,3,4)
          }else{
            aneuploid_region = c(1,3,4)
          }
          for(i in 1:length(aneuploid_region)){
            if(chromosome_status$mosaic_copy >= aneuploid_region[i] - parameter$mosaic_fd_cutoff & 
                chromosome_status$mosaic_copy <= aneuploid_region[i] + parameter$mosaic_fd_cutoff){
              chromosome_status$type = 'aneuploid'
            }
          }
        }
      }
      if(chromosome_status$type == 'cnv' | chromosome_status$type == 'cnv and low_ratio_mosaic'){
        
        true_cnv = c()
        for(i in 1:chromosome_status$num_cnv){
          shift_value = abs(chromosome_status$copy_number[i] - chromosome_status$yh)
          if(shift_value < parameter$cnv_mosaic_cn_cutoff) 
            true_cnv = c(true_cnv,FALSE)
          else if(chromosome_status$cnv_length[i] <= parameter$minimum_cnv_length) true_cnv = c(true_cnv,FALSE)
          else if(chromosome_status$all_bins[i]-chromosome_status$NA_bins[i] > parameter$NA_max_ratio_1) true_cnv = c(true_cnv,FALSE)
          else if(chromosome_status$all_bins[i]-chromosome_status$NA_bins[i] > parameter$NA_max_ratio_2) true_cnv = c(true_cnv,FALSE)
          else if(chromosome_status$all_bins[i]-chromosome_status$NA_bins[i] > parameter$NA_max_ratio_3) true_cnv = c(true_cnv,FALSE)
          else if(chromosome_status$cnv_length[i] <= parameter$cnv_mosaic_length_cutoff & shift_value < parameter$cnv_fd_cutoff)
            true_cnv = c(true_cnv,FALSE)
          else true_cnv = c(true_cnv,TRUE)
        }
        
        chromosome_status$spos = chromosome_status$spos[true_cnv]
        chromosome_status$epos = chromosome_status$epos[true_cnv]
        chromosome_status$copy_number = chromosome_status$copy_number[true_cnv]
        chromosome_status$NA_bins = chromosome_status$NA_bins[true_cnv]
        chromosome_status$all_bins = chromosome_status$all_bins[true_cnv]
        chromosome_status$cnv_length = chromosome_status$cnv_length[true_cnv]
        chromosome_status$cyto = chromosome_status$cyto[true_cnv]
        chromosome_status$ctype = chromosome_status$ctype[true_cnv]
        chromosome_status$num_cnv = length(chromosome_status$copy_number)
        
        chromosome_status$cnv_mosaic = c()
        if(chromosome_status$num_cnv > 0){
          for(i in 1:chromosome_status$num_cnv){
            shift_value = abs(chromosome_status$copy_number[i] - chromosome_status$yh)
            chromosome_status$shift_value = shift_value
            if(shift_value >= parameter$cnv_mosaic_cn_cutoff & shift_value < 1 - parameter$cnv_mosaic_cn_cutoff &
               chromosome_status$cnv_length[i] > parameter$cnv_mosaic_length_cutoff){
              chromosome_status$cnv_mosaic = c(chromosome_status$cnv_mosaic,'(mos)')
            }else{
              chromosome_status$cnv_mosaic = c(chromosome_status$cnv_mosaic,'')
            }
          }
        }
      }
    }
    
    return(chromosome_status)
  }
  
  if(is.null(cnv_info$copy)){
    return(c(paste(cnv_info$sampleName,'46,XN','--',sep='\t'),'none'))
  }
  
  temp = stringr::str_extract_all(cnv_info$sampleName, "Y-|Y_|P-|P_|YM_|YM-|QM_|QM-|RM_|RM-|RZ_|RZ-|DM_|DM-|D_|D-|XM_|XM-|X_|X-")[[1]]
  if(length(temp) == 0){
    min_mosaic_cutoff = parameter$liuchanwu_mosaic_cutoff
  }else{
    min_mosaic_cutoff = parameter$yangshui_mosaic_cutoff
  }
  
  chrs = parameter$chrs
  Final_result_cnv = c()
  Detail_result_cnv = c()
  Final_result_aneuploid = c()
  Detail_result_aneuploid = c()
  aneuploid_num = c()
  low_ratio_detail = c()
  low_ratio_type = c()
  low_ratio_ratio = c()
  
  for(i in 1:length(chrs)){
    chromosome <- chromosome_status(cnv_info,parameter$uselen[i],min_mosaic_cutoff[i],chrs[i]);
    for(j in 1:length(cnv_info$chr)){
      if(cnv_info$chr[j] == chrs[i]) chromosome$add_cnv(cnv_info,j)
    }
    if(chromosome$num_cnv == 0) next
    chromosome$deter_chromosome_mosaic()
    chromosome$deter_output()
    if(chromosome$type == 'cnv' | chromosome$type == 'cnv and low_ratio_mosaic'){
      # manage cnv and low ratio mosaic information
      if(chromosome$type == 'cnv and low_ratio_mosaic'){
        copy = chromosome$mosaic_copy 
        chr = chromosome$chr
        ctype = chromosome$mosaic_ctype
        yh = chromosome$yh
        low_ratio_detail = c(low_ratio_detail,paste(chr,'_',copy,sep=''))
        low_ratio_type = c(low_ratio_type,ctype)
        low_ratio_ratio = c(low_ratio_ratio,paste(format(abs(copy-yh)*100,nsmall=1),'%',sep=''))
      }
      if(length(chromosome$copy_number) == 0) next
      for(j in 1:chromosome$num_cnv){
        chr = chromosome$chr
        spos = chromosome$spos[j]
        epos = chromosome$epos[j]
        copy = chromosome$copy_number[j]
        ctype = chromosome$ctype[j]
        cnv_mosaic = chromosome$cnv_mosaic[j]
        cyto = chromosome$cyto[j]
        cnv_length = chromosome$cnv_length[j]
        options(scipen = 200)
        Detail_result_cnv = c(Detail_result_cnv,paste(chr,':',spos,'-',epos,'_',copy,'_',ctype,sep=''))
        copy_final = format(round(cnv_length/1e6,2),nsmall = 2)
        Final_result_cnv = c(Final_result_cnv,paste(cyto,'(',ctype,', ',copy_final,'Mb)',cnv_mosaic,sep=''))
      }
    }else{
      # manage aneuploid or mosaic information
      copy = chromosome$mosaic_copy 
      chr = chromosome$chr
      yh = chromosome$yh
      cyto_chr = stringr::str_extract_all(chr, "[XY\\d]+")[[1]]
      sign = ifelse(copy - yh > 0,'+','-')
      if(chromosome$type == 'high_ratio_mosaic'){
        Detail_result_aneuploid = c(Detail_result_aneuploid,paste(chr,'_',copy,sep=''))
        Final_result_aneuploid = c(Final_result_aneuploid,paste(sign,cyto_chr,'(mos)',sep = ''))
        aneuploid_num = c(aneuploid_num,chromosome$ctype)
      }else{
        Detail_result_aneuploid = c(Detail_result_aneuploid,paste(chr,'_',copy,sep=''))
        Final_result_aneuploid = c(Final_result_aneuploid,paste(sign,cyto_chr,sep = ''))
        aneuploid_num = c(aneuploid_num,chromosome$ctype)
      }
    }
  }
  
  if(!is.null(low_ratio_detail)){
    low_ratio_detail = paste(low_ratio_detail,collapse = ";")
    low_ratio_type = paste(low_ratio_type,collapse = ";")
    low_ratio_ratio = paste(low_ratio_ratio,collapse = ";")
    low_ratio_sentence = paste(cnv_info$sampleName,low_ratio_detail,low_ratio_type,low_ratio_ratio,sep='\t')
  }else{
    low_ratio_sentence = 'none'
  }
  
  Final_sentence = paste(c(Final_result_aneuploid,Final_result_cnv),collapse = ";")
  Detail_sentence = paste(c(Detail_result_aneuploid,Detail_result_cnv),collapse = ";")
  
  if(Detail_sentence == ''){
    # no cnv, mosaic, aneuploid
    return(c(paste(cnv_info$sampleName,'46,XN','--',sep='\t'),low_ratio_sentence))
  }else{
    if(length(Final_result_aneuploid) == 0){
      Final_sentence = paste('46,XN',Final_sentence,sep=';')
    }else{
      num = 46
      chang = length(Final_result_aneuploid)
      
      for(i in 1:length(Final_result_aneuploid)){
        num = num + aneuploid_num[i]
      }
      
      if(length(stringr::str_extract_all(Final_result_aneuploid[chang],'X|Y')[[1]])>0){
        Final_result_aneuploid = c(Final_result_aneuploid[chang],Final_result_aneuploid[-chang])
        aneuploid_num = c(aneuploid_num[chang],aneuploid_num[-chang])
        if(length(stringr::str_extract_all(Final_result_aneuploid[chang],'X|Y')[[1]])>0){
          Final_result_aneuploid = c(Final_result_aneuploid[chang],Final_result_aneuploid[-chang])
          aneuploid_num = c(aneuploid_num[chang],aneuploid_num[-chang])
        }
          
        min_for = ifelse(chang < 2,1,2)
        X = 0
        Y = 0
        for(i in 1:min_for){
          if(length(grep('X',Final_result_aneuploid[i]))>0){
            X = X + aneuploid_num[i]
          }
          if(length(grep('Y',Final_result_aneuploid[i]))>0){
            Y = Y + aneuploid_num[i]
          }
        }
        if(cnv_info$sex == 'F'){
          X = X + 2
        }else{
          X = X + 1
          Y = Y + 1
        }
        num_sentence = paste(num,paste(c(rep('X',X),rep('Y',Y)),collapse = ""),sep=',')
      }else{
        num_sentence = paste(num,'XN',sep=',')
      }
      Final_sentence_begin = paste(c(num_sentence,Final_result_aneuploid),collapse = ",")
      Final_sentence = paste(c(Final_sentence_begin,Final_result_cnv),collapse = ";")
    }
    return(c(paste(cnv_info$sampleName,Final_sentence,Detail_sentence,sep='\t'),low_ratio_sentence))
  }
}

if(IT_version){
  file_cnv = as.list(paste(inpath,list.files(inpath,".cnv"),sep='/'))
}else{
  file_cnv = as.list(paste(inpath,list.files(inpath,".cnv.txt"),sep='/'))
}

file_cnv_list = list()

cytoBand = read_cytoBand()
reference = as.matrix(read.table(reference_file, row.names=1, stringsAsFactors=F, header=F))
#parameter$uselen = sapply(1:24,function(x) length(reference[x,reference[x,]>0]))

parameter$uselen = c(28585,28055,23717,20524,19412,25683,19045,16336,3373,16809,18040,17186,11718,11427,
                     11343,12219,12277,9633,9805,8487,4761,5527,12576,1185)

all_cnv_info = lapply(file_cnv,read_cnv_result,cytoBand)
output_cnv_merge_out(all_cnv_info)
cnv_info_output = lapply(all_cnv_info,manage_output,parameter)
output_files(cnv_info_output)




