#!/bin/bash

###------通过Mutsigcv获得显著的driver gene并绘制瀑布图
file_path=/share/Onc_KYproject/SharingPlatform/Analysis/WES/12.driver_gene/01.Mutsigcv
file_path_data=/share/Onc_KYproject/SharingPlatform/Analysis/WES/12.driver_gene/01.Mutsigcv/test
Rscript=/share/work1/sum4489/software/miniconda3/envs/R-3.6.1/bin/Rscript


###----Step1:去除假阳性位点 输入：流程跑的结果result_Somatic.xls,PoN.xls；result_Somatic.xls是每个样本下的snv结果合并文件；PoN.xls是共有的  输出：mutation_delFP.xls 
##awk会根据分隔符将行分成若干个字段，$0为整行，$1为第一个字段，$2 为第2个地段，依此类推 
##OFS='-'是代表着前面的$2,$3,$4,$5,$6,$13,$15以-结合输出，OFS列输出分隔符，以\t分列$0
###NR和FNR都是记录读取数据的行数；NR,表示awk开始执行程序后所读取的数据行数.FNR,与NR功用类似,不同的是awk每打开一个新文件,FNR便从0重新累计.
##所以当是读取的文件是PoN.xls.xls这个文件的内容的时候，把它的文件内容放到a数组的第一列里，
#当在读取第二个文件时，即读取result_Somatic.xls文件时，判断，这个result_Somatic.xls的第一列是否在这个假阳性位点文件里，是的话就不输出，不是的话输出result_Somatic.xls内容。、
##这里$2,$3,$4,$5,$6,$13,$15分别指的是Chr,Start,End,Ref,Alt,Gene.refGene,ExonicFunc.refGene，如果这几列在你的文件里不是在对应的这几列，这里需要更改。
###得到result_Somatic.xls的所有表头
awk '{print $2,$3,$4,$5,$6,$13,$15"\t"$0}' OFS='-' ${file_path_data}/result_Somatic.xls |awk 'NR==FNR {a[$1];next} NR>FNR {if($1 in a) next; else print $0}' ${file_path_data}/PoN.xls - |cut -f 2-138 >> ${file_path_data}/mutation_delFP.xls

head -n 1 ${file_path_data}/result_Somatic.xls >${file_path_data}/title.xls

##给mutation_delFP.xls添加表头
cat ${file_path_data}/title.xls ${file_path_data}/mutation_delFP.xls >> ${file_path_data}/result_Somatic_Filter.xls

###-----step2:转化为Mutsigcv可读取的格式.maf
perl ${file_path}/maf_format.pl --i ${file_path_data}/result_Somatic_Filter.xls --o ${file_path_data}/result_Somatic_Filter.maf

###-----step3:运行Mutsigcv
###Mutsigcv输入文件除了maf文件是我们自己的，其他的都是可以从mutsigcv软件中提取
#其中exome_full192.coverage.txt 覆盖度文件，gene.covariates.txt 协变量文件，mutation_type_dictionary_file.txt 突变类型字典文件,chr_files_hg19参考基因组
#./mutsigcv这个是输出文件都以mutsigcv开头

file_MutSigCV_path=/share/Onc_KYproject/SharingPlatform/Scripts/mutsigcv/software/

${file_MutSigCV_path}/MutSigCV_1.4/run_MutSigCV.sh ${file_MutSigCV_path}/v81 ${file_path_data}/result_Somatic_Filter.maf ${file_MutSigCV_path}/exome_full192.coverage.txt ${file_MutSigCV_path}/gene.covariates.txt  ${file_path_data}/mutsigcv ${file_MutSigCV_path}/mutation_type_dictionary_file.txt ${file_MutSigCV_path}/chr_files_hg19

###-----step4:提取结果文件中的基因名和P值
#FS( Field Separator ) : 列分割符。决定了怎么将一行划分为几段
awk -v FS='\t' -v OFS='\t' '{print $1,$14}' ${file_path_data}/mutsigcv.sig_genes.txt > ${file_path_data}/mutsigcv.sig_genes.new.txt

###定义num（样本量）
##expr求表达式变量的值
num=`cat ${file_path_data}/TMB.xls | wc -l`
sample_num=`expr $num - 1`

###----step5:获取瀑布图的输入文件F1_MutGene_Matrix.txt, F2_MutGene_Stats.txt,F3_NumMut_Sorted.txt 这里TopGene取30
perl ${file_path}/toF12_MutGene_All.pl --i ${file_path_data}/result_Somatic_Filter.xls --tmb ${file_path_data}/TMB.xls --sig ${file_path_data}/mutsigcv.sig_genes.new.txt --num $sample_num -topgene 30  --o ${file_path_data}

###----step6:画瀑布图
$Rscript ${file_path}/MutationOncoplot.R --file1 ${file_path_data}/F1_MutGene_Matrix.txt --file2 ${file_path_data}/F2_MutGene_Stats.txt --file3 ${file_path_data}/F3_NumMut_Sorted.txt --out ${file_path_data}/MutationOncoplot_ConsiderDriverGene.pdf







#==================================================================
#Author:sum4489
#email:sum4489@berryoncology.com
#Date： 2021-02-08
#FileName： run.sh
#===================================================================

