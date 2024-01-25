/share/work1/wangrr/local/miniconda3/lib/R/bin/Rscript t.R
#基因注释，结果文件my_interst_gene_list.xls
perl 03.merge_cnv_gistic.pl 

#hebing amp del
qsub -cwd -l vf=2G run_gisticChromPlot.sh
