awk 'NR==FNR{a[$1]=$1}NR>FNR{if(a[$3]==$3)print $0}' cancer_gene_census_OncTSG_name.txt CNV_region_genelist.txt > `pwd`/result/CNV_driver_genelist.txt

/share/public/software/R/3.5.3/lib64/R/bin/Rscript 02.CNV_drivergenes.data.r
/share/public/software/R/3.5.3/lib64/R/bin/Rscript genes.conf_95.R 0.05
