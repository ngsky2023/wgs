perl a-FigureFile_OCCC.pl
#awk -F '\t' '$3=="gene"' hg19.gff |awk -F '[\t=;]' '{print $1,$4,$5,$10}' OFS='\t' |awk -F '\t' 'NR==FNR {a[$2];next} NR>FNR {if($4 in a) print $1,$4,$2,$3}' OFS='\t' amp_ogtsg_rank.txt - >target_gene.txt
#awk -F '\t' '$3=="gene"' hg19.gff |awk -F '[\t=;]' '{print $1,$4,$5,$10}' OFS='\t' |awk -F '\t' 'NR==FNR {a[$2];next} NR>FNR {if($4 in a) print $1,$4,$2,$3}' OFS='\t' del_ogtsg_rank.txt - >>target_gene.txt
#awk -F '\t' '$8>0.5 || $13<(-0.15)' wgs_Gene_Chart.txt >wgs_Gene_Chart_top.txt
#perl target_gene_top.pl
