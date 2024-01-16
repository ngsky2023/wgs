#

###1.Removing non-biallelic variants , Setting variant id ,  Removing dup variants
/share/work2/jiangdz/software/plink2/plink2 --pfile $dir/all_phase3_ns --snps-only --make-bed --out 1000g --threads 8 --max-alleles 2 --set-all-var-ids 'chr'@:#:\$r:\$a --rm-dup force-first && echo done1

#remove High LD
/share/work2/jiangdz/software/plink2/plink --bfile /home/jiangdz/workdir/software/plink2/test_select/1000g --make-set /home/jiangdz/workdir/software/plink2/test_v4/High_LD.txt --write-set --out /home/jiangdz/workdir/software/plink2/test_select/LD/hild  -allow-extra-chr && echo done3
/share/work2/jiangdz/software/plink2/plink --bfile /home/jiangdz/workdir/software/plink2/test_select/1000g --make-bed --exclude /home/jiangdz/workdir/software/plink2/test_select/LD/hild.set  --out /home/jiangdz/workdir/software/plink2/test_select/LD/1000g.hild  -allow-extra-chr && echo done4

#LD
#/share/work2/jiangdz/software/plink2/plink --bfile /home/jiangdz/workdir/software/plink2/test_select/LD/1000g.hild --indep-pairwise 50 10 0.2 -out /home/jiangdz/workdir/software/plink2/test_select/LD/1000g-LD --allow-extra-chr --make-bed --threads 16 --maf 0.05 --geno 0.1  && echo done1
#/share/work2/jiangdz/software/plink2/plink  --bfile /home/jiangdz/workdir/software/plink2/test_select/LD/1000g.hild --extract /home/jiangdz/workdir/software/plink2/test_select/LD/1000g-LD.prune.in --make-bed --out /home/jiangdz/workdir/software/plink2/test_select/LD/1000g-LD-0.1-0.05 --allow-extra-chr --const-fid family_id  && echo done2

#select 
#注释1000g-LD-0.1-0.05.bim.vcf生成1000g-LD-0.1-0.05.bim.anno.hg19_multianno.txt  
#cut -f 2 1000g-LD-0.1-0.05.bim |sed 's/:/\t/g' |awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$1":"$2":"$3":"$4,$3,$4,60,"PASS","MQ=25"}' |grep -v PAR |cat vcf_head.txt - > 1000g-LD-0.1-0.05.bim.vcf

for i in 14
do 
awk -v b=$i 'BEGIN{FS=OFS="\t"}NR>1{if($b>=0.85){print $17}}' 1000g-LD-0.1-0.05.bim.anno.hg19_multianno.txt >keep_site.pca
/share/work2/jiangdz/software/plink2/plink -bfile ../Merge/merged-0.05-0.1 --extract keep_site.pca --make-bed -out merged-0.05-0.1-keep --allow-extra-chr 
cut -f 2 merged-0.05-0.1-keep.bim > keep_site_overlap.pca
/share/work2/jiangdz/software/plink2/plink -bfile 1000g-LD-0.1-0.05 --extract keep_site_overlap.pca --make-bed -out 1000g-keep --allow-extra-chr 
/share/work2/jiangdz/software/plink2/plink --bfile 1000g-keep --bmerge merged-0.05-0.1-keep  --pca -out pca_10000 --allow-extra-chr && echo done2 
awk 'BEGIN{FS="[ \t]"}FNR==NR{a[$1]=$5}NR>FNR{if(a[$2]){$1=a[$2];print $0}else{print $0}}' /home/jiangdz/workdir/software/plink2/test/all_phase3_ns.psam pca_10000.eigenvec > pca_10000.eigenvec.add
python3 pca.py pca_10000.eigenvec.add pca_10000.eigenval pca_$i
done

