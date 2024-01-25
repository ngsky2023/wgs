#Download 1000 Genomes Phase3 Genotype Data: 
#https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1
#https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1
#https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1

bin=`dirname $0`
NewName=$bin/NewChrName.txt
plink2=/share/work2/jiangdz/software/plink2/plink2
plink=/share/work2/jiangdz/software/plink2/plink
1000g_plink=/share/work2/jiangdz/software/plink2/test/


cd $1000g_plink 
###1.Remove non-biallelic variants , Set variant id ,  Remove dup variants , Keep snp variants for  1000g 
$plink2 --pfile all_phase3_ns --snps-only --make-bed --out 1000g --threads 8 --max-alleles 2 --set-all-var-ids 'chr'@:#:\$r:\$a --rm-dup force-first && echo done1

###2.Remove High LD for 1000g snps
$plink --bfile 1000g --make-set $bin/High_LD.txt --write-set --out hild  -allow-extra-chr && echo done2
$plink --bfile 1000g --make-bed --exclude hild.set  --out 1000g.hild  -allow-extra-chr && echo done3

###3.LD for 1000g snps 
$plink --bfile 1000g.hild --indep-pairwise 50 10 0.2 -out 1000g-LD --allow-extra-chr --make-bed --threads 16 --maf 0.05 --geno 0.1  && echo done4
$plink  --bfile 1000g.hild --extract 1000g-LD.prune.in --make-bed --out 1000g-LD-0.1-0.05 --allow-extra-chr --const-fid family_id  && echo done5

###4.Bim2vcf and Annotation
cut -f 2 1000g-LD-0.1-0.05.bim |sed 's/:/\t/g' |awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$1":"$2":"$3":"$4,$3,$4,60,"PASS","MQ=25"}' |grep -v PAR |cat $bin/vcf_head.txt - > 1000g-LD-0.1-0.05.bim.vcf
perl /share/work1/capsmart_rd/v3.15.0/capSMART2.0/csmart/mutdect/scripts/Annovar/2016Feb01/convert2annovar.pl --includeinfo --format vcf4  --outfile 1000g-LD-0.1-0.05.bim.avinput 1000g-LD-0.1-0.05.bim.vcf && \
perl /share/work1/capsmart_rd/v3.15.0/capSMART2.0/csmart/mutdect/scripts/Annovar/2016Feb01/table_annovar.pl --buildver hg19 --thread 8 --remove --otherinfo --protocol refGene,popfreq_max_20150413,avsnp150,1000g2015aug_all,1000g2015aug_eas -operation g,f,f,f,f  -nastring . 1000g-LD-0.1-0.05.bim.avinput /share/Onc_Soft_DB/database/annovar/humandb --outfile 1000g-LD-0.1-0.05.bim.anno && echo done6

###5.Process of PUCH Cohort
for i in `ls run/snv-indel/germ/strelka/*`
do
 sample=`basename $i`

 #Select variant in chr1-22,chrX,chrY,chrM for PUCH Cohort, and DP>=10 && GQ>=30
 echo "/share/public/software/bcftools-1.9/bin/bcftools view -v snps -i \"FORMAT/DP>=10 & FORMAT/GQ>=30\" -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM $i/variants.vcf.gz |/share/public/software/bcftools-1.9/bin/bcftools annotate --rename-chrs $NewName  - -Oz -o $i/variants.chr.vcf.gz --threads 8 && echo done" >  plink_preprocess_$sample.sh 

 #Remove non-biallelic variants , Set variant id ,  Remove dup variants, Add sex for PUCH Cohort
 echo "/share/work2/jiangdz/software/plink2/plink2 --vcf $i/variants.chr.vcf.gz --make-bed --out $i/$sample --threads 8 --max-alleles 2 --set-all-var-ids 'chr'@:#:\\\$r:\\\$a --rm-dup force-first  --fam $i/$sample.fam --split-par hg19 && echo done" >>plink_preprocess_$sample.sh
 qsub -cwd -q oncrd.q plink_preprocess_$sample.sh
done

$plink --merge-list $bin/beds.txt --make-bed --threads 16 -out merged --allow-extra-chr && echo done
$plink --bfile merged -out merged-0.05-0.1 --allow-extra-chr --make-bed  --threads 16 --maf 0.05 --geno 0.1

###6.Select SNPs
awk 'BEGIN{FS=OFS="\t"}NR>1{if($14>=0.85){print $17}}' 1000g-LD-0.1-0.05.bim.anno.hg19_multianno.txt >keep_site.pca
$plink -bfile merged-0.05-0.1 --extract keep_site.pca --make-bed -out merged-0.05-0.1-keep --allow-extra-chr
cut -f 2 merged-0.05-0.1-keep.bim > keep_site_overlap.pca
$plink -bfile 1000g-LD-0.1-0.05 --extract keep_site_overlap.pca --make-bed -out 1000g-keep --allow-extra-chr
###7.PCA
$plink --bfile 1000g-keep --bmerge merged-0.05-0.1-keep  --pca -out pca --allow-extra-chr 
awk 'BEGIN{FS="[ \t]"}FNR==NR{a[$1]=$5}NR>FNR{if(a[$2]){$1=a[$2];print $0}else{print $0}}' all_phase3_ns.psam pca.eigenvec > pca.eigenvec.add
python3 $bin/pca.py pca.eigenvec.add pca.eigenval pca
