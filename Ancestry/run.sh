bin=`dirname $0`
NewName=$bin/NewChrName.txt
plink2=/share/work2/jiangdz/software/plink2/plink2
plink=/share/work2/jiangdz/software/plink2/plink
1000g_plink=/share/work2/jiangdz/software/plink2/test/

###1.Remove non-biallelic variants , Set variant id ,  Remove dup variants , Keep snp variants for  1000g 
$plink2 --pfile $1000g_plink/all_phase3_ns --snps-only --make-bed --out $1000_plink/1000g --threads 8 --max-alleles 2 --set-all-var-ids 'chr'@:#:\$r:\$a --rm-dup force-first && echo done1

###2.Remove High LD for 1000g snps
$plink --bfile $1000_plink/1000g --make-set $bin/High_LD.txt --write-set --out $1000g_plink/hild  -allow-extra-chr && echo done2
$plink --bfile $1000_plink/1000g --make-bed --exclude $1000g_plink/hild.set  --out $1000g_plink/1000g.hild  -allow-extra-chr && echo done3

###3.LD 
$plink --bfile $1000g_plink/1000g.hild --indep-pairwise 50 10 0.2 -out $1000g_plink/1000g-LD --allow-extra-chr --make-bed --threads 16 --maf 0.05 --geno 0.1  && echo done4
$plink  --bfile $1000g_plink/1000g.hild --extract $1000g_plink/1000g-LD.prune.in --make-bed --out $1000g_plink/1000g-LD-0.1-0.05 --allow-extra-chr --const-fid family_id  && echo done


#1.
for i in `ls run/snv-indel/germ/strelka/*`
do
 sample=`basename $i`

 ###Select variant in chr1-22,chrX,chrY,chrM
 echo "/share/public/software/bcftools-1.9/bin/bcftools view -v snps -i \"FORMAT/DP>=10 & FORMAT/GQ>=30\" -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM $i/variants.vcf.gz |/share/public/software/bcftools-1.9/bin/bcftools annotate --rename-chrs $NewName  - -Oz -o $i/variants.chr.vcf.gz --threads 8 && echo done" >  plink_preprocess_$sample.sh 

 ###Remove non-biallelic variants , Set variant id ,  Remove dup variants, Add sex
 echo "/share/work2/jiangdz/software/plink2/plink2 --vcf $i/variants.chr.vcf.gz --make-bed --out $i/$sample --threads 8 --max-alleles 2 --set-all-var-ids 'chr'@:#:\\\$r:\\\$a --rm-dup force-first  --fam $i/$sample.fam --split-par hg19 && echo done" >>plink_preprocess_$sample.sh

 qsub -cwd -q oncrd.q plink_preprocess_$sample.sh
done


###merge and LD
ld=/home/jiangdz/workdir/software/plink2/test_v4/LD
#ls /share/work2/jiangdz/software/plink2/test/*/*.bed |sed 's/.bed//' > beds.txt 
echo "/share/work2/jiangdz/software/plink2/plink -bfile /home/jiangdz/workdir/software/plink2/test_v4/1000g --merge-list /home/jiangdz/workdir/software/plink2/test_v4/beds.txt --make-bed --threads 8 -out merged --allow-extra-chr && echo done " > LD/merge.sh
echo "/share/work2/jiangdz/software/plink2/plink --bfile $ld/merged --indep-pairwise 50kb 5 0.2 -out $ld/merged-LD --allow-extra-chr  --make-bed --threads 8 --maf 0.1 --geno 0.05  && echo done1" > LD/LD.sh
echo "/share/work2/jiangdz/software/plink2/plink  --bfile $ld/merged --extract $ld/merged-LD.prune.in --make-bed --out $ld/merged-LD-0.1-0.05 --allow-extra-chr --const-fid PUCH  && echo done2" >> LD/LD.sh

###PCA

echo "/share/work2/jiangdz/software/plink2/plink -bfile $ld/merged-LD-0.1-0.05 --make-set $ld/../High_LD.txt --write-set --out $ld/hild  -allow-extra-chr && echo done " > PCA/pca.sh
echo "/share/work2/jiangdz/software/plink2/plink -bfile $ld/merged-LD-0.1-0.05 --make-bed --exclude $ld/hild.set  --out $ld/merged-LD-0.1-0.05.hild  -allow-extra-chr && echo done1 " >> PCA/pca.sh
echo "/share/work2/jiangdz/software/plink2/plink -bfile $ld/merged-LD-0.1-0.05.hild --pca -out pca --allow-extra-chr && echo done2 " >> PCA/pca.shD
