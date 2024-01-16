sample=$1
Tbam=$2
Nbam=$3
#smoove
export PATH=/share/work1/muting/local/lumpy-sv/bin/:$PATH
export PATH=/share/work1/muting/local/:$PATH
export PATH=/share/work1/muting/pkg/hts/:$PATH
export PATH=/share/public/software/bcftools/1.4/:$PATH
export LD_LIBRARY_PATH=/share/work1/muting/pkg/hts:$LD_LIBRARY_PATH
/share/work1/wangrr/local/simple/bin/smoove call -x --name $sample --exclude /share/work1/wangrr/DB/hg19/btu356_LCR-hs37d5.bed --fasta /share/work1/wangrr/DB/hg19/hg19AddVirus.fa -p 2 --genotype $Tbam $Nbam

sh ./$sample-lumpy-cmd.sh | grep -v IMPRECISE | svtyper -i - -B $Tbam,$Nbam | /share/public/software/pigz/2.3.4/pigz -p 4 > $sample.gt.vcf.gz

rm *.orig.bam
