#!/bin/bash
V=''
n=0
out=''
for vcf in $@
do
	n=$[$n+1]
	if [ $n -eq 1 ]
	then
		out=$vcf
	else
		V=" $V -V $vcf"
	fi
done
/share/public/software/java/1.8.0_162/bin/java -cp /share/public/software/gatk/3.8/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R /share/work1/wangrr/DB/hg19/hg19AddVirus.fa --assumeSorted $V -out $out
