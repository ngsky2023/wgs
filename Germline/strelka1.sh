
source /share/Onc_KYproject/hanwj4457/analysis/project/bashrc
bam=$1
sample=$2
#/share/public/software/strelka/2.8.4/bin/configureStrelkaGermlineWorkflow.py --bam=$bam --referenceFasta=/share/work1/wangrr/DB/hg19/hg19AddUnlocalizedVirus.fasta --runDir=$sample  &&
cd $sample &&
#./runWorkflow.py -m local -j 20
#perl /share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/strelka/pon/filter_germ.pl -i results/variants/variants.vcf.gz > variants.vcf 2> filter.log
/share/public/software/perl/5.24.1/bin/perl /share/Onc_Soft_DB/database/annovar/Annovar/convert2annovar.pl --includeinfo --format vcf4 --allsample --outfile variants.annovar variants.vcf

/share/public/software/perl/5.24.1/bin/perl /share/Onc_Soft_DB/database/annovar/Annovar/table_annovar.pl --buildver hg19 --thread 20 --remove --otherinfo -nastring . --protocol refGene,wgEncodeGencodeBasicV19,HGMD,cytoBand,avsnp150,clinvar,cosmic91,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,1000g2015aug_amr,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,exac03,popfreq_max_20150413,gnomad211_exome,Chinese_freq,cadd,dann,gerp++gt2,dbnsfp35a,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,wgEncodeBroadHmmGm12878HMM,tfbsConsSites,phastConsElements46way,phastConsElements100way,spidex,dbscsnv11,icgc21,nci60 -operation g,g,f,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,r,f,f,f,f --outfile variants variants.annovar.$sample.avinput /share/Onc_Soft_DB/database/annovar/humandb/

/share/public/software/perl/5.24.1/bin/perl /share/work1/yangrt/workdir/research/git_release3/GeneticTumor/Target_Pipeline/Annotation//add_hgvs_0216.pl /share/public/software/Onc_Soft/database/dbsnp/dbsnp141.total.hgvs variants.hg19_multianno.txt variants.hg19_multianno.txt.new && mv variants.hg19_multianno.txt.new variants.hg19_multianno.txt
