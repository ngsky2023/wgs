sample=$1
dbsnpHGVS=$2
annovarDB=$3
annovar=$4
strelka=$5
genome=$6

bin=`dirname $0`
realpath=`pwd`/$1
outdir=`dirname $realpath`
samples=`cat $1|cut -f 2`

if [ ! -e $outdir/strelka ] ;then
    mkdir -p $outdir/strelka
fi 

cd $outdir/strelka

for sample in $samples
do 
    bam=$outdir/analysis/$sample/2.Realign/$sample.bam
    $strelka --bam=$bam --referenceFasta=$genome --runDir=$sample  &&
    cd $sample &&
    ./runWorkflow.py -m local -j 20
    perl $bin/filter_germ.pl -i results/variants/variants.vcf.gz > variants.vcf 2> filter.log
    /share/public/software/perl/5.24.1/bin/perl $annovar/convert2annovar.pl --includeinfo --format vcf4 --allsample --outfile variants.annovar variants.vcf
    /share/public/software/perl/5.24.1/bin/perl $annovar/table_annovar.pl --buildver hg19 --thread 20 --remove --otherinfo -nastring . --protocol refGene,wgEncodeGencodeBasicV19,HGMD,cytoBand,avsnp150,clinvar,cosmic91,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,1000g2015aug_amr,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,exac03,popfreq_max_20150413,gnomad211_exome,Chinese_freq,cadd,dann,gerp++gt2,dbnsfp35a,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,wgEncodeBroadHmmGm12878HMM,tfbsConsSites,phastConsElements46way,phastConsElements100way,spidex,dbscsnv11,icgc21,nci60 -operation g,g,f,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,r,f,f,f,f --outfile variants variants.annovar.$sample.avinput $annovarDB
    
    /share/public/software/perl/5.24.1/bin/perl $bin/../SnvIndel/add_hgvs_0216.pl $dbsnpHGVS variants.hg19_multianno.txt variants.hg19_multianno.txt.new && mv variants.hg19_multianno.txt.new variants.hg19_multianno.txt
done
