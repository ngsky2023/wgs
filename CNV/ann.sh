for ss in *.cnv.vcf
do
sample=$(basename $ss | sed 's/.cnv.vcf//')
echo "/share/public/software/perl/5.24.1/bin/perl /share/public/database/Annovar/2017Jul16/convert2annovar.pl --includeinfo --format vcf4 -allsample --outfile $sample $sample.cnv.vcf
/share/public/software/perl/5.24.1/bin/perl /share/public/database/Annovar/2017Jul16/table_annovar.pl -buildver hg19 --remove -otherinfo --protocol refGene,wgEncodeGencodeBasicV19,cytoBand,genomicSuperDups,wgRna,targetScanS,dgvMerged,gwasCatalog,wgEncodeBroadHmmGm12878HMM,tfbsConsSites -operation g,g,r,r,r,r,r,r,r,r -nastring . $sample.TUMOR.avinput /share/public/database/Annovar/2017Jul16/humandb/ --outfile $sample
/share/public/pipeline/RD_Onc/WGS/V0.9/SCNAscore_rank_stat.pl -i $sample.hg19_multianno.txt > $sample.score"
done
