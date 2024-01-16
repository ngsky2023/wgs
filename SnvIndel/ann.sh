vcf=$1
sample=$2
/share/public/software/perl/5.24.1/bin/perl /share/public/database/Annovar/2017Jul16/convert2annovar.pl --includeinfo --format vcf4 --outfile $sample.txt $vcf &&
/share/public/software/perl/5.24.1/bin/perl /share/public/pipeline/RD_Onc/WGS/V0.9.1/table_annovar.pl --buildver hg19 --thread 2 --remove --otherinfo --protocol refGene,wgEncodeGencodeBasicV19,cytoBand,avsnp150,clinvar_20170130,cosmic70,popfreq_max_20150413,genomicSuperDups -operation g,g,r,f,f,f,f,r -nastring . $sample.txt /share/public/database/Annovar/2017Jul16/humandb/ --outfile $sample\_ann
/share/public/software/perl/5.24.1/bin/perl /share/public/pipeline/RD_Onc/WGS/V0.9.1/add_hgvs_0216.pl /share/work1/wangrr/DB/hg19/dbsnp141.total.hgvs $sample\_ann.hg19_multianno.txt $sample\_ann.hg19_multianno.txt.new &&
mv $sample\_ann.hg19_multianno.txt.new $sample\_ann.hg19_multianno.txt
/share/public/software/perl/5.24.1/bin/perl /share/public/pipeline/RD_Onc/WGS/V0.9.1/anno_table.Tp.pl -i $sample\_ann.hg19_multianno.txt -o $sample.xls -n $sample -s MuTect2 -a /share/work1/wangrr/DB/hg19/9606.geneid2go_kegg.xls
