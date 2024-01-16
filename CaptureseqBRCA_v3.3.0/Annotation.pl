#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';
use Cwd;

my $pwd = getcwd;
my $usage =
"
Usage:

  Options:
  -i|input      <FILE>   Input file of SNP,formate VCF.
  -g|gene_info  <FILE>   gene_info for annotation, default is '/home/leiyoubing/work/software/Gene/gene_info'.
  -bed          <FILE>   The bed-file of target region.
  -s|sample     <STR>    Sample name.
  -t|type       <STR>    variant type, 'snp' or 'indel'.
  -id|taxid     <NUM>    Tax id, such as mouse is 10090, default is human '9606'.
  -w|soft       <STR>    Software of call SNV, GATK or samtools? default is 'GATK'.
  -anno         <STR>    The all path of annovar, default ia '/home/zhangyu/workdir/exon_seq/Annotate/annovar'.
  -o|outdir     <DIR>    Output file direcory, if this is not exists created, default(./).
  -v|hgvs	<DIR>    HGVS format database
  -h|help                Help

For example:
        perl $0 -i sample.snp.vcf -s sample -t snp -o ./ -bed S04380110_Regions.bed
";

my ($input,$gene_info,$bed,$sample,$type,$taxid,$soft,$annovar,$outdir,$gene2go,$hgvs,$interpretation,$help);
GetOptions(
  "i|input=s"=>\$input,
  "g|gene_info=s"=>\$gene_info,
  "bed=s"=>\$bed,
  "s|sample=s"=>\$sample,
  "interpretation=s"=>\$interpretation,
  "t|type=s"=>\$type,
  "id|taxid=s"=>\$taxid,
  "w|soft=s"=>\$soft,
  "anno=s"=>\$annovar,
  "o|outdir=s"=>\$outdir,
  "v|hgvs=s"=>\$hgvs,
  "h|help"=>\$help
);
$outdir ||= $pwd;
$taxid ||= 9606;
$soft ||= 'GATK';
$annovar ||= "/home/zhangyu/workdir/exon_seq/Annotate/annovar/2016Feb01/annovar";
$gene_info ||= '/share/work2/staff/jiangdezhi/DataBase/GO_KEGG_Ann/Annotation/9606.geneid2go_kegg.xls';
$hgvs ||= "/home/zhangyu/workdir/database/dbSNP/download/dbsnp141.total.hgvs";
if (!$input || !$sample || !$type || !$bed || $help){
        die "$usage\n";
}

system("mkdir -m 755 -p $outdir/$sample/4.Annot") unless (-e "$outdir/$sample/4.Annot");

open SH,">$outdir/$sample/4.Annot/annotate.$type.$sample.sh" || die "$!";
print SH "#!/bin/bash\n";
print SH "echo ===$type annotate start at : `date` ===\n";
if($interpretation !~ /Gynecological/){
	print SH "perl $annovar/convert2annovar.pl --includeinfo --format vcf4 --allsample --outfile $outdir/$sample/4.Annot/$sample\_$type.annovar  $input\n\n";  
	print SH "perl $annovar/table_annovar.pl --buildver hg19 --thread 3 --remove --otherinfo --protocol refGene,wgEncodeGencodeBasicV19,mitimpact24,cytoBand,avsnp147,clinvar_20170130,cosmic70,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,1000g2015aug_amr,esp6500siv2_all,esp6500siv2_ea,esp6500siv2_aa,exac03,popfreq_max_20150413,cadd,dann,gerp++gt2,dbnsfp33a,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,wgEncodeBroadHmmGm12878HMM,tfbsConsSites,phastConsElements46way,phastConsElements100way,spidex -operation g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,r,f -nastring . $outdir/$sample/4.Annot/$sample\_$type.annovar.$sample.avinput $annovar/humandb --outfile $outdir/$sample/4.Annot/$sample\_$type.annovar\n";
	print SH "perl $Bin/add_hgvs_0216.pl $hgvs $outdir/$sample/4.Annot/$sample\_$type.annovar.hg19_multianno.txt $outdir/$sample/4.Annot/$sample\_$type.annovar.hg19_multianno.txt.new\n";
	print SH "mv $outdir/$sample/4.Annot/$sample\_$type.annovar.hg19_multianno.txt.new $outdir/$sample/4.Annot/$sample\_$type.annovar.hg19_multianno.txt\n";
	#print SH "perl $Bin/annovar_statistic.pl -a $outdir/$sample/Annot/$sample\_$type.annovar.hg19_multianno.txt -s $outdir/$sample/Annot/$sample -t $type\n\n";
	print SH "perl $Bin/anno_table.Tp.pl -i $outdir/$sample/4.Annot/$sample\_$type.annovar.hg19_multianno.txt -o $outdir/$sample/4.Annot/$sample\_$type.annovar.xls -n $sample  -s Mixture -a $gene_info -bed $bed\n\n";
	print SH "echo \"$sample $type annotate finish\" >$outdir/$sample/4.Annot/annotate_$type.$sample.mark\n";
}else{
	 print SH "echo \"$sample $type annotate finish\" >$outdir/$sample/4.Annot/annotate_$type.$sample.mark\n";
}
print SH "echo ===$type annotate end at : `date` ===\n";
close SH;
