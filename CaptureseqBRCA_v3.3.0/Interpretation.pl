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
  -s|sample     <STR>    Sample name.
  -t|type       <STR>    variant type, 'snp' or 'indel'.
  -k|Xknown	<DIR>	 known hot spot for known disease 
  -d|Xindel	<DIR>    known INDEL for known disease 
  -o|outdir     <DIR>    Output file direcory, if this is not exists created, default(./).
  -h|help                Help

For example:
        perl $0 -i \$input -s sample_name  -o \$outdir -d \$Xindel -k \$Xknown  
";

my ($input,$gene_info,$bed,$sample,$type,$Xknown,$Xindel,$outdir,$gene2go,$help);
my ($Pdepth,$Pratio,$Vdepth,$Vratio,$Bdepth,$Bratio,);
my ($iPdepth,$iPratio,$iVdepth,$iVratio,$iBdepth,$iBratio,);
GetOptions(
        "i|input=s"=>\$input,
        "s|sample=s"=>\$sample,
        "t|type=s"=>\$type,
        "k|Xknown=s"=>\$Xknown,
        "d|Xindel=s"=>\$Xindel,
        "o|outdir=s"=>\$outdir,
        "h|help"=>\$help,
        "Pdp:i"=>\$Pdepth,
	"Pr:i"=>\$Pratio,
	"Vdp:i"=>\$Vdepth,
	"Vr:i"=>\$Vratio,
	"Bdp:i"=>\$Bdepth,
	"Br:i"=>\$Bratio,
        "iPdp:i"=>\$iPdepth,
	"iPr:i"=>\$iPratio,
	"iVdp:i"=>\$iVdepth,
	"iVr:i"=>\$iVratio,
	"iBdp:i"=>\$iBdepth,
	"iBr:i"=>\$iBratio,
);
$outdir ||= $pwd;
#$Xknown ||= "/home/zhangyu/workdir/exon_seq/Annotate/genetic_must/x-linkage/final/x-link_v1.sum.txt";
my $clinvar_P = "/share/public/database/Gynecological_cancer_backup/IPDB/ver4/DB_clinvar20180326_GRCh37_P.txt";
$Xindel ||= "/home/zhangyu/workdir/exon_seq/Annotate/genetic_must/x-linkage/final/x-link_v1.indel_transfor.txt";
$type ||="Xlinked";

if (!$input || !$sample||!$type  || $help){
        die "$usage\n";
}
my $panel_n;
if ($outdir =~ /\/.*\_(\d+)genes.*\/analysis/){
        $panel_n = $1;
}else{
        $panel_n = 2;
}
if($type eq "Xlinked"){
	system("mkdir -p $outdir/$sample/5.Interpretation/$type") unless (-e "$outdir/$sample/5.Interpretation/$type");
	open SH,">$outdir/$sample/5.Interpretation/$type/Interpretation_$sample.sh" || die "$!";
	print SH "#!/bin/bash\n";
	print SH "echo ===$sample interpretation start at : `date` ===\n";
	print SH "cat $input/$sample\_snp.annovar.xls $input/$sample\_indel.annovar.xls > $outdir/$sample/5.Interpretation/$type/$sample.annovar.xls\n";
	print SH "perl $Bin/SearchMust.pl $outdir/$sample/5.Interpretation/$type/$sample.annovar.xls $Xknown $outdir/$sample/5.Interpretation/$type/report.out\n";
	print SH "perl $Bin/SearchIndel.pl $outdir/$sample/5.Interpretation/$type/$sample.annovar.xls $Xindel $Xknown $outdir/$sample/5.Interpretation/$type/report.indel\n";
	print SH "if \[ \-f \"$outdir/$sample/5.Interpretation/$type/report.indel.info\" \]\;\n";
	print SH "then\n";
	print SH "cat $outdir/$sample/5.Interpretation/$type/report.out $outdir/$sample/5.Interpretation/$type/report.indel.info \>$outdir/$sample/5.Interpretation/$type/report\n";
	print SH "else\n";
	print SH "mv $outdir/$sample/5.Interpretation/$type/report.out $outdir/$sample/5.Interpretation/$type/report\n";
	print SH "fi\n";
	print SH "echo \"$sample interpretation finish\" >$outdir/$sample/5.Interpretation/$type/Interpretation_$sample.mark\n";
	print SH "echo ===$sample interpretation end at : `date` ===\n";
	close SH;

}elsif($type eq "ACMG"){
	system("mkdir -p $outdir/$sample/5.Interpretation/$type") unless (-e "$outdir/$sample/5.Interpretation/$type");
        open SH,">$outdir/$sample/5.Interpretation/$type/Interpretation_$sample.sh" || die "$!";
        print SH "#!/bin/bash\n";
        print SH "echo ===$sample interpretation start at : `date` ===\n";
        print SH "cat $input/$sample\_snp.annovar.xls $input/$sample\_indel.annovar.xls > $outdir/$sample/5.Interpretation/$type/$sample.annovar.xls\n";
        print SH "perl $Bin/Clinical_Level_single.pl -i $outdir/$sample/5.Interpretation/$type/$sample.annovar.xls -o $outdir/$sample/5.Interpretation/$type/$sample\_Clinical.xls\n";
        print  SH "perl $Bin/ACMG_level.pl $outdir/$sample/5.Interpretation/$type/$sample\_Clinical.xls  $outdir/$sample/5.Interpretation/$type/$sample\_level.xls\n";
	print SH "perl $Bin/for_Report_single.pl $outdir/$sample/5.Interpretation/$type/$sample\_level.xls  $outdir/$sample/5.Interpretation/$type/$sample\_summary.xls $outdir/$sample/5.Interpretation/$type/$sample\_VUS.xls\n";
	print SH "echo ===$sample interpretation end at : `date` ===\n";
}elsif($type eq "ACMG_Gyne"){
        system("mkdir -m 755 -p $outdir/$sample/5.Interpretation/$type") unless (-e "$outdir/$sample/5.Interpretation/$type");
        open SH,">$outdir/$sample/5.Interpretation/$type/Interpretation_$sample.sh" || die "$!";
        print SH "#!/bin/bash\n";
        print SH "echo ===$sample interpretation start at : `date` ===\n";
        #print SH "cat $input/$sample\_snp.annovar.xls $input/$sample\_indel.annovar.xls > $outdir/$sample/5.Interpretation/$type/$sample.annovar.xls\n";
        print SH "transvar ganno --refseq  --vcf $outdir/$sample/3.Variants/snp/$sample.snp.vcf > $outdir/$sample/3.Variants/snp/$sample.snp.transvar.vcf\n";
        print SH "transvar ganno --refseq  --vcf $outdir/$sample/3.Variants/indel/$sample.indel.vcf > $outdir/$sample/3.Variants/indel/$sample.indel.transvar.vcf\n";
        print SH "perl $Bin/Annotation_transvar_v2.pl -i $outdir/$sample/4.Annot/$sample\_snp.annovar.xls  -tr $outdir/$sample/3.Variants/snp/$sample.snp.transvar.vcf -o $outdir/ -s $sample -log $outdir/$sample/5.Interpretation/$type/$sample\_snp_transvar.log\n";
        print SH "perl $Bin/Annotation_transvar_v2.pl -i $outdir/$sample/4.Annot/$sample\_indel.annovar.xls  -tr $outdir/$sample/3.Variants/indel/$sample.indel.transvar.vcf -o $outdir/ -s $sample -log $outdir/$sample/5.Interpretation/$type/$sample\_indel_transvar.log\n";
        #print SH "perl $Bin/Annotation_transvar.pl -o $outdir/ -s $sample -log $outdir/$sample/5.Interpretation/$type/$sample\_transvar.log\n";
        print SH "perl $Bin/annovar.3mode.pl  $outdir/$sample/5.Interpretation/$type/$sample\_snp.annovar.3mode.tmp.xls  $outdir/$sample/5.Interpretation/$type/$sample\_snp.annovar.3mode.xls\n";
        print SH "perl $Bin/annovar.3mode.pl  $outdir/$sample/5.Interpretation/$type/$sample\_indel.annovar.3mode.tmp.xls  $outdir/$sample/5.Interpretation/$type/$sample\_indel.annovar.3mode.xls\n";
        print SH "rm  $outdir/$sample/5.Interpretation/$type/$sample\_snp.annovar.3mode.tmp.xls && rm  $outdir/$sample/5.Interpretation/$type/$sample\_indel.annovar.3mode.tmp.xls\n";
        #print SH "cat $outdir/$sample/4.Annot/$sample\_snp.annovar.xls $outdir/$sample/5.Interpretation/$type/$sample\_indel.annovar.3mode.xls > $outdir/$sample/5.Interpretation/$type/$sample\_annovar.xls\n";
        print SH "perl $Bin/annovar_merge.pl -f1 $outdir/$sample/5.Interpretation/$type/$sample\_snp.annovar.3mode.xls -f2 $outdir/$sample/5.Interpretation/$type/$sample\_indel.annovar.3mode.xls -o $outdir/$sample/5.Interpretation/$type/$sample\_annovar.xls\n";
        print SH "perl $Bin/Clinical_Level_single.pl -i $outdir/$sample/5.Interpretation/$type/$sample\_annovar.xls -o $outdir/$sample/5.Interpretation/$type/$sample\_Clinical.xls\n";
        print  SH "perl $Bin/ACMG_level.pl $outdir/$sample/5.Interpretation/$type/$sample\_Clinical.xls  $outdir/$sample/5.Interpretation/$type/$sample\_level.xls\n";
        print SH "perl $Bin/for_Report_single.pl $outdir/$sample/5.Interpretation/$type/$sample\_level.xls  $outdir/$sample/5.Interpretation/$type/$sample\_P.xls $outdir/$sample/5.Interpretation/$type/$sample\_VUS.xls\n";
        #print SH "perl $Bin/for_Report_screen.pl $outdir/$sample/5.Interpretation/$type/$sample\_level.xls  $outdir/$sample/5.Interpretation/$type/$sample.level1 $outdir/$sample/5.Interpretation/$type/$sample.level2 $outdir/$sample/5.Interpretation/$type/$sample.level3\n";
        print SH "grep -v \"\#\" $outdir/$sample/3.Variants/snp/$sample.snp.vcf >$outdir/$sample/5.Interpretation/$type/$sample.snp\n";
        print SH "grep -v \"\#\" $outdir/$sample/3.Variants/indel/$sample.indel.vcf >$outdir/$sample/5.Interpretation/$type/$sample.indel\n";
        print SH "cat $outdir/$sample/5.Interpretation/$type/$sample.snp $outdir/$sample/5.Interpretation/$type/$sample.indel >$outdir/$sample/5.Interpretation/$type/$sample.snp_indel.txt\n";
        print SH "perl $Bin/SearchMust.Gynecological_tumor.vcf.pl $outdir/$sample/5.Interpretation/$type/$sample.snp_indel.txt  $clinvar_P  $outdir/$sample/5.Interpretation/$type/$sample.vcf.report.xls\n";
        #print SH "perl $Bin/S1explain.pl -i $outdir/$sample/5.Interpretation/$type/$sample.vcf.report.xls\n";
        #print SH "perl $Bin/S1explain_verscreen.pl -i $outdir/$sample/5.Interpretation/$type/$sample.level1 -o $outdir/$sample/5.Interpretation/$type/$sample.level1.info\n";
        #print SH "perl $Bin/S1explain_verscreen.pl -i $outdir/$sample/5.Interpretation/$type/$sample.level2 -o $outdir/$sample/5.Interpretation/$type/$sample.level2.info\n";
        #print SH "perl $Bin/S1explain_verscreen2.pl -i $outdir/$sample/5.Interpretation/$type/$sample.level3 -o $outdir/$sample/5.Interpretation/$type/$sample.level3.info\n";
        #print SH "cat $input/$sample\_snp.annovar.$sample.avinput $input/$sample\_indel.annovar.$sample.avinput > $outdir/$sample/5.Interpretation/$type/$sample.annovar.avinput\n";
        #print SH "perl $Bin/SearchMust.Gynecological_tumor.tmp.pl $outdir/$sample/5.Interpretation/$type/$sample.annovar.avinput  $Xknown  $outdir/$sample/5.Interpretation/$type/$sample\_report.xls\n";
        print SH "perl $Bin/level_1base.pl -f $outdir/$sample/5.Interpretation/$type/$sample\_level.xls -v $outdir/$sample/4.Annot/$sample\_indel.annovar.$sample.avinput -o $outdir/$sample/5.Interpretation/$type/$sample\_level_1base.xls\n";
        #print SH "perl $Bin/result_screen_v2.pl $outdir/$sample/5.Interpretation/$type/$sample\_level.xls $outdir/$sample/4.Annot/$sample\_indel.annovar.BRCA_$sample.avinput  $outdir/$sample/3.Variants/snp/$sample.atlas.snp.vcf $outdir/$sample/3.Variants/indel/$sample.atlas.indel.vcf  $sample\_result_v2.xls  $panel_n\n";
        print SH "perl $Bin/result_screen_v2.pl -l $outdir/$sample/5.Interpretation/$type/$sample\_level.xls -a $outdir/$sample/4.Annot/$sample\_indel.annovar.$sample.avinput  -as $outdir/$sample/3.Variants/snp/$sample.atlas.snp.vcf -ai $outdir/$sample/3.Variants/indel/$sample.atlas.indel.vcf  -o $sample\_result_v2.xls -Pdp $Pdepth -Pr $Pratio -Vdp $Vdepth -Vr $Vratio -Bdp $Bdepth -Br $Bratio -n $panel_n\n";
        #print SH "perl $Bin/result_info_v2.pl -i $sample\_result_v2.xls -o $sample.result.info -n $panel_n\n";
        print SH "perl $Bin/result_info_v2.pl -i $sample\_result_v2.xls -o $sample.result.info -iPdp $iPdepth -iPr $iPratio -iVdp $iVdepth -iVr $iVratio -iBdp $iBdepth -iBr $iBratio -n $panel_n\n";
	print SH "echo \"$sample interpretation finish\" >$outdir/$sample/5.Interpretation/$type/Interpretation_$sample.mark\n";   #Z17T00142_indel.annovar.BRCA_Z17T00142.avinput
	print SH "echo ===$sample interpretation end at : `date` ===\n";
        close SH;
}else{
	system("mkdir -p $outdir/$sample/5.Interpretation/$type") unless (-e "$outdir/$sample/5.Interpretation/$type");

        open SH,">$outdir/$sample/5.Interpretation/$type/Interpretation_$sample.sh" || die "$!";
        print SH "#!/bin/bash\n";
        print SH "echo ===$sample interpretation start at : `date` ===\n";
        print SH "grep -v \"\#\" $outdir/$sample/3.Variants/snp/$sample.snp.vcf >$outdir/$sample/5.Interpretation/$type/$sample.snp\n";
	print SH "grep -v \"\#\" $outdir/$sample/3.Variants/indel/$sample.indel.vcf >$outdir/$sample/5.Interpretation/$type/$sample.indel\n";
	print SH "cat $outdir/$sample/5.Interpretation/$type/$sample.snp $outdir/$sample/5.Interpretation/$type/$sample.indel >$outdir/$sample/5.Interpretation/$type/$sample.snp_indel.txt\n";
	print SH "perl $Bin/SearchMust.Gynecological_tumor.vcf.pl $outdir/$sample/5.Interpretation/$type/$sample.snp_indel.txt  $Xknown  $outdir/$sample/5.Interpretation/$type/$sample.vcf.report.xls\n";

	#print SH "cat $input/$sample\_snp.annovar.$sample.avinput $input/$sample\_indel.annovar.$sample.avinput > $outdir/$sample/5.Interpretation/$type/$sample.annovar.avinput\n";
        #print SH "perl $Bin/SearchMust.Gynecological_tumor.tmp.pl $outdir/$sample/5.Interpretation/$type/$sample.annovar.avinput  $Xknown  $outdir/$sample/5.Interpretation/$type/$sample\_report.xls\n";
        print SH "echo \"$sample interpretation finish\" >$outdir/$sample/5.Interpretation/$type/Interpretation_$sample.mark\n";
        print SH "echo ===$sample interpretation end at : `date` ===\n";
        close SH;
}
