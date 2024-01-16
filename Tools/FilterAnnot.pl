#!/usr/bin/env perl
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
use FindBin '$Bin';

my $usage =
"
Usage:
	 
Options:
	-input    <FILE>  The input sample list (same with the one used at the input of ExomeSomatic pipeline), required.
	-inputdir <DIR>   The directory which contains outputs from ExomeSomatic pipeline, required.
	-outdir   <DIR>   Output file direcory, if this is not exists created, required.
	-tid      <FILE>  Major transcript ID information, default is '/share/work1/staff/zhangyu/database/clinvar/20160307/variant_summary.GRCh37.gene_NM.txt'.
	-genesum  <FILE>  NCBI gene summary information, default is '/share/work1/staff/lijianbai/HumanRefs/NCBI_Gene/gene_RefSeqGene_Summary.txt'.
	-clinvar  <FILE>  ClinVar records for mutations, default is '/share/work1/staff/zhangyu/database/clinvar/20160307/clinvar_20160302.vcf'.
	-geneprob <FILE>  The file contains TUSON probability for each gene to be drivers, and ranks accordingly, default is '/share/work1/staff/lijianbai/HumanRefs/TUSON_Lists/Probability_TUSON.txt'.
	-popfreq  <NUM>   The maximum basic population frequency allowed for it's being called as somatic, default is 0.05 for somatic mode, and 0.001 for putative somatic mode.
	-mode     <STR>   The filtration mode, either 'somatic' or 'putative_somatic', required. 
	-pname    <STR>   Project name, which will be used as the prefix of output files, required.
	-var_freq <NUM>   Variant allele frequency, default is 0.01. 
	-totalD   <NUM>   Minimum total depth allowed to keep the variant, default is 10. 
	-varD     <NUM>   Minimum variant allele depth allowed to keep the variant, default is 4.
	-rand     <NUM>   Number of randomized reads or not randomized. Default is 'NO'. 

For example:
perl $0 -input InputFile -outdir OutputDIR -tid MajorTranscriptID_GD -genesum NCBI_GeneSummary -clinvar ClinVarRecords.vcf -popfreq MaxPopFreq -mode FiltrationMode -pname ProjectName -var_freq 0.01 -totalD 10 -varD 4
CAUTION: for multi tumor samples from single patient without normal control, run under putative somatic mode seperately for each tumor. 
";

my ($input,$inputdir,$outdir,$tid,$genesum,$clinvar,$geneprob,$popfreq,$mode,$var_freq,$totalD,$varD,$rand);
GetOptions(
	"input:s"=>\$input,
	"inputdir:s"=>\$inputdir,
	"outdir:s"=>\$outdir,
	"tid:s"=>\$tid,
	"genesum:s"=>\$genesum,
	"clinvar:s"=>\$clinvar,
	"geneprob:s"=>\$geneprob,
	"popfreq:s"=>\$popfreq,
	"mode:s"=>\$mode,
	"pname:s"=>\$pname,
	"var_freq:s"=>\$var_freq,
	"totalD:s"=>\$totalD,
	"varD:s"=>\$varD,
	"rand:s"=>\$rand
);
chomp $pname;

$tid ||="/share/public/database/clinvar/current/variant_summary.GRCh37.gene_NM.txt";
$genesum ||="/share/work2/baijian/HumanRefs/gene_RefSeqGene_Summary.txt";
$clinvar ||="/share/public/database/clinvar/current/clinvar.vcf";
$geneprob ||="/share/work2/baijian/HumanRefs/TUSON_Lists/Probability_TUSON.txt";
$var_freq ||=0.01;
$totalD ||=10;
$varD ||=4;
$rand ||="NO";

if($mode eq "somatic"){
	$popfreq ||=0.05;
}else{
	$popfreq ||=0.001;
}

die "$usage\n" if(!$input || !$inputdir || !$outdir || !$mode);
#################################################################################################################################
open(FR,"$input")||die "$!";
$i=0;
while(defined($a=<FR>)){
	chomp $a;
	@b=split(/\t/,$a);
	$samp{$b[0]} = 1;
}
close(FR);
@pairs=sort keys %samp;
$i=0;$a="";@b=();$n=$t=$pair="";

%samp=();
###########################################################################################################################################
open(GP,"$geneprob")||die "$!";
while(defined($gp=<GP>)){
	chomp $gp;
	@gp=split(/\t/,$gp);
	next if $gp[0] eq "Gene";
	$ogscore{$gp[0]}=$gp[11];
	$ogr{$gp[0]}=$gp[12];
	$tsgscore{$gp[0]}=$gp[5];
	$tsgr{$gp[0]}=$gp[6];
}
close(GP);
$gp="";@gp=();
## Gene summary from NCBI
open(SUM,"$genesum")||die "$!";
while(defined($sum=<SUM>)){
	chomp $sum;
	@sum=split(/\t/,$sum);
	$sum{$sum[2]}=$sum[4];
}
close(SUM);
$sum="";@sum=();
###########################################################################################################################################
$raw_var="$pname"."_Somatic_Raw.xls";
$j=0;
open(FW,">$outdir/$raw_var")||die "$!";
for my $sample (@pairs){
	my $pp = "$inputdir/snv-indel/annotation/$sample/chr*/$sample.chr*.xls";
	for my $file (<$pp>){
		open(MUT,$file)||die "$!";
		while(defined($a=<MUT>)){
			chomp $a;
			if($j==0 and $a=~/SampleName/){
				print FW "$a\n";
			}elsif(not $a=~/SampleName/){
				print FW "$a\n";
			}
		}
		close(MUT);$a="";
		$j++;
	}
}
close(FW);$j=0;

open(FR,"$outdir/$raw_var")||die "$!";
while(defined($raw=<FR>)){
	chomp $raw;
	@b=split(/\t/,$raw);
	next if ($b[0] eq "SampleName");
	next if ($b[43]>=$popfreq and $ogr{$b[12]}>200 and $tsgr{$b[12]}>200);
	next if ($b[43]>10*$popfreq and ($ogr{$b[12]}<=200 or $tsgr{$b[12]}<=200));
	$locus="$b[1]\t$b[2]\t$b[3]";
	$samplocus="$b[0]\t$b[1]\t$b[2]\t$b[3]";
	$exist{$locus}=1;
	if(not $exist{$samplocus}){
		$count{$locus}++;
		$exist{$samplocus}=1;
	}
	$countsamp{$b[0]}=1;
}
close(FR);
$raw=$locus=$samplocus="";@b=();
@countsamp=keys %countsamp;
#################################################################################################################################################
# Filter 1: normal coverage >= 8, normal mut read <=1(or <0.005), tumor coverage >=10, tumor mut reads >=2; loci with val allele contain NonRef call
# Meanwhile, output regions of interest (ROI) bed files and SeattleInput files
# while($j<@pairs){
# 	open(SNV,"$inputdir/somatic_$pairs[$j]/2.merged/snp/somatic.merged1.vcf")||die "$!";
# 	open(INDEL,"$inputdir/somatic_$pairs[$j]/2.merged/indel/somatic.merged1.vcf")||die "$!";
# 	while(defined($a=<SNV>)){
# 		next if ($a=~/^\#/);
# 		chomp $a;
# 		@b=@n=@t=();
# 		@b=split(/\t/,$a);
# 		@n=split(/:/,$b[9]);
# 		@t=split(/:/,$b[10]);
# 		$rdn=$n[3]+$n[4];
# 		$rdt=$t[3]+$t[4];
# 		$vafn=$n[4]/$rdn;
# 		$vaft=$t[4]/$rdt;
# 		next if ($rdn<8 or $rdt<10 or ($n[4]>2 and $vafn>=0.01) or $t[4]<=2);
# 		next if ($b[4]=~/NON_REF/);
# 		next if ($b[4]=~/,/);
# 		$loci="somatic\_$pairs[$j]\t$b[0]\t$b[1]";
# 		$snv{$loci}=1;
# 	}
# 	close(SNV);
# 	$a="";@b=@n=@t=();$rdn=$rdt=$vafn=$vaft=$start=$end=0;

# 	while(defined($a=<INDEL>)){
# 		next if ($a=~/^\#/);
# 		chomp $a;
# 		@b=@n=@t=();
# 		@b=split(/\t/,$a);
# 		@n=split(/:/,$b[9]);
# 		@t=split(/:/,$b[10]);
# 		$rdn=$n[3]+$n[4];
# 		$rdt=$t[3]+$t[4];
# 		$vafn=$n[4]/$rdn;
# 		$vaft=$t[4]/$rdt;
# 		next if ($rdn<8 or $rdt<10 or ($n[4]>2 and $vafn>=0.01) or $t[4]<=2);
# 		next if ($b[4]=~/NON_REF/);
# 		next if ($b[4]=~/,/);
# 		$loci="somatic\_$pairs[$j]\t$b[0]\t$b[1]";
# 		$indel{$loci}=1;
# 	}
# 	close(INDEL);
# 	$a="";@b=@n=@t=();$rdn=$rdt=$vafn=$vaft=$start=$end=0;
# 	## then filter the specific loci in varscan raw snp and indel outputs
# 	open(SNV_VAR,"$inputdir/somatic_$pairs[$j]/0.variant/VarScan/varscan_snp.raw.vcf")||die "$!";
# 	open(INDEL_VAR,"$inputdir/somatic_$pairs[$j]/0.variant/VarScan/varscan_indel.raw.vcf")||die "$!";
# 	while(defined($a=<SNV_VAR>)){
# 		next if ($a=~/^\#/);
# 		chomp $a;
# 		@b=@n=@t=();
# 		@b=split(/\t/,$a);
# 		@n=split(/:/,$b[9]);
# 		@t=split(/:/,$b[10]);
# 		$rdn=$n[3]+$n[4];
# 		$rdt=$t[3]+$t[4];
# 		$vafn=$n[4]/$rdn;
# 		$vaft=$t[4]/$rdt;
# 		$loci="somatic\_$pairs[$j]\t$b[0]\t$b[1]";
# 		next if (not $snv{$loci});
# 		if($snv{$loci} and ($rdn<8 or $rdt<10 or ($n[4]>2 and $vafn>=0.01) or $t[4]<=2)){
# 			$snv{$loci}="";
# 		}
# 	}
# 	close(SNV_VAR);
# 	$a=$loci="";@b=@n=@t=();$rdn=$rdt=$vafn=$vaft=$start=$end=0;

# 	while(defined($a=<INDEL_VAR>)){
# 		next if ($a=~/^\#/);
# 		chomp $a;
# 		@b=@n=@t=();
# 		@b=split(/\t/,$a);
# 		@n=split(/:/,$b[9]);
# 		@t=split(/:/,$b[10]);
# 		$rdn=$n[3]+$n[4];
# 		$rdt=$t[3]+$t[4];
# 		$vafn=$n[4]/$rdn;
# 		$vaft=$t[4]/$rdt;
# 		$loci="somatic\_$pairs[$j]\t$b[0]\t$b[1]";
# 		next if (not $indel{$loci});
# 		if($indel{$loci} and ($rdn<8 or $rdt<10 or ($n[4]>2 and $vafn>=0.01) or $t[4]<=2)){
# 			$indel{$loci}="";
# 		}
# 	}
# 	close(INDEL_VAR);
# 	$a=$loci="";@b=@n=@t=();$rdn=$rdt=$vafn=$vaft=$start=$end=0; 
# 	$j++;
# }
# $j=0;

$filter_var="$pname"."_Somatic_Filter.xls";
open(ANNRAW,"$outdir/$raw_var")||die "$!";
open(ANN,">$outdir/$filter_var")||die "$!";
while(defined($a=<ANNRAW>)){
	chomp $a;
	@b=split(/\t/,$a);
	if($a=~/SampleName/){
		print ANN "$a\n";
	}else{
		$info="$b[1]\t$b[2]\t$b[4]\t$b[5]";
		$locus="$b[1]\t$b[2]\t$b[3]";
		next if not $exist{$locus};
		if($mode eq "somatic"){
			print ANN "$a\n";
		}elsif($mode eq "putative_somatic"){
			$freq=$count{$locus}/($#countsamp+1);
			next if ($count{$locus}>=2 and ($freq>=0.1 or $b[9]>=0.4) and $ogr{$b[12]}>500);
			next if $b[9]>=0.9;
			next if ($b[9]>=0.3 and ((($tsgr{$b[12]}>500 or (not $tsgr{$b[12]})) and ($ogr{$b[12]}>500 or (not $ogr{$b[12]}))) and $b[43]>0));
			next if ($b[9]>=0.3 and ((($tsgr{$b[12]}>500 or (not $tsgr{$b[12]})) and ($ogr{$b[12]}>500 or (not $ogr{$b[12]}))) and ($b[23 ne "."])));
			next if ($b[9]>=0.3 and (($b[49] eq "T") and ($b[51] eq "B") or (($b[49] eq "T") and ($b[53] eq "B")) or (($b[51] eq "B") and ($b[53] eq "B"))) and ($b[43]>0 or $b[44]<=1 or (($b[49] eq "T") and ($b[51] eq "B") and ($b[53] eq "B"))));
			# next if ($b[9]<=0.05 and ($tsgr{$b[12]}>500 or (not $tsgr{$b[12]})) and ($ogr{$b[12]}>500 or (not $ogr{$b[12]})) and $b[43]>=0 and (not $sum{$b[12]}=~/tumor/) and (not $sum{$b[12]}=~/cancer/) and (not $sum{$b[12]}=~/carcinoma/) and (not $sum{$b[12]}=~/sarcoma/) and (not $sum{$b[12]}=~/adenoma/) and (not $sum{$b[12]}=~/lymphoma/) and (not $sum{$b[12]}=~/leukemia/));
			print ANN "$a\n";
  		}
	}
}
close(ANNRAW);close(ANN);
$a="";@b=();$pos="";$loci=$info="";%snv=%indel=();
####################################################################################################################################################################################################################
## Transcript ID information
open(TI,"$tid")||die "$!";
while(defined($ti=<TI>)){
	chomp $ti;
	@ti=split(/\t/,$ti);
	$ti{$ti[0]}=$ti[1];
}
close(TI);
$ti="";@ti=();
## ClinVar records
open(CLV,"$clinvar")||die "$!";
$clnsig{0}="CLNSIG=Uncertain significance";
$clnsig{1}="CLNSIG=not provided";
$clnsig{2}="CLNSIG=Benign";
$clnsig{3}="CLNSIG=Likely benign";
$clnsig{4}="CLNSIG=Likely pathogenic";
$clnsig{5}="CLNSIG=Pathogenic";
$clnsig{6}="CLNSIG=drug response";
$clnsig{7}="CLNSIG=histocompatibility";
$clnsig{8}="CLNSIG=other";
while(defined($clv=<CLV>)){
	chomp $clv;
	next if $clv=~/^#/;
	@clv=split(/\t/,$clv);
	$refl=length $clv[3];
	if($clv[4]=~/,/){
		@var=split(/,/,$clv[4]);
		$varl=length $var[0];
	}else{
		$varl=length $clv[4];
	}
	if($refl==$varl){
		$pos=$clv[1];
		if(not $clv[4]=~/,/){
			$loci="chr$clv[0]\t$pos\t$clv[3]\t$clv[4]";
		}else{
			$loci1="chr$clv[0]\t$pos\t$clv[3]\t$var[0]";
			$loci2="chr$clv[0]\t$pos\t$clv[3]\t$var[1]";
		}
	}elsif($refl>$varl){
		$pos=$clv[1]+1;
		$loci="chr$clv[0]\t$pos\tDel";
	}elsif($refl<$varl){
		$pos=$clv[1];
		$loci="chr$clv[0]\t$pos\tIns";
	}
	@info=split(/\;/,$clv[7]);
	$i=0;
	while($i<@info){
		if($info[$i]=~/^CLNSIG/){
			@clnsig=split(/=/,$info[$i]);
			$clnsig=$clnsig{$clnsig[1]};
			$clnsig="CLNSIG=Uncertain significance" if not $clnsig{$clnsig[1]};
		}elsif($info[$i]=~/^CLNDBN/){
			$clndbn=$info[$i];
		}elsif($info[$i]=~/^CLNREVSTAT/){
			$clnrevstat=$info[$i];
		}elsif($info[$i]=~/^CLNACC/){
			$clnacc=$info[$i];
		}elsif($info[$i]=~/^CLNDSDB/){
			$clndsdb=$info[$i];
		}elsif($info[$i]=~/^CLNDSDBID/){
			$clndsdbid=$info[$i];
		}
		$i++;
	}
	$clv{$loci}="$clnsig;$clndbn;$clnrevstat;$clnacc;$clndsdb;$clndsdbid";
}
close(CLV);
$clv=$loci=$clnsig=$clndbn=$clnrevstat=$clnacc=$clndsdb=$clndsdbid="";@clv=@var=@info=();$i=0;$loci1=$loci2="";
#######################################################################################################################################
$rep_var="$pname"."_Somatic_Report.xls";
$repdetail_var="$pname"."_Somatic_ReportDetail.xls";
open(FILTER,"$outdir/$filter_var")||die "$!";
open(REPORT,">$outdir/$rep_var")||die "$!";
open(REPORTDET,">$outdir/$repdetail_var")||die "$!";
print REPORT "SampleID\tGeneSymbol\tMajorTranscriptID\tExon#\tcDNA_Change\tAA_Change\tRegion\tMutFunc\tReadDepth_Total\tReadDepth_Var\tVarAlleleFrac\tGeneDescrip\tGeneSummary_NCBI\tClinVar\tTSG_Score\tTSG_Rank\tOG_Score\tOG_Rank\n";
while(defined($a=<FILTER>)){
	chomp $a;
	@b=split(/\t/,$a);
	if($b[0] eq "SampleName"){
		print REPORTDET "$a\n";
	}
	next if ((not $b[11]=~/exonic/) and (not $b[11]=~/splicing/));
	next if (($b[14] eq "synonymous SNV") and (not $b[11]=~/splicing/));
	next if ($b[43]>$popfreq);
	next if (($b[8]/($b[7]+$b[8]))<$var_freq or ($b[7]+$b[8])<$totalD or $b[8]<$varD);
	print REPORTDET "$a\n";
	@transc=split(/,/,$b[15]);
	$rdtotal=$b[7]+$b[8];
	$rdvar=$b[8];
	$vaf=$rdvar/$rdtotal;
	$i=0;
	@info=();
	$control=0;
	while($i<@transc){
		@info=split(/:/,$transc[$i]);
		if(($ti{$info[0]} eq $info[1]) and $ti{$info[0]}){
			$ti=$info[1];
			$exon=$info[2];
			$cdna=$info[3];
			$aa=$info[4];
			$control=1;
		}
		$i++;
	}
	if($control==0 and (not $b[15]=~/^\./)){
		@info=split(/:/,$transc[0]);
		if($#info==4){
			$ti=$info[1];
			$exon=$info[2];
			$cdna=$info[3];
			$aa=$info[4];
		}elsif($#info==2 and $info[1]=~/^exon/){
			$ti=$info[0];
			$exon=$info[1];
			$cdna=$info[2];
			$aa="\.";
		}elsif($#info==2 and (not $info[1]=~/^exon/)){
			$ti="\.";$exon="\.";$cdna="\.";$aa="\.";
		}
	}elsif($control==0 and ($b[15]=~/^\./)){
		$ti="\.";$exon="\.";$cdna="\.";$aa="\.";
	}
	if(($b[4] ne "-") and ($b[5] ne "-")){
		$loci="$b[1]\t$b[2]\t$b[4]\t$b[5]";
		if($clv{$loci}){
			$clv=$clv{$loci};
		}else{
			$clv="NoRecord";
		}
	}elsif(($b[4] eq "-") and ($b[5] ne "-")){
		$loci="$b[1]\t$b[2]\tIns";
		if($clv{$loci}){
			$clv=$clv{$loci};
		}else{
			$clv="NoRecord";
		}
	}elsif(($b[4] ne "-") and ($b[5] eq "-")){
		$loci="$b[1]\t$b[2]\tDel";
		if($clv{$loci}){
			$clv=$clv{$loci};
		}else{
			$clv="NoRecord";
		}
	}
	if(not $sum{$b[12]}){
		$summary=".";
	}else{
		$summary=$sum{$b[12]};
	}
	print REPORT "$b[0]\t$b[12]\t$ti\t$exon\t$cdna\t$aa\t$b[11]\t$b[14]\t$rdtotal\t$rdvar\t$vaf\t$b[91]\t$summary\t$clv\t";
	if($tsgr{$b[12]}){
		print REPORT "$tsgscore{$b[12]}\t$tsgr{$b[12]}\t$ogscore{$b[12]}\t$ogr{$b[12]}\n";
	}else{
		print REPORT "NA\tNA\tNA\tNA\n";
	}
}
close(FILTER);close(REPORT);
$a=$ti=$exon=$cdna=$aa="";@b=();@transc=();$rdtotal=$rdvar=$vaf=$i=0;@info=();%sum=%clv=();
