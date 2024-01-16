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

my (@input,$outdir,$tid,$genesum,$clinvar,$geneprob,$popfreq,$mode,$var_freq,$totalD,$varD,$rand);
GetOptions(
	"input=s{1,}"=>\@input,
	"outdir:s"=>\$outdir,
	"tid:s"=>\$tid,
	"historyDir:s"=>\$historyDir,
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

$tid ||="$Bin/gene_NM.list";
$genesum ||="/share/work1/wangrr/DB/hg19/gene_RefSeqGene_Summary.txt";
$clinvar ||="/share/public/database/clinvar/20170619/current/clinvar.vcf";
$geneprob ||="$Bin/Probability_TUSON.txt";
$var_freq ||=0.03;
$totalD ||=8;
$varD ||=4;
$rand ||="NO";

if($mode eq "somatic"){
	$popfreq ||=0.05;
}else{
	$popfreq ||=0.001;
}

die "$usage\n" if($#input < 0  || !$outdir);
######################################################################
my %hotspots;
open IN , "/share/work1/wangrr/DB/hg19/cancer_hotspots.tsv";
<IN>;
while (<IN>){
	chomp;
	my ($gene , $pos , $m) = (split /\t/ , $_)[0,1,2];
	my @m = split /\|/ , $m;
	my @mm;
	for my $m (@m){
		$m =~ s/:\d+//;
		$m =~ s/\*/X/;
		$hotspots{"$gene\t$pos$m"} = 1;
	}
}
close IN;
open IN , "/share/work1/wangrr/DB/363/must.txt";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$F[1] =~ s/\*/X/;
	$hotspots{"$F[0]\t$F[1]"} = 1;
	for my $s (split /\// , $F[6]){
		$hotspots{"$F[0]\t$s"} = 1;
	}
}
close IN;
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

#################################################################################################################################
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
my $syn_var="$pname"."_Somatic_Synonymous.xls";
open(REPORT,">$outdir/$rep_var")||die "$!";
open(REPORTDET,">$outdir/$repdetail_var")||die "$!";
open(SYN,">$outdir/$syn_var")||die "$!";
my $fn = 0;
print REPORT "SampleID\tGeneSymbol\tMajorTranscriptID\tExon#\tcDNA_Change\tAA_Change\tRegion\tMutFunc\tReadDepth_Total\tReadDepth_Var\tVarAlleleFrac\tGeneDescrip\tGeneSummary_NCBI\tClinVar\tTSG_Score\tTSG_Rank\tOG_Score\tOG_Rank\n";
for my $file (@input){
	open(FILTER,"$file")||die "$!";
	$fn++;
	while(defined($a=<FILTER>)){
		$ti="\.";$exon="\.";$cdna="\.";$aa="\.";
		chomp $a;
		@b=split(/\t/,$a);
		if($b[0] eq "SampleName" and $fn == 1){
			print REPORTDET "$a\n";
		}
		my $ishotspots = 0;
		for my $ss (split /,/ , $b[15]){
			my ($trans , $exon , $chgvs , $phgvs);
			for my $sc (split /:/ , $ss){
				if ($sc =~ /NM_/){
					$trans = $sc;
				}elsif ($sc =~ /exon/){
					$exon = $sc;
				}elsif ($sc =~ /c\./){
					$chgvs = $sc;
				}elsif ($sc =~ /p\./){
					$phgvs = $sc;
				}
			}
			if ($phgvs =~ /p\.(\w\d+.)/){
				if (exists $hotspots{"$b[12]\t$1"}){
					$ishotspots = 1;
					last;
				}
			}
			if (exists $hotspots{"$b[12]\t$chgvs"}){
				$ishotspots = 1;
				last;
			}
            if ($b[12] eq 'EGFR' and $exon eq 'exon19' and $chgvs =~ /del/){
                    $ishotspots = 1;
                    last;
            }
		}
	
		next if ((not $b[11]=~/exonic/) and (not $b[11]=~/splicing/));
		next if ($b[14] eq "unknown");
		next if  $b[11] eq 'ncRNA_exonic';
		next if ($b[43]>$popfreq);
		next if (($b[7]+$b[8])<$totalD or $b[8]<$varD);
		if ($ishotspots and $var_freq > 0.01){
			next if (($b[8]/($b[7]+$b[8]))<0.01);
		}else{
			next if (($b[8]/($b[7]+$b[8]))<$var_freq);
		}


		if ($mode eq "putative_somatic" and $ishotspots == 0){
			next if $b[9]>=0.9;
			next if ($b[9]>=0.2 and ((($tsgr{$b[12]}>500 or (not $tsgr{$b[12]})) and ($ogr{$b[12]}>500 or (not $ogr{$b[12]}))) and $b[43]>0));
			next if ($b[9]>=0.2 and ((($tsgr{$b[12]}>500 or (not $tsgr{$b[12]})) and ($ogr{$b[12]}>500 or (not $ogr{$b[12]}))) and ($b[23] ne ".")));
			next if (($ogr{$b[12]}>500 or (not $ogr{$b[12]})) and $b[9]>=0.3 and (($b[49] eq "T") and ($b[51] eq "B") or (($b[49] eq "T") and ($b[53] eq "B")) or (($b[51] eq "B") and ($b[53] eq "B"))) and ($b[43]>0 or $b[44]<=1 or (($b[49] eq "T") and ($b[51] eq "B") and ($b[53] eq "B"))));
			next if $b[23] ne "." and $b[43] > 0.005 and $b[25] eq '.';

		}
		if (($b[14] eq "synonymous SNV") and (not $b[11]=~/splicing/)){
			print SYN "$a\n";
			next;
		}

		@transc=split(/[,;]/,$b[15]);
		$rdtotal=$b[7]+$b[8];
		$rdvar=$b[8];
		$vaf=$rdvar/$rdtotal;
		$i=0;
		@info=();
		$control=0;
		if ($b[12] =~ /[,;]/){
			$b[12] = (sort {length($a)<=>length($b)} (split /[,;]/ , $b[12]))[0];
		}
	
		while($i<@transc){
			@info=split(/:/,$transc[$i]);
			if ((exists $ti{$info[0]} and $ti{$info[0]} eq $info[1]) or (exists $ti{$b[12]} and $ti{$b[12]} eq $info[0])){
				if($#info==4){
					$ti=$info[1];
					$exon=$info[2];
					$cdna=$info[3];
					$aa=$info[4];
				}elsif($#info<4 and $info[1]=~/^exon/){
					$ti=$info[0];
					$exon=$info[1];
					$cdna=$info[2];
					$aa="\.";
				}
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
			}elsif($#info<4 and $info[1]=~/^exon/){
				$ti=$info[0];
				$exon=$info[1];
				$cdna=$info[2];
				$aa="\.";
			}elsif($#info<4 and (not $info[1]=~/^exon/)){
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
		if ($b[14] eq 'synonymous SNV' and $aa =~ /(\w)\d+(\w)/){
			next if $1 eq $2;
			print STDERR "$1\t$2\t$a\n";
		}
	
		next if $cdna eq '' or $cdna eq '.' or $cdna =~ /UTR/;
		print REPORT "$b[0]\t$b[12]\t$ti\t$exon\t$cdna\t$aa\t$b[11]\t$b[14]\t$rdtotal\t$rdvar\t$vaf\t$b[91]\t$summary\t$clv\t";
		print REPORTDET "$a\n";
		if($tsgr{$b[12]}){
			print REPORT "$tsgscore{$b[12]}\t$tsgr{$b[12]}\t$ogscore{$b[12]}\t$ogr{$b[12]}\t";
		}else{
			print REPORT "NA\tNA\tNA\tNA\t";
		}
		print REPORT "$ishotspots\t$b[1]\t$b[2]\t$b[3]\t$b[4]\t$b[5]\n";
	}
}
close(FILTER);close(REPORT);
$a=$ti=$exon=$cdna=$aa="";@b=();@transc=();$rdtotal=$rdvar=$vaf=$i=0;@info=();%sum=%clv=();
