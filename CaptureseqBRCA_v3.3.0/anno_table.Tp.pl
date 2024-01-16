#modified by zhangyu, 2016.0331
#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage =
"
Usage:

  Options:
  -i|input      <FILE>   Input file of the result from annovar.
  -n|name	<STR>	 Sample name.
  -o|out        <FILE>   Output file.
  -a|anno       <FILE>   The Go, KEGG, GeneInfo file.
  -s|soft       <STR>    Software of call SNV/Indel, GATK or Samtools or Mixture? default is 'Mixture'.
  -bed          <FILE>   The bed-file of target region.
  -h|help                Help

For example:
        perl $0 -i T_snp.annovar.hg19_multianno.txt -o T_snp.annovar.xls -a 9606.geneid2go_kegg.xls -s GATK -bed S04380110_Regions.bed
";

my ($input,$name,$outfile,$anno,$soft,$bed,$help);
GetOptions(
  "i|input=s"=>\$input,
  "n|name=s"=>\$name,
  "o|out=s"=>\$outfile,
  "a|anno=s"=>\$anno,
  "s|soft=s"=>\$soft,
  "bed=s"=>\$bed,
  "h|help"=>\$help
);
$soft ||= 'Mixture';
$anno ||= '/share/public/database/KEGG/9606.geneid2go_kegg.xls';
$name ||= 'commit';

if (!$input || !$outfile || !$bed || $help){
        die "$usage\n";
}
if (($soft !~ /GATK/i) & ($soft !~ /samtools/i) & ($soft !~ /Mixture/i) & ($soft !~ /Somatic/i)){
	die "-soft erro...\n$usage\n";
}

##read go_kegg_geneInfo
my %descript;
open INFO,"$anno" || die "$!";
while(<INFO>){
	chomp;
	my @tmp = split /\t/;
	my $gene_symbol=$tmp[1];
	my $des=$tmp[3];
	my $bp=$tmp[4];
	my $cc=$tmp[5];
	my $mf=$tmp[6];
	my $pathway=$tmp[8];
	$descript{$gene_symbol} = [$des,$bp,$cc,$mf,$pathway];
}
close INFO;

##read bed
my %bed_region = ();
open BED,"$bed" || die "$!";
while(my $line = <BED>){
	chomp($line);
	my @arr = split /\t+/,$line;
	next if (@arr < 3);
	$bed_region{$arr[0]}{$arr[1]} = $arr[2];
}
close BED;

open IN,"$input" || die "$!";
open TABLE,">$outfile" || die "$!";
my $head = <IN>;
chomp $head;
my @atm=split /\t/,$head;
my $num=$#atm;
print TABLE "SampleName\t";
for (my $i=0;$i<=4;$i++){
	print TABLE "$atm[$i]\t";
}
print TABLE "Genotype\tDepth_ref\tDepth_alt\tAlt_ratio\tTarget_Flank\t";
for (my $i=5;$i<$#atm;$i++){
	print TABLE "$atm[$i]\t";
}
print TABLE "Gene_description\tGO\tKEGG\n";

while(my $line = <IN>){
	chomp($line);
	#print "$line\n";
	my @msg = split /\t/,$line;#print "$msg[-1]\t$msg[-2]\n";
	my ($type,$d1,$d2) = ('',0,0);
	if ($soft =~ /GATK/i){                                       #GT:AD:GQ:PL:TP  0/1:1,3:27:111,0,27:18  0|0:6,0:18:0,18,234:18  1|0:7,2:58:58,0,251:18
		if ($msg[-1] =~ /^(\d\/\d):(\d+),(\d+)/ || $msg[-1] =~ /^(\d\|\d):(\d+),(\d+)/){
			$type = $1;
			$d1 = $2;
			$d2 = $3;
			$type =~ s/\|/\//g;
			if ($type eq '0/1' || $type eq '1/0'){
				$type = 'het';
			}else{
				$type = 'hom';
			}
		}elsif($msg[-3] =~ /^(\d\/\d):(\d+),(\d+)/ || $msg[-3] =~ /^(\d\|\d):(\d+),(\d+)/){
			$type = $1;
			$d1 = $2;
			$d2 = $3;
			$type =~ s/\|/\//g;
                        if ($type eq '0/1' || $type eq '1/0'){
                                $type = 'het';
                        }else{
                                $type = 'hom';
                        }
		}
	}elsif($soft =~ /samtools/i){
		if ($msg[-1] =~ /^(\d\/\d):/ || $msg[-1] =~ /^(\d\|\d):/){
			$type = $1;
			if ($type eq '0/1' || $type eq '1/0'){
				$type = 'het';
			}else{
				$type = 'hom';
			}
		}
		if ($msg[-3] =~ /DP4=(\d+),(\d+),(\d+),(\d+)/i){
			$d1 = $1 + $2;
			$d2 = $3 + $4;
		}
	}elsif($soft =~ /Mixture/i){
                if($msg[-2] =~ /GT:VR:RR:DP:GQ/){
                        my @atm=split /:/,$msg[-1];
                        $type = $atm[0];
			$d1 = $atm[2];
			$d2 = $atm[1];
                        $type =~ s/\|/\//g;
                }elsif($msg[-2] =~ /GT:GL:GOF:GQ:NR:NV/){
                        my @atm=split /:/,$msg[-1];
                        $type = $atm[0];
                        $type =~ s/\|/\//g;
                        $d1 = $atm[4]-$atm[5];
                        $d2 = $atm[5];
                }elsif(($msg[-2] =~ /GT:AD:GQ:PL/) || ($msg[-2] =~ /GT:AD:DP:GQ:PL/)){
                        if($msg[-1] =~ /^(\d\/\d):(\d+),(\d+)/ || $msg[-1] =~ /^(\d\|\d):(\d+),(\d+)/){
                                $type = $1;
				$d1 = $2;
				$d2 = $3;
                                $type =~ s/\|/\//g;
                        }
                }
		if ($type eq '0/1' || $type eq '1/0'){
                        $type = 'het';
                }else{
                        $type = 'hom';
                }
        }elsif($soft =~ /Somatic/){
		if($msg[-2] =~ /GT:GQ:DP:RD:AD:FREQ/){
                        my @atm=split /:/,$msg[-1];
                        $type = $atm[0];
                        $type =~ s/\|/\//g;$type=~s/1\/0/0\/1/g;
                        $d1 = $atm[3];
                        $d2 = $atm[4];
                }if ($type eq '0/1' || $type eq '1/0'){
                        $type = 'het';
                }else{
                        $type = 'hom';
                }
	}#print "$d1\t$d2\n";
	my $ratio=".";
	if(($d1+$d2)>0){
		$ratio=sprintf ("%.2f",$d2/($d1+$d2));
	}
	my $Target_flank = 'Flank';
	foreach my $start(keys %{$bed_region{$msg[0]}}){
		my $end = $bed_region{$msg[0]}{$start};
		if (($msg[1] >= $start)&& ($msg[1] <= $end)){
			$Target_flank = 'Target';
			last;
		}
	}

	#Region\tGene\tFunction\tAAChange\tclinvar_20140211\t1000g2012apr_all\t1000g2014oct_eas\tesp6500si_all\tExAC_ALL\tExAC_EAS\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\n
	my ($Region,$Gene,$Function,$AAChange) = ('','','','');
	if ($msg[5] !~ /exonic/i){
		if($msg[10] !~ /exonic/i){
			#if($msg[15] =~ /exonic/i){
			#	print "CCDS: $_\n";
			#}
			$Region = $msg[5];
			$Gene = $msg[6];
			$Function = $msg[8];
			$AAChange = $msg[9];
		}else{
			$Region = $msg[10];
			$Gene = $msg[11];
			$Function = $msg[13];
			$AAChange = $msg[14];
		}
	}else{
		$Region = $msg[5];
		$Gene = $msg[6];
		$Function = $msg[8];
		$AAChange = $msg[9];
	}

	my $symbol = $Gene;
	my @tmp = split /,|;/,$symbol;
	$symbol = $tmp[0];
	$symbol = $1 if ($symbol =~ /(.+?)\s*\(/);
	
	my ($description,$go,$kegg);
	if (exists $descript{$symbol}){
		$description=$descript{$symbol}->[0];
		$go=$descript{$symbol}->[1]." & ".$descript{$symbol}->[2]." & ".$descript{$symbol}->[3];
		$kegg=$descript{$symbol}->[4];
		#print TABLE "$msg[0]\t$msg[1]\t$msg[2]\t$msg[3]\t$msg[4]\t$type\t$d1\t$d2\t$msg[28]\t$msg[-5]\t$Target_flank\t$Region\t$Gene\t$Function\t$AAChange\t$msg[54]\t$msg[22]\t$msg[24]\t$msg[55]\t$msg[56]\t$msg[59]\t$msg[29]\t$msg[30]\t$msg[31]\t$msg[32]\t$msg[33]\t$msg[34]\t$msg[37]\t$msg[38]\t$descript{$symbol}\n";
	}else{
		$description=".";
		$go="-";
		$kegg="-";
		#print TABLE "$msg[0]\t$msg[1]\t$msg[2]\t$msg[3]\t$msg[4]\t$type\t$d1\t$d2\t$msg[28]\t$msg[-5]\t$Target_flank\t$Region\t$Gene\t$Function\t$AAChange\t$msg[54]\t$msg[22]\t$msg[24]\t$msg[55]\t$msg[56]\t$msg[59]\t$msg[29]\t$msg[30]\t$msg[31]\t$msg[32]\t$msg[33]\t$msg[34]\t$msg[37]\t$msg[38]\t.\n";
	}
	next if($d2 <= 1);#print "$name\tss\n";
	print  TABLE "$name\t";	
	for (my $i=0; $i<=4;$i++){
		print TABLE  "$msg[$i]\t";
	}
	print TABLE "$type\t$d1\t$d2\t$ratio\t$Target_flank\t";
	for (my $i=5; $i<$num;$i++){
		print TABLE "$msg[$i]\t";
	}
	print TABLE "$description\t$go\t$kegg\n";
}
close IN;
close TABLE;

