#!usr/bin/perl
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
use FindBin '$Bin';

my $usage =
"
Usage:
	 
Options:
	-tumor    <FILE>  The depth data of tumor samples (the .sample_interval_summary output of GATK DepthOfCoverage), required.
	-normal   <FILE>  The depth data of normal pool (the .sample_interval_summary output of GATK DepthOfCoverage), required.
	-bed      <FILE>  The bed file defines the targeted regions with gene assignment column, required.
	-output   <DIR>   Output directory, required.
	-project  <STR>   The name you want to assign to this analysis (as the prefix of output files), required. 
	-num_norm <NUM>   Number of normal samples used in the analysis, required. 
	-help

For example:
perl $0 -tumor /share/work1/staff/lijianbai/Tests/GeneCNV_Scripts/GeneCNV_V2/TumorCoverage/Tumor_62Gene.sample_interval_summary -normal /share/work1/staff/lijianbai/Tests/GeneCNV_Scripts/GeneCNV_V2/BloodCoverage/Blood_62Gene.sample_interval_summary -bed /home/lijianbai/workdir/BEDs/IDT39_InheritedCancer_v1_ext.bed -output /share/work1/staff/lijianbai/Tests/GeneCNV_Scripts/GeneCNV_V2/ -project BCH_62gene
";

my ($tumor,$normal,$bed,$output,$project,$num_norm,$help);
GetOptions(
	"tumor:s"=>\$tumor,
	"normal:s"=>\$normal,
	"bed:s"=>\$bed,
	"output:s"=>\$output,
	"project:s"=>\$project,
	"num_norm:s"=>\$num_norm,
	"help:s"=>\$help
);
chomp $project;
die "$usage\n" if(!$tumor || !$normal || !$bed || !$output || !$project || !$num_norm ||$help);
if($num_norm<300){
	$factor=300/$num_norm;
}else{
	$factor=1;
}
####################################################################################################
print "Converting data format and sorting tumor coverage data along with normal pool ...\n";
open(BED,"$bed")||die "$!";
while($region=<BED>){
	chomp $region;
	@region=split(/\t/,$region);
	$pos="$region[0]\t$region[2]";
	$gene{$pos}=$region[3];
}
close(BED);
$region=$pos="";@region=();$start=0;

open(NORM,"$normal")||die "$!";
open(NORM_MOD,">$output/$project\_NormalCoverage.txt")||die "$!";
$line=0;
while($norm=<NORM>){
	chomp $norm;
	$line++;
	@norm=split(/\t/,$norm);
	$col=3;
	if($line==1){
		print NORM_MOD "Target\tGene";
		while($col<@norm){
			if($norm[$col]=~/_mean_cvg/){
				@sample=split(/_mean_cvg/,$norm[$col]);
				print NORM_MOD "\t$sample[0]";
				$index{$col}=1;
			}
			$col++;
		}
		print NORM_MOD "\n";
	}else{
		@region=split(/[:-]/,$norm[0]);
		$pos="$region[0]\t$region[2]";
		print NORM_MOD "$norm[0]\t$gene{$pos}";
		while($col<@norm){
			if($index{$col}){
				print NORM_MOD "\t$norm[$col]";
			}
			$col++;
		}
		print NORM_MOD "\n";
	}
}
close(NORM);close(NORM_MOD);
$line=$col=0;$norm=$pos="";@norm=@sample=@region=();%index=();

open(TUM,"$tumor")||die "$!";
while($tum=<TUM>){
	chomp $tum;
	@tum=split(/\t/,$tum);
	$line++;
	$col=3;
	if($line==1){
		push @header,"Target\tGene";
		while($col<@tum){
			if($tum[$col]=~/_mean_cvg/){
				@sample=split(/_mean_cvg/,$tum[$col]);
				push @header,"\t$sample[0]";
				$index{$col}=1;
			}
			$col++;
		}
		$header=join("",@header);
	}else{
		@region=split(/[:-]/,$tum[0]);
		$pos="$region[0]\t$region[2]";
		@tumor=();
		push @tumor,"$tum[0]\t$gene{$pos}";
		while($col<@tum){
			if($index{$col}){
				push @tumor,"\t$tum[$col]";
			}
			$col++;
		}
		$tumor{$tum[0]}=join("",@tumor);
	}
}
close(TUM);
$line=$col=0;$tum="";@tum=@sample=@tumor=@header=@region=();%index=();

open(NORM,"$normal")||die "$!";
open(TUM_MOD,">$output/$project\_TumorCoverage.txt")||die "$!";
while($norm=<NORM>){
	chomp $norm;
	@norm=split(/\t/,$norm);
	if($norm[0] eq "Target"){
		print TUM_MOD "$header\n";
	}else{
		if($tumor{$norm[0]}){
			print TUM_MOD "$tumor{$norm[0]}\n";
		}else{
			die "Region $norm[0] does not exist in tumor data!\n";
		}
	}
}
close(NORM);close(TUM_MOD);
$line=$col=0;$norm=$header="";@norm=();%tumor=();
##################################################################################################
print "Normalizing depth using mean depth, and ranking tumor depth against normalized normal pool ...\n";
system("/share/public/software/R/3.4.0/bin/Rscript $Bin/Ranking.R $output/$project\_TumorCoverage.txt $output/$project\_NormalCoverage.txt $output/$project\_TumorNormalized_Rank.txt");
##################################################################################################
print "Calculating the possibility of with/without CNV for each gene ...\n";
open(RANK,"$output/$project\_TumorNormalized_Rank.txt")||die "$!";
open(PA,">$output/$project\_TumorNormalized_PossibilityAmp.txt")||die "$!";
open(PD,">$output/$project\_TumorNormalized_PossibilityDel.txt")||die "$!";
while($rank=<RANK>){
	chomp $rank;
	@rank=split(/\t/,$rank);
	$col=3;
	if($rank[0] eq "Target"){
		print PA "Gene\tNumTargets";
		print PD "Gene\tNumTargets";
		while($col<@rank){
			@sample=split(/_/,$rank[$col]);
			$sample{$col}=$sample[0];
			print PA "\t$sample[0]";
			print PD "\t$sample[0]";
			push(@samples,$sample[0]);
			$col=$col+2;
		}
		print PA "\n";
		print PD "\n";
	}else{
		$target{$rank[1]}++;
		while($col<@rank){
			if($sample{$col}){
				$genesamp="$rank[1]\t$sample{$col}";
				if($rank[$col]==1){
					$rank[$col]=0.99;
				}else{
					$rank[$col]=$rank[$col];
				}
				if(not $possibilityA{$genesamp}){
					$possibilityA{$genesamp}=$rank[$col];
					$rpA{$genesamp}=$factor*(1-$rank[$col]);
					$possibilityD{$genesamp}=1-$rank[$col];
					$rpD{$genesamp}=$factor*$rank[$col];
				}else{
					$possibilityA{$genesamp}=$possibilityA{$genesamp}*$rank[$col];
					$rpA{$genesamp}=$rpA{$genesamp}*$factor*(1-$rank[$col]);
					$possibilityD{$genesamp}=$possibilityD{$genesamp}*(1-$rank[$col]);
					$rpD{$genesamp}=$rpD{$genesamp}*$factor*$rank[$col];
				}
			}
			$col=$col+2;
		}
	}
}
close(RANK);
$rank=$genesamp="";$col=0;@rank=@sample=();
@genes=sort keys %target;
$i=$j=0;
while($i<@genes){
	print PA "$genes[$i]\t$target{$genes[$i]}";
	print PD "$genes[$i]\t$target{$genes[$i]}";
	$j=0;
	while($j<@samples){
		$genesamp="$genes[$i]\t$samples[$j]";
		$pA=$possibilityA{$genesamp}/($possibilityA{$genesamp}+$rpA{$genesamp});
		$pD=$possibilityD{$genesamp}/($possibilityD{$genesamp}+$rpD{$genesamp});
		print PA "\t$pA";
		print PD "\t$pD";
		$j++;
	}
	print PA "\n";
	print PD "\n";
	$i++;
}
close(PA);close(PD);
$i=$j=$pA=$pD=0;@genes=@samples=();%target=%sample=%possibilityA=%rpA=%possibilityD=%rpD=();
