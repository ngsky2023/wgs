#!/usr/bin/env perl
use Getopt::Long;
use FindBin '$Bin';

my $usage =
"
Usage:
	 
Options:
	-inputvcf       <FILE>  The raw vcf file to be filtered, required.
	-outputvcf      <NAME>  Prefix of filtered vcf file, required.
	-Caller         <NAME>  Name of variant caller, either varscan or mutect2, mutect2single, required. 
	-T_depth        <NUM>   Minimum read depth in tumor for a locus to be considered as informative, default is 8. 
	-T_depth_Var    <NUM>   Minimum number of variant allele reads in tumor for a locus to be considered as variant, default is 4.
	-N_depth        <NUM>	Minimum read depth in normal for a locus to be considered as informative, default is 8.
	-N_depth_Var    <NUM>	Maximum number of variant allele reads in normal allowed for a somatic call, default is 2.
	-T_Var_Strand   <NUM>	Minimum number of variant allele reads on each strand in tumor for a locus to be considered as variant, default is 2.
	-T_Var_Frac     <NUM>   Minimum fraction of variant allele reads in tumor for a locus to be considered as variant, default is 0.01.
	-SomaticP       <NUM>   Somatic p value cutoff for one tailed Fisher exact test (which only will be applied to varscan caller). 

For example:
perl $0 -inputvcf input.vcf -outputvcf filtered -Caller varscan -T_depth 8 -T_depth_Var 4 -N_depth 8 -N_depth_Var 2 -T_Var_Strand 2 -T_Var_Frac 0.01 -SomaticP 0.05

<You will get 3 outputs from the above command line: filtered.vcf, filtered_SNV.vcf, and filtered_INDEL.vcf>
";
my ($inputvcf,$outputvcf,$Caller,$T_depth,$T_depth_Var,$N_depth,$N_depth_Var,$T_Var_Strand,$T_Var_Frac,$somaticp);
$T_depth ||=8;
$T_depth_Var ||=4;
$N_depth ||=8;
$N_depth_Var ||=2;
$T_Var_Strand ||=2;
$T_Var_Frac ||=0.01;
$N_Var_Frac ||=0.1;
$somaticp ||=0.05;


GetOptions(
	"inputvcf:s"=>\$inputvcf,
	"outputvcf:s"=>\$outputvcf,
	"Caller:s"=>\$Caller,
	"T_depth:s"=>\$T_depth,
	"T_depth_Var:s"=>\$T_depth_Var,
	"N_depth:s"=>\$N_depth,
	"N_depth_Var:s"=>\$N_depth_Var,
	"T_Var_Strand:s"=>\$T_Var_Strand,
	"T_Var_Frac:s"=>\$T_Var_Frac,
	"N_Var_Frac:s"=>\$N_Var_Frac,
	"SomaticP:s"=>\$somaticp
);
die "$usage\n" if(!$inputvcf || !$outputvcf || !$Caller);
chomp ($inputvcf,$outputvcf,$Caller,$T_depth,$T_depth_Var,$N_depth,$N_depth_Var,$T_Var_Strand,$T_Var_Frac);

my %must;
open IN , "$Bin/must_csmart.txt";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$F[0] = "chr$F[0]" if $F[0] !~ /chr/;
	for my $p (($F[1]-1)..$F[2]){
		$must{"$F[0]\t$p"} = 1;
	}
}
close IN;

#####################################################################################################################
if ($inputvcf =~ /\.gz$/){
	open(IN,"zcat $inputvcf|")||die "$!";
}else{
	open(IN,"$inputvcf")||die "$!";
}
my $combined="$outputvcf".".vcf";
my $snv="$outputvcf"."_SNV.vcf";
my $indel="$outputvcf"."_INDEL.vcf";
open(OUT,">$combined")||die "$!";
open(SNV,">$snv")||die "$!";
open(INDEL,">$indel")||die "$!";
my $in="";
while(defined($in=<IN>)){
	my @in=@ref=@var=@spv=@info=@param=@n=@t=@nd=@td=@nstrand=@tstrand=();
	my $ndp=$tdp=$nvar=$tvar=$tvstr=$tvf=$nvf=$p=0;
	my $ecnt=$hcnt=$maxed=$mined=$rpa=$ru="";
	chomp $in;
	if($in=~/^##/){
		print OUT "$in\n";
		print SNV "$in\n";
		print INDEL "$in\n";
		next;
	}elsif ($in =~ /^#/){
		$in = join("\t" , qw/#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT	TUMOR	NORMAL/);
		print OUT "$in\n";
		print SNV "$in\n";
		print INDEL "$in\n";
		next;
	}else{
		@in=split(/\t/,$in);
		@ref=split(//,$in[3]);
		@var=split(//,$in[4]);
		if($Caller eq "varscan"){
			@spv=split(/SPV=/,$in[7]);
			next if $spv[1]>=$somaticp;
			@n=split(/:/,$in[9]);
			next if $n[0] ne "0/0";
			@t=split(/:/,$in[10]);
			@nstrand=split(/,/,$n[6]);
			@tstrand=split(/,/,$t[6]);
			$ndp=$n[2];
			$tdp=$t[2];
			next if ($tdp==0 or $ndp==0);
			$nvar=$nstrand[2]+$nstrand[3];
			$tvar=$tstrand[2]+$tstrand[3];
			if($tstrand[2]<=$tstrand[3]){
				$tvstr=$tstrand[2];
			}else{
				$tvstr=$tstrand[3];
			}
			$tvf=$tvar/$tdp;
			$nvf=$nvar/$ndp;
		}elsif ($Caller eq "mutect2"){
			if($in[6] eq "PASS"){
				print OUT "$in\n";
				if($#ref==$#var){
					print SNV "$in\n";
				}else{
					print INDEL "$in\n";
				}
				next;
			}else{
				@param=split(/;/,$in[7]);
				$p=0;
				my $ecnt=$hcnt=$maxed=$mined=$rpa=0;
				my $ru="";
				while($p<=$#param){
					@info=split(/[=,]/,$param[$p]);
					if($info[0] eq "ECNT"){
						$ecnt=$info[1];  ## Number of events in the specific haplotype, to check if the variant is located in clustered variants
					}elsif($info[0] eq "HCNT"){
						$hcnt=$info[1];  ## Number of haplotypes supporting this variant
					}elsif($info[0] eq "MAX_ED"){
						$maxed=$info[1];  ## Maximum distance between events in this active region
					}elsif($info[0] eq "MIN_ED"){
						$mined=$info[1];  ## Minimum distance between events in this active region
					}elsif($info[0] eq "RPA"){
						$rpa=$info[1];  ## Number of times tandem repeat unit is repeated
					}elsif($info[0] eq "RU"){
						$ru=$info[1];  ## The sequence of repeated unit
					}
					$p++;
				}
				next if ($ecnt>=4);
				next if ($ecnt>=2 and $mined<=20);
				next if ($in[6]=~/clustered_events/ and $in[6]=~/normal/);
				$in[3]=~s/^.//;$in[4]=~s/^.//;
				next if ($in[6]=~/triallelic_site/ or $in[6]=~/homologous_mapping_event/ or $in[6]=~/multi_event_alt_allele_in_normal/ or $in[6]=~/panel_of_normals/);
				@n=split(/:/,$in[10]);
				@t=split(/:/,$in[9]);
				@nd=split(/,/,$n[1]);
				@td=split(/,/,$t[1]);
				$ndp=$nd[0]+$nd[1];
				$tdp=$td[0]+$td[1];
				next if ($tdp==0 or $ndp==0);
				$nvar=$nd[1];
				$tvar=$td[1];
				if($t[3]<=$t[4]){
					$tvstr=$t[3];
				}else{
					$tvstr=$t[4];
				}
				next if ($rpa>=3 and (($ru eq $in[3]) or ($ru eq $in[4])) and (($tvar/$tdp)<0.1 or $nvar>=1));
				next if ($rpa>5 and (($tvar/$tdp)<0.1 or $nvar>=1));
				next if ($in[6]=~/t_lod_fstar/ and $tvar<=10 and $nvar>=1);
			}
			$tvf=$tvar/$tdp;
			$nvf=$nvar/$ndp;
		}elsif ($Caller eq "gatk4"){
			@param=split(/;/,$in[7]);
			$p=0;
			my $ecnt=$hcnt=$maxed=$mined=$rpa=0;
			my $ru="";
			while($p<=$#param){
				@info=split(/[=,]/,$param[$p]);
				if($info[0] eq "ECNT"){
					$ecnt=$info[1];  ## Number of events in the specific haplotype, to check if the variant is located in clustered variants
				}elsif($info[0] eq "HCNT"){
					$hcnt=$info[1];  ## Number of haplotypes supporting this variant
				}elsif($info[0] eq "MAX_ED"){
					$maxed=$info[1];  ## Maximum distance between events in this active region
				}elsif($info[0] eq "MIN_ED"){
					$mined=$info[1];  ## Minimum distance between events in this active region
				}elsif($info[0] eq "RPA"){
					$rpa=$info[1];  ## Number of times tandem repeat unit is repeated
					my @rpa = split /,/ , $rpa;
					$rpa = (sort {$b<=>$a} @rpa)[0];
				}elsif($info[0] eq "RU"){
					$ru=$info[1];  ## The sequence of repeated unit
				}
				$p++;
			}
			my ($rs , $as) = @in[3,4];
			my $rus = '';
			if (length($rs) > length($as)){
				$rus = $rs;
				$rus =~ s/$as//;
			}elsif (length($rs) < length($as)){
				$rus = $as;
				$rus =~ s/$rs//;
			}
			if ($in[9] =~ /^0\/0/){
				($in[9],$in[10]) = ($in[10],$in[9]);
				$in = join("\t" , @in);
			}
			@t=split(/:/,$in[9]);
			@n=split(/:/,$in[10]);
			@nd=split(/,/,$n[1]);
			@td=split(/,/,$t[1]);
			if ($#td == 2){
			        if ($td[1] > $td[2]){
			                $in[4] =~ s/,.+$//;
			        }else{
			                @td[1,2] = @td[2,1];
			                @nd[1,2] = @nd[2,1];
			                $in[4] =~ s/^.+,//;
			        }
			        $in = join("\t" , @in);
			}elsif ($#td > 2){
			        next;
			}
			$ndp=$nd[0]+$nd[1];
			$tdp=$td[0]+$td[1];
			$nvar=$nd[1];
			$tvar=$td[1];
			$t[4] =~ s/^\d+,//;
			$t[5] =~ s/^\d+,//;
			if($t[4]<=$t[5]){
				$tvstr=$t[4];
			}else{
				$tvstr=$t[5];
			}
			next if $tdp == 0 or $ndp == 0;
			$tvf=$tvar/$tdp;
			$nvf=$nvar/$ndp;
			#next if ($rpa>=3 and ($ru eq $rus or $ru x 2 eq $rus or $ru x 3 eq $rus) and (($tvar/$tdp)<0.1 or $nvar>=1));
			if($in[6] eq "PASS" or exists $must{"$in[0]\t$in[1]"}){
				$tvstr=10 if $tvstr < 10 and exists $must{"$in[0]\t$in[1]"};
			}else{
				#next if ($ecnt>=2 and $mined<=20);
				next if ($in[6] =~ /contamination|panel_of_normals/);
				next if ($in[6] =~ /fragment_length|mapping_quality|strand_artifact|base_quality|read_position/);
				next if ($tdp==0 or $ndp==0);
			}
		}elsif ($Caller eq "gatk4single"){
			if($in[6] eq "PASS"){
				print OUT "$in\n";
				if($#ref==$#var){
					print SNV "$in\n";
				}else{
					print INDEL "$in\n";
				}
				next;
			}else{
				@param=split(/;/,$in[7]);
				$p=0;
				my $ecnt=$hcnt=$maxed=$mined=$rpa=0;
				my $ru="";
				while($p<=$#param){
					@info=split(/[=,]/,$param[$p]);
					if($info[0] eq "ECNT"){
						$ecnt=$info[1];  ## Number of events in the specific haplotype, to check if the variant is located in clustered variants
					}elsif($info[0] eq "HCNT"){
						$hcnt=$info[1];  ## Number of haplotypes supporting this variant
					}elsif($info[0] eq "MAX_ED"){
						$maxed=$info[1];  ## Maximum distance between events in this active region
					}elsif($info[0] eq "MIN_ED"){
						$mined=$info[1];  ## Minimum distance between events in this active region
					}elsif($info[0] eq "RPA"){
						$rpa=$info[1];  ## Number of times tandem repeat unit is repeated
					}elsif($info[0] eq "RU"){
						$ru=$info[1];  ## The sequence of repeated unit
					}
					$p++;
				}
				next if ($ecnt>=4);
				#next if ($ecnt>=2 and $mined<=20);
				next if ($in[6]=~/clustered_events/ and ($in[6]=~/germline_risk/ or $in[6]=~/normal/));
				$in[3]=~s/^.//;$in[4]=~s/^.//;
				next if ($in[6]=~/multiallelic/ or $in[6]=~/homologous_mapping_event/ or $in[6]=~/multi_event_alt_allele_in_normal/ or $in[6]=~/panel_of_normals/);

				@t=split(/:/,$in[9]);
				@td=split(/,/,$t[1]);
				$ndp=$N_depth+1;
				$tdp=$td[0]+$td[1];
				next if ($tdp==0 or $ndp==0);
				$nvar=0;
				$tvar=$td[1];
				if($t[3]<=$t[4]){
					$tvstr=$t[3];
				}else{
					$tvstr=$t[4];
				}
				next if ($rpa>=3 and (($ru eq $in[3]) or ($ru eq $in[4])) and (($tvar/$tdp)<0.1 or $nvar>=1));
				next if ($rpa>5 and (($tvar/$tdp)<0.1 or $nvar>=1));
				next if ($in[6]=~/t_lod/ and $tvar<=10 and $nvar>=1);
			}
			$tvf=$tvar/$tdp;
			$nvf=$nvar/$ndp;
		}elsif ($Caller eq "mutect2single"){
			if($in[6] eq "PASS"){
				print OUT "$in\n";
				if($#ref==$#var){
					print SNV "$in\n";
				}else{
					print INDEL "$in\n";
				}
				next;
			}else{
				@param=split(/;/,$in[7]);
				$p=0;
				my $ecnt=$hcnt=$maxed=$mined=$rpa=0;
				my $ru="";
				while($p<=$#param){
					@info=split(/[=,]/,$param[$p]);
					if($info[0] eq "ECNT"){
						$ecnt=$info[1];  ## Number of events in the specific haplotype, to check if the variant is located in clustered variants
					}elsif($info[0] eq "HCNT"){
						$hcnt=$info[1];  ## Number of haplotypes supporting this variant
					}elsif($info[0] eq "MAX_ED"){
						$maxed=$info[1];  ## Maximum distance between events in this active region
					}elsif($info[0] eq "MIN_ED"){
						$mined=$info[1];  ## Minimum distance between events in this active region
					}elsif($info[0] eq "RPA"){
						$rpa=$info[1];  ## Number of times tandem repeat unit is repeated
					}elsif($info[0] eq "RU"){
						$ru=$info[1];  ## The sequence of repeated unit
					}
					$p++;
				}
				next if ($ecnt>=4);
				next if ($ecnt>=2 and $mined<=20);
				next if ($in[6]=~/clustered_events/ and $in[6]=~/normal/);
				$in[3]=~s/^.//;$in[4]=~s/^.//;
				next if ($in[6]=~/triallelic_site/ or $in[6]=~/homologous_mapping_event/ or $in[6]=~/multi_event_alt_allele_in_normal/ or $in[6]=~/panel_of_normals/);
				@t=split(/:/,$in[9]);
				@td=split(/,/,$t[1]);
				$ndp=$N_depth+1;
				$tdp=$td[0]+$td[1];
				next if ($tdp==0 or $ndp==0);
				$nvar=0;
				$tvar=$td[1];
				if($t[3]<=$t[4]){
					$tvstr=$t[3];
				}else{
					$tvstr=$t[4];
				}
				next if ($rpa>=3 and (($ru eq $in[3]) or ($ru eq $in[4])) and (($tvar/$tdp)<0.1 or $nvar>=1));
				next if ($rpa>5 and (($tvar/$tdp)<0.1 or $nvar>=1));
				next if ($in[6]=~/t_lod_fstar/ and $tvar<=10 and $nvar>=1);
			}
			$tvf=$tvar/$tdp;
			$nvf=$nvar/$ndp;
		}
		if($tdp>=$T_depth and $tvar>=$T_depth_Var and $ndp>=$N_depth and $nvar<=$N_depth_Var and $tvstr>=$T_Var_Strand and $tvf>=$T_Var_Frac and $nvf<$tvf/3){
		#if($tdp>=$T_depth and $tvar>=$T_depth_Var and $ndp>=$N_depth and $tvstr>=$T_Var_Strand and $tvf>=$T_Var_Frac and $nvf<$N_Var_Frac and $nvf<$tvf/3){
			print OUT "$in\n";
			if($#ref==$#var){
				print SNV "$in\n";
			}else{
				print INDEL "$in\n";
			}
		}else{
			print STDERR "filter\t$in\n";
		}
	}
}
close(IN);close(OUT);close(SNV);close(INDEL);
$in="";@in=@spv=@info=@param=@n=@t=@nd=@td=@nstrand=@tstrand=@ref=@var=();$ecnt=$hcnt=$maxed=$mined=$rpa=$ru="";
