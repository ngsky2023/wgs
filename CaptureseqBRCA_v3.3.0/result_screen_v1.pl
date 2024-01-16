#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

#Modified:yangrutao

die "perl $0 <in.level> <*indel.annovar.BRCA_*.avinput> <atlas.snp.vcf> <atlas.indel.vcf> <result.xls>" unless @ARGV==5;
`sort -k6 -rn -o $ARGV[0] $ARGV[0]`;
my $info="";
my $pmid="";
my @result;
my $Result;
my %hgvsc;
my %transcript;
my %hs_hgvsp;
my %hs_hgvsp_simp;
my %hs_start;
my %hs_ref;
my %hs_alt;
my %aa = ();  
my %aa123=();
my %atlas_VR;
my %atlas_RR;
my $clinvar_hgvs="/share/work2/yangrutao/workdir/research/databases/55gene.list";
my $main_transcript="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/variant_summary.GRCh37.gene_NM.txt";
my $hgvs="/share/work2/yangrutao/workdir/research/databases/HGVS_NM_NP_clinvar20170202.xls";
my $aa_file="/share/work2/yangrutao/workdir/research/databases/aa.list.csv" ;
load_hgvsc_from_file($clinvar_hgvs, \%hgvsc);
load_aa_from_file($aa_file,\%aa,\%aa123);
load_main_transcript_file($main_transcript,\%transcript);
load_hgvs_file($hgvs,\%hs_hgvsp,\%hs_hgvsp_simp);
load_atlas_file($ARGV[2],\%atlas_VR,\%atlas_RR);
load_atlas_file($ARGV[3],\%atlas_VR,\%atlas_RR);
#print Dumper(%aa123);exit;
#my %func=(
#		  "frameshift deletion" => "移码缺失变异",
#		  "frameshift insertion" => "移码插入变异",
#		  "frameshift substitution"=> "移码替换变异",
#		  "nonframeshift deletion"=> "非移码缺失变异",
#		  "nonframeshift insertion"=> "非移码插入变异",
#		  "nonframeshift substitution"=> "非移码替换变异",
#		  "nonsynonymous SNV"=> "非同义单核苷酸位点变异",
#		  "stopgain"=> "提前获得终止密码子",
#		  "stoploss"=> "终止密码子丢失",
#		  "synonymous SNV"=> "同义单核苷酸位点变异",
#		  "unknown"=> "未知变异",
#		  );
open(IN,"<$ARGV[0]") or die $!;
while(<IN>)
{
	chomp;
	next if ($_ =~ /^#/);
	next if ($_ =~ /^Level/);
	my @temp=split(/\t/,$_,);
	my $clisig = $temp[0];
	unless ($clisig ~~ @result){
		push (@result,$clisig);
	}
}
close IN;
my $result_key;
foreach my $key (@result){
	$result_key .= "$key+";
}
if ($result_key =~ /Pathogenic/){
	$Result = "Pathogenic";
}elsif($result_key =~ /Likely pathogenic/){
	$Result = "Likely pathogenic";
}elsif($result_key =~ /Uncertain significance/){
	$Result = "Uncertain significance";
}elsif($result_key =~ /Likely benign/){
	$Result = "Likely benign";
}elsif($result_key =~ /Benign/){
	$Result = "Benign";
}

open IN1,"<$ARGV[0]";
open IN2,"<$ARGV[1]";
open OUT,">$ARGV[4]";
print OUT "#Sample\tResult\tChr\tPos\tdbSNP\tRef\tAlt\tGene\tDepth_ref\tDepth_alt\tAlt_ratio\tTranscript\tExon\tHGVS(c.)\tNucleoprotein\tHGVS(p.)\tZygosity\tFunc\.refGene\tType\tName\tCliSig\tEvidence\tPMID\tAAchange\tChr\tStart\tEnd\tCytoBand\tMark\n";
while (<IN2>){
	chomp;
	next if ($_ =~ /^SampleName/);
	my @tmp =split /\t/;
	$tmp[0] =~ s/^#//;
	$hs_start{"$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]"}=$tmp[6];
	$hs_ref{"$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]"}=$tmp[8];
	$hs_alt{"$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]"}=$tmp[9];
		
}
close IN2;

while(<IN1>){
	chomp;
	next if ($_ =~ /^Level/);
	next if ($_ =~ /^Uncertain significance\tPVS0PS0PM1PP0BA0BS0BP0_pm2\tSampleName/);
	my @atm=split /\t/,$_;
	my ($function,$gene,$chr,$start,$end,$genedetail,$ExonicFunc,$AAchange,$genotype,$frequence,$clinvar,$level,$evidence)=($atm[13],$atm[14],$atm[3],$atm[4],$atm[5],$atm[15],$atm[16],$atm[17],$atm[8],$atm[28],$atm[26],$atm[0],$atm[1]);
	my ($sample,$trscript,$exon,$hc,$hp,$ref,$alt,$ref_depth,$alt_depth,$alt_ratio)=($atm[2],$atm[96],$atm[97],$atm[98],$atm[99],$atm[6],$atm[7],$atm[9],$atm[10],$atm[11],);
    my ($cytoBand,$dbsnp,$mark)=($atm[24],$atm[25],$atm[100]);
	if ($gene =~ /^MRE11A$/){
		$gene= "MRE11";
	}
	if ($atm[17] ne "." && $atm[98] eq "." && $atm[99] eq "."){
		#print "$atm[17]\n";
		my @tmp=split(/,/,$atm[17]);
		for (my $i=0;$i<=$#tmp;$i++){	
			my @exp=split(/:/,$tmp[$i],);
			for (my $l=0;$l<=$#exp;$l++){
				my $NM=$exp[$l] if ($exp[$l]=~ /NM/);
				my $en=$exp[$l] if ($exp[$l]=~ /exon/);
				my $c=$exp[$l] if ($exp[$l]=~ /c\./);
				my $p=$exp[$l] if ($exp[$l]=~ /p\./);
				if ($tmp[17] =~ /$trscript/){
					$hc = $c;
					$hp = $p;
				}
			}
		}
	}
	if ($atm[98] =~ /^c.([ACGT])(\d+)([ACGTYRN])$/) {
            $atm[98] ="c\.$2$1>$3";
    }
	my $hgvsp;
	if($hp){
     	if($hp=~/p.[a-zA-Z]{1}\d/)
     	{
     	$hgvsp=hgvp_aa1to3($hp);
		}
	}
	my $hc_name=$transcript{$gene}{$trscript} . ":" . $atm[98];
	my $hgvs_p;
	if (defined $hs_hgvsp{$gene}{$transcript{$gene}{$trscript}}){
		$hgvs_p = $hs_hgvsp{$gene}{$transcript{$gene}{$trscript}};
	}else {
		my $trscript_simp =$transcript{$gene}{$trscript};
		if ($trscript_simp =~ /^(NM\_\d+)\.\d+$/){
			$trscript_simp = $1 ;
			$hgvs_p = $hs_hgvsp_simp{$gene}{$trscript_simp};
		}
	}
	my $name;
	if ($hgvs_p ne "." && $hp ne "." ){
		$name = $transcript{$gene}{$trscript} . "\($gene\)" . ":". $atm[98] . "\($hgvsp\)";
	}else{
		$name = $transcript{$gene}{$trscript} . ":" . "\($gene\)" . $atm[98] ;
	}
	if(exists($hgvsc{$atm[14]}{$hc_name}{"info"})){
		if($hgvsc{$atm[14]}{$hc_name}{"pmid"}){
			$pmid=$hgvsc{$atm[14]}{$hc_name}{"pmid"};
		}else{
			$pmid = ".";
		}
	}else{
			$pmid = ".";
	}
	if ($ExonicFunc !~/SNV/){
		my $string="$chr\t$start\t$ref\t$alt";
		if (defined $hs_start{$string}){
			$start = $hs_start{$string};
			$ref = $hs_ref{$string};
			$alt = $hs_alt{$string};
		}
		if (defined $atlas_VR{$string}  &&  defined $atlas_RR{$string}){
			$ref_depth = $atlas_RR{$string};
			$alt_depth = $atlas_VR{$string};
			$alt_ratio = $alt_depth / ($ref_depth + $alt_depth);
			$alt_ratio = sprintf "%.2f",$alt_ratio;
		}
	}
	if ($ExonicFunc =~/SNV/){
		my $string="$chr\t$start\t$ref\t$alt";
		if (defined $atlas_VR{$string}  &&  defined $atlas_RR{$string}){
			$ref_depth = $atlas_RR{$string};
			$alt_depth = $atlas_VR{$string};
			$alt_ratio = $alt_depth / ($ref_depth + $alt_depth);
			$alt_ratio = sprintf "%.2f",$alt_ratio;
		}
	}
	if($function !~/intronic|UTR/){
		if ($name =~ /MRE11A/){
			$name =~ s/MRE11A/MRE11/g;
		}
		print OUT "$atm[2]\t$Result\t$chr\t$start\t$dbsnp\t$ref\t$alt\t$gene\t$ref_depth\t$alt_depth\t$alt_ratio\t$transcript{$gene}{$trscript}\t$exon\t$hc\t$hgvs_p\t$hgvsp\t$genotype\t$function\t$ExonicFunc\t$name\t$level\t$evidence\t$pmid\t$AAchange\t$chr\t$start\t$end\t$cytoBand\t$mark\n";

    }
}
close IN1;
`sort -k6 -n -o $ARGV[0] $ARGV[0]`;
#---------------------------------
sub load_hgvsc_from_file{
	my ($fhgvs, $hgvsc) = @_;
	open(FAA,$fhgvs);
	my $line=<FAA>;
	while($line=<FAA>)
	{
		chomp($line);
		my @temp=split(/\t/,$line,);
		$hgvsc->{$temp[0]}{$temp[1]}{"info"}=$line;
		$hgvsc->{$temp[0]}{$temp[1]}{"pmid"}=$temp[3];
	}
	close FAA;
}

sub load_main_transcript_file{
	my ($ftranscript,$transcript) = @_ ;
	open(TRS,$ftranscript);
	while(<TRS>){
		chomp;
		next if ($_ =~ /^\s*$/);
		my @tmp = split /\t/,$_;
		if ($tmp[0] =~ /^MRE11A$/){
			$tmp[0] = "MRE11";
		}
		$transcript{$tmp[0]}{$tmp[1]}=$tmp[1] . "\." . $tmp[2]; 
	}
	close TRS;
}

sub load_hgvs_file{
	my ($fhgvs,$hs_hgvsp,$hs_hgvsp_simp) = @_ ;
	open(HGVS,$fhgvs);
	while(<HGVS>){
		chomp;
		next if ($_ =~ /^\s*$/);
		my @tmp = split /\t/,$_;
		if ($tmp[1] =~ /^MRE11A$/){
			$tmp[1] = "MRE11";
		}
		$hs_hgvsp{$tmp[1]}{$tmp[2]}=$tmp[3];
		if ($tmp[2] =~ /^(NM\_\d+)\.\d+$/) {
			$hs_hgvsp_simp{$tmp[1]}{$1} = $tmp[3];
		}
	}
	close HGVS;
}

sub load_aa_from_file{   # load_aa_from_file($aa_file,\%aa,\%aa123);

	my ($faa, $aaname, $aa123) = @_;
	open(FAA,$faa);
	my $line=<FAA>;
	while($line=<FAA>)
	{
	   my @tmp=split /\,/,$line;
	   #$aaname{$tmp[0]}=$tmp[2];
		$aaname->{$tmp[0]}=$tmp[2];
	    $aa123->{$tmp[1]}=$tmp[0];
	}
	close(FAA);

}

sub hgvp_aa1to3{
	my ($hgvp)=@_;
	#print $hgvp,"\n";
	$hgvp=~/p.(.*)/;
	my $info=$1;
	my $info3="";
	for (my $i=0;$i<length($info);$i++)
	{
	my $char=substr($info,$i,1);
	if($char=~/[A-Z]/){
		$info3.=$aa123{$char};
	}
	else
	{$info3.=$char;}
	}

	$hgvp="p.".$info3;
	return $hgvp;

}

sub load_atlas_file{

	my ($atlas,$atlas_VR,$atlas_RR) = @_;
	open(FA,$atlas);
	while(<FA>)
	{
		chomp;
		next if ($_ =~ /^\s*$/);
		next if ($_ =~ /^#/);
		my @tmp=split(/\t|:/,$_);
		my $atlas_string="$tmp[0]\t$tmp[1]\t$tmp[3]\t$tmp[4]";
		$atlas_VR{$atlas_string}=$tmp[14];
		$atlas_RR{$atlas_string}=$tmp[15];
	}
	close(FA);
}