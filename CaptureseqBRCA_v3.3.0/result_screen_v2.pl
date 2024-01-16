#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
#Modified:yangrutao
my $version="2.0";
#die "perl $0 <in.level> <*indel.annovar.BRCA_*.avinput> <atlas.snp.vcf> <atlas.indel.vcf> <result.xls> <panel N> " unless @ARGV==6;
my $panel_n;
my ($level,$avinput,$atlas_snp,$atlas_indel,$out,);
my ($Pdepth,$Pratio,$Vdepth,$Vratio,$Bdepth,$Bratio,);
GetOptions(
				"help|?" =>\&USAGE,
				"l:s"=>\$level,
                                "a:s"=>\$avinput,
                                "as:s"=>\$atlas_snp,
                                "ai:s"=>\$atlas_indel,
				"n:i"=>\$panel_n,
				"o:s"=>\$out,
				"Pdp:i"=>\$Pdepth,
				"Pr:i"=>\$Pratio,
				"Vdp:i"=>\$Vdepth,
				"Vr:i"=>\$Vratio,
				"Bdp:i"=>\$Bdepth,
				"Br:i"=>\$Bratio,
				) or &USAGE;
&USAGE unless ($level && $avinput && $out && $panel_n);
#`sort -k6 -rn -o $ARGV[0] $ARGV[0]`;
`sort -k6 -rn -o $level $level`;
my @gene_names;
my $info="";
my $pmid="";
my @result;
my $Result;
my %hgvsc;
my %transcript;
my %hs_hgvsp;
my %Alleleid;
my %hs_hgvsp_simp;
my %hs_start;
my %hs_ref;
my %hs_alt;
my %aa = ();  
my %aa123=();
my %atlas_VR;
my %atlas_RR;
my $panel_genes = "/share/public/database/Gynecological_cancer_backup/panel_genes/" . $panel_n . "panel_genename.cfg"; 
my $clinvar_hgvs="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/Allgenes_Clinvar20180326.list";
my $clinvar_summary="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/variant_summary_cita_GRCh37_NM_20180326.txt";
my $main_transcript="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/variant_summary.GRCh37.gene_MainNM_ClinVar20171031.txt"; 
my $hgvs="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/HGVS_NM_NP_clinvar20170202.xls";
my $aa_file="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/aa.list.csv" ;
load_panel_gene_file($panel_genes,@gene_names);
load_hgvsc_from_file($clinvar_hgvs, \%hgvsc);
load_clinvar_summary_file($clinvar_summary,\%Alleleid);
load_aa_from_file($aa_file,\%aa,\%aa123);
load_main_transcript_file($main_transcript,\%transcript);
load_hgvs_file($hgvs,\%hs_hgvsp,\%hs_hgvsp_simp);
#load_atlas_file($ARGV[2],\%atlas_VR,\%atlas_RR);
#load_atlas_file($ARGV[3],\%atlas_VR,\%atlas_RR);
load_atlas_file($atlas_snp,\%atlas_VR,\%atlas_RR);
load_atlas_file($atlas_indel,\%atlas_VR,\%atlas_RR);
#print Dumper(%Alleleid);exit;
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

#open IN2,"<$ARGV[1]";
open IN2,"<$avinput" or die $!;
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

#open(IN,"<$ARGV[0]") or die $!;
open(IN,"<$level") or die $!;
while(<IN>)
{
	chomp;
	next if ($_ =~ /^#/);
	next if ($_ =~ /^Level/);
	next if ($_ =~ /^Uncertain significance\tPVS0PS0PM1PP0BA0BS0BP0_pm2\tSampleName/);
	next if ($_ =~ /^Benign\tPVS0PS0PM0PP0BA1BS0BP0_ba1\tSampleName/);
	my @temp=split(/\t/,$_,);
	my $clisig = $temp[0];
	my $gene = $temp[14];
	my $Func_refGene = $temp[13];
	my ($chr,$start,$end,$ref,$alt,$genedetail,$ExonicFunc,)=($temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$temp[15],$temp[16],);
        my ($ref_depth,$alt_depth,$alt_ratio)=($temp[9],$temp[10],$temp[11],);
        if ($ExonicFunc !~/SNV/){
                my $string="$chr\t$start\t$ref\t$alt";
                if (defined $hs_start{$string}){
                        $start = $hs_start{$string};
                        $ref = $hs_ref{$string};
                        $alt = $hs_alt{$string};
                }
                my $string2="$chr\t$start\t$ref\t$alt";
                if (defined $atlas_VR{$string2}  &&  defined $atlas_RR{$string2}){
                        $ref_depth = $atlas_RR{$string2};
                        $alt_depth = $atlas_VR{$string2};
                        $alt_ratio = $alt_depth / ($ref_depth + $alt_depth);
                        $alt_ratio = sprintf "%.2f",$alt_ratio;
                }
        }
        elsif ($ExonicFunc =~/SNV/){
                my $string="$chr\t$start\t$ref\t$alt";
                if (defined $atlas_VR{$string}  &&  defined $atlas_RR{$string}){
                        $ref_depth = $atlas_RR{$string};
                        $alt_depth = $atlas_VR{$string};
                        $alt_ratio = $alt_depth / ($ref_depth + $alt_depth);
                        $alt_ratio = sprintf "%.2f",$alt_ratio;
                }
        }
        #if ((grep /^$gene$/,@gene_names) && (not ($Func_refGene =~ /intronic/  && $temp[17] eq "." )) && (not ($Func_refGene =~ /intronic/ && $temp[17] !~ /NM_/)) && $Func_refGene !~ /UTR/  &&  $ref_depth+$alt_depth >30  &&  $alt_ratio >0.2){    #capture
        if ((grep /^$gene$/,@gene_names) && not ($Func_refGene =~ /intronic/  && $temp[17] eq "." ) && not ($Func_refGene =~ /intronic/ && $temp[17] !~ /NM_/) && $Func_refGene !~ /UTR/ ){ 
	#if ((grep /^$gene$/,@gene_names) && (not ($Func_refGene =~ /intronic/  && $temp[17] eq "." )) && (not ($Func_refGene =~ /intronic/ && $temp[17] !~ /NM_/)) && $Func_refGene !~ /UTR/ ){    #&&  $ref_depth+$alt_depth >50  &&  $alt_ratio >0.2
		if ($clisig =~ /Pathogenic/i && $ref_depth+$alt_depth >$Pdepth  &&  $alt_ratio >$Pratio){ 
                        unless (grep /$clisig/, @result){
                        	push (@result,$clisig);
                        }
                }elsif($clisig =~ /Uncertain/ && $ref_depth+$alt_depth >$Vdepth  &&  $alt_ratio >$Vratio){ 
                        unless (grep /$clisig/, @result){
                        	push (@result,$clisig);
                        }
                }elsif($clisig =~ /Benign/i && $ref_depth+$alt_depth >$Bdepth  &&  $alt_ratio >$Bratio){ 
                        unless (grep /$clisig/, @result){
                        	push (@result,$clisig);
                        }
                }
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

#open IN1,"<$ARGV[0]";
#open OUT,">$ARGV[4]";
open IN1,"<$level";
open OUT,">$out";
print OUT "#Sample\tResult\tChr\tPos\tdbSNP\tRef\tAlt\tGene\tDepth_ref\tDepth_alt\tAlt_ratio\tTranscript\tExon\tHGVS(c.)\tNucleoprotein\tHGVS(p.)\tZygosity\tFunc\.refGene\tType\tName\tCliSig\tEvidence\tPMID\tAAchange\tChr\tStart\tEnd\tCytoBand\tMark\n";

while(<IN1>){
	chomp;
	next if ($_ =~ /^Level/);
	next if ($_ =~ /^Uncertain significance\tPVS0PS0PM1PP0BA0BS0BP0_pm2\tSampleName/);
	next if ($_ =~ /^Benign\tPVS0PS0PM0PP0BA1BS0BP0_ba1\tSampleName/);
	my @atm=split /\t/,$_;
        my $clisig = $atm[0];
	my ($function,$gene,$chr,$start,$end,$genedetail,$ExonicFunc,$AAchange,$genotype,$frequence,$clinvar,$level,$evidence)=($atm[13],$atm[14],$atm[3],$atm[4],$atm[5],$atm[15],$atm[16],$atm[17],$atm[8],$atm[28],$atm[26],$atm[0],$atm[1]);
	my ($sample,$trscript,$exon,$hc,$hp,$ref,$alt,$ref_depth,$alt_depth,$alt_ratio)=($atm[2],$atm[130],$atm[131],$atm[132],$atm[133],$atm[6],$atm[7],$atm[9],$atm[10],$atm[11],);
        my ($cytoBand,$dbsnp,$mark)=($atm[24],$atm[25],$atm[134]);
	if ($gene =~ /^MRE11A$/){
		$gene= "MRE11";
	}
	if ($atm[17] ne "." && $atm[132] eq "." && $atm[133] eq "."){
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
	if ($atm[132] =~ /^c.([ACGT])(\d+)([ACGTYRN])$/) {
            $atm[132] ="c\.$2$1>$3";
    }
	my $hgvsp;
	if(defined $hp && $hp ne "."){
     	if($hp=~/p.[a-zA-Z]{1}\d/)
     	{
     	$hgvsp=hgvp_aa1to3($hp);
		}
	}else{
		$hgvsp = ".";
	}
	my $hc_name=$transcript{$gene}{$trscript} . ":" . $atm[132];
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
		$name = $transcript{$gene}{$trscript} . "\($gene\)" . ":". $atm[132] . "\($hgvsp\)";
	}else{
		$name = $transcript{$gene}{$trscript} . ":" . "\($gene\)" . $atm[132] ;
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
                my $string2="$chr\t$start\t$ref\t$alt";
		if (defined $atlas_VR{$string2}  &&  defined $atlas_RR{$string2}){
			$ref_depth = $atlas_RR{$string2};
			$alt_depth = $atlas_VR{$string2};
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
        if ($mark !~ /clinvar/){
                if (defined $Alleleid{$transcript{$gene}{$trscript}}{$gene}{$hc}){
                        $mark = "clinvar " . $Alleleid{$transcript{$gene}{$trscript}}{$gene}{$hc} . ";$mark";
                }
        } 
	if($clisig =~/Pathogenic/i && $function !~/UTR/ &&  (not ($function =~ /intronic/ && $atm[17] eq ".")) && (not ($function =~ /intronic/ && $atm[17] !~ /NM_/)) && $ref_depth+$alt_depth >$Pdepth  &&  $alt_ratio >$Pratio){    #capture 
        #if($function !~/UTR/ &&  (not ($function =~ /intronic/ && $atm[17] eq ".")) && (not ($function =~ /intronic/ && $atm[17] !~ /NM_/))){    #$ref_depth+$alt_depth >50  &&  $alt_ratio >0.2
		if ($name =~ /MRE11A/){
			$name =~ s/MRE11A/MRE11/g;
		}
		print OUT "$atm[2]\t$Result\t$chr\t$start\t$dbsnp\t$ref\t$alt\t$gene\t$ref_depth\t$alt_depth\t$alt_ratio\t$transcript{$gene}{$trscript}\t$exon\t$hc\t$hgvs_p\t$hgvsp\t$genotype\t$function\t$ExonicFunc\t$name\t$level\t$evidence\t$pmid\t$AAchange\t$chr\t$start\t$end\t$cytoBand\t$mark\n";
        }elsif($clisig =~/Uncertain/ && $function !~/UTR/ &&  (not ($function =~ /intronic/ && $atm[17] eq ".")) && (not ($function =~ /intronic/ && $atm[17] !~ /NM_/)) && $ref_depth+$alt_depth >$Vdepth  &&  $alt_ratio >$Vratio){    #capture 
		if ($name =~ /MRE11A/){
			$name =~ s/MRE11A/MRE11/g;
		}
		print OUT "$atm[2]\t$Result\t$chr\t$start\t$dbsnp\t$ref\t$alt\t$gene\t$ref_depth\t$alt_depth\t$alt_ratio\t$transcript{$gene}{$trscript}\t$exon\t$hc\t$hgvs_p\t$hgvsp\t$genotype\t$function\t$ExonicFunc\t$name\t$level\t$evidence\t$pmid\t$AAchange\t$chr\t$start\t$end\t$cytoBand\t$mark\n";
        }elsif($clisig =~/Benign/i && $function !~/UTR/ &&  (not ($function =~ /intronic/ && $atm[17] eq ".")) && (not ($function =~ /intronic/ && $atm[17] !~ /NM_/)) && $ref_depth+$alt_depth >$Bdepth  &&  $alt_ratio >$Bratio){    #capture 
		if ($name =~ /MRE11A/){
			$name =~ s/MRE11A/MRE11/g;
		}
		print OUT "$atm[2]\t$Result\t$chr\t$start\t$dbsnp\t$ref\t$alt\t$gene\t$ref_depth\t$alt_depth\t$alt_ratio\t$transcript{$gene}{$trscript}\t$exon\t$hc\t$hgvs_p\t$hgvsp\t$genotype\t$function\t$ExonicFunc\t$name\t$level\t$evidence\t$pmid\t$AAchange\t$chr\t$start\t$end\t$cytoBand\t$mark\n";
        }
    
}
close IN1;
#`sort -k6 -n -o $ARGV[0] $ARGV[0]`;
`sort -k6 -n -o $level $level`;
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
		if ($temp[3] ne "-"){
			$hgvsc->{$temp[0]}{$temp[1]}{"pmid"}=$temp[3];
		}
	}
	close FAA;
}
sub load_clinvar_summary_file{
	my ($fclinv, $allelid) = @_;
	open(CL,$fclinv);
	my $line=<CL>;
	while($line=<CL>)
	{
		chomp($line);
		my @tmp=split(/\t/,$line);
                my ($trans,$gen,$hgvsc);
                if ($tmp[2] =~ /(NM\_\d+\.\d+)\((.*?)\):(c.*?)\s/){
                        $trans = $1;
                        $gen = $2;
                        $hgvsc = $3; 
                }elsif ($tmp[2] =~ /(NM\_\d+\.\d+)\((.*?)\):(c.*?)$/){
                        $trans = $1;
                        $gen = $2;
                        $hgvsc = $3; 
                }
		$allelid->{$trans}{$gen}{$hgvsc} = $tmp[0];
	}
	close CL;
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
sub load_panel_gene_file{
	my ($cfg_55, $gene_name) = @_;
	open(PAN,$cfg_55);
	while(<PAN>)
	{
		chomp;
		next if ($_ =~/^\s+$/);
		unless (grep /<$_>/ ,@gene_names){
			push (@gene_names,$_ );
		}
	}
	close(PAN);
}

sub USAGE {
	my $usage=<<"USAGE";
Program: result_screen_v2.pl
Version: $version
Contact: lizhifeng\@berrygenomics.com
Modified:yangrutao  

Description:
Suitable for whole exome sequencing.
Output result.xls

Usage:
  Options:
   -l     level file
   -a     *indel.annovar.BRCA_*.avinput
   -as    atlas.snp.vcf
   -ai    atlas.indel.vcf
   -o     result.xls
   -Pdp   P/LP variant depth(ref+alt)
   -Pr    P/LP variant alt_ratio
   -Vdp   VUS variant depth(ref+alt)
   -Vr    VUS variant alt_ratio
   -Bdp   B/LB variant depth(ref+alt)
   -Br    B/LB variant alt_ratio
   -h     Help

USAGE
	print $usage;
	exit;
}