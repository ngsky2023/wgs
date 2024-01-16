use strict;
use warnings;

die "perl $0 <in.vcf> <in.sum> <out>" unless @ARGV==3;

open IN1,"<$ARGV[0]";
my %hash;
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	my ($Chr,$Pos,$Ref,$Alt)=($atm[0],$atm[1],$atm[3],$atm[4]);
	my $info=$atm[-1];
	my $Genotype;
	if($info =~ /^(\d\/\d):/ || $info =~ /^(\d\|\d):/){
		$Genotype=$1;
		$Genotype=~ s/\|/\//g;
		if ($Genotype eq '0/1' || $Genotype eq '1/0'){
                                $Genotype = 'het';
                        }else{
                                $Genotype = 'hom';
                        }
	}	
	$hash{$Chr}{$Pos}{$Ref}{$Alt}=$Genotype;
}
close IN1;

my %nm;
open NP,"</share/public/database/Gynecological_cancer_backup/IPDB/ver4/HGVS_NM_NP_clinvar20170202.xls" or die $!;
while (<NP>){
        chomp;
        next  if ($_ =~ /^\s+$/);
        my @dat1=split /\t/,$_;
        if ($dat1[1] =~ /^MRE11A$/) {
                $dat1[1] = "MRE11";
        }
        $nm{$dat1[1]}{$dat1[2]}=$dat1[3];
}
close NP;

open IN2,"<$ARGV[1]";
open OUT,">$ARGV[2]";
print OUT "Gene\tHGVS(c.)\tHGVS(p.)\tGenotype\n";
while(<IN2>){
	chomp;
	my @atm=split /\t/;
	my ($Chr,$Pos,$Ref,$Alt,$Gene,$Level,$Disease,$Type,$c,$g,$Citation)=($atm[0],$atm[1],$atm[3],$atm[4],$atm[5],$atm[9],$atm[7],$atm[17],$atm[14],$atm[15],$atm[46]);
	$Chr="chr".$Chr;
	my $trans;
	if ($c =~ /(NM_.*?)\:c\..*/){
		$trans = $1;
	}
	if(exists $hash{$Chr}{$Pos}{$Ref}{$Alt}){
		print OUT "$Gene\t$c\t$nm{$Gene}{$trans}\:$g\t$hash{$Chr}{$Pos}{$Ref}{$Alt}\n";
	}
}
close IN2;
close OUT;
