use strict;
use warnings;

die "perl $0 <in.annotate> <in.sum> <out>" unless @ARGV==3;

open IN1,"<$ARGV[0]";
my %hash;
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	my ($Chr,$Start,$End,$Ref,$Alt,$Genotype,$Alt_ratio,$snp)=($atm[0],$atm[1],$atm[2],$atm[3],$atm[4],$atm[5],$atm[8],$atm[58]);
	if($Chr eq "chrX" && $Alt_ratio >0){
		$hash{$Start}{$End}{$Ref}{$Alt}=[$Genotype,$snp];
	}
}
close IN1;

open IN2,"<$ARGV[1]";
open OUT,">$ARGV[2]";
print OUT "Gene\tHGVS(c.)\tHGVS(p.)\tGenotype\tPosition\tType\tDisease\tLevel\tCitation\n";
while(<IN2>){
	chomp;
	my @atm=split /\t/;
	my ($Chr,$Start,$End,$Ref,$Alt,$Gene,$Level,$Disease,$Type,$c,$g,$Citation)=($atm[0],$atm[1],$atm[2],$atm[3],$atm[4],$atm[5],$atm[7],$atm[8],$atm[9],$atm[14],$atm[15],$atm[17]);
	if($Level eq "Red"){
		$Level = "Pathogenic";
	}else{
		$Level = "Likely pathogenic";
	}
	if(exists $hash{$Start}{$End}{$Ref}{$Alt}){
		print OUT "$Gene\t$c\t$g\t$hash{$Start}{$End}{$Ref}{$Alt}->[0]\tChrX:$Start-$End\t$Type\t$Disease\t$Level\t$Citation\n";
	}
}
close IN2;
close OUT;
