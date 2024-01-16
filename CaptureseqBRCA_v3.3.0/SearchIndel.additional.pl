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
	$hash{$Chr}{$Pos}{$Ref}{$Alt}=$_;
}
close IN1;

open IN2,"<$ARGV[1]";
open OUT,">>$ARGV[2]";
while(<IN2>){
	chomp;
	my @atm=split /\t/;
	#my ($Chr,$Pos,$Ref,$Alt,$Gene,$Level,$Disease,$Type,$c,$g,$Citation)=($atm[0],$atm[1],$atm[2],$atm[3],$atm[4],$atm[6],$atm[5],$atm[14],$atm[11],$atm[12],$atm[9]);  #55genes
	my ($Chr,$Pos,$Ref,$Alt,$Gene,$Level,$Disease,$Type,$c,$g,$Citation)=($atm[0],$atm[1],$atm[3],$atm[4],$atm[5],$atm[9],$atm[7],$atm[17],$atm[14],$atm[15],$atm[46]);  #57genes
	$Chr="chr".$Chr;
	if(exists $hash{$Chr}{$Pos}{$Ref}{$Alt}){
		print OUT "$hash{$Chr}{$Pos}{$Ref}{$Alt}\n";
	}
}
close IN2;
close OUT;
