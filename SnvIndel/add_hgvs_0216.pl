use strict;
use warnings;

die "perl $0 <in1.key> <in2.infor> <out>" unless @ARGV==3;

my $lst=shift;
my $total=shift;
my $out=shift;

my %hash;

open IN1,"<$lst";
while(<IN1>){
	chomp;
	my @atm=split /\t/;
	my $id="rs".$atm[0];
	$hash{$id}=[$atm[1],$atm[2]];
}
close IN1;

open IN2, "<$total";
open OUT, ">$out";
my $title=<IN2>;
print OUT "$title";
while(<IN2>){
	chomp;
	my @tmp=split /\t/;
	my $num=$#tmp;
	my $gene=$tmp[6];
	my ($detail,$AAchange,$snp)=($tmp[7],$tmp[9],$tmp[17]);
	if($AAchange eq "."){
		if(exists $hash{$snp}){
			$AAchange="$hash{$snp}->[0]:$hash{$snp}->[1]";
		}elsif($detail ne "." && $detail !~ /dist/){
			$AAchange=$detail;
		}
	}
	print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[8]\t$AAchange\t";
	for(my $i =10;$i<$#tmp;$i++){
		print OUT "$tmp[$i]\t";
	}
	print OUT "$tmp[$num]\n";
}
close IN2;
close OUT;
