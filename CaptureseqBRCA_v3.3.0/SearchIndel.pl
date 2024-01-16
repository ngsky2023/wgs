use strict;
use warnings;

die "perl $0 <in.annotate> <in.indel> <in.sum> <out.prefix>" unless @ARGV==4;

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
open OUT,">$ARGV[3].AlleleID";
while(<IN2>){
	my $line1=$_;
	my $line2=<IN2>;
	chomp($line1,$line2);
	my ($Start1,$End1,$Ref1,$Alt1,$Gene,$AlleleID,$Level,$Disease,$Type,$c,$g,$Citation)=(split /\t/,$line1)[1,2,3,4,5,6,7,8,9,14,15,17];
	my ($Start2,$End2,$Ref2,$Alt2)=(split /\t/,$line2)[1,2,3,4];
	if($Level eq "Red"){
		$Level = "Pathogenic";
	}else{
		$Level = "Likely pathogenic";
	}
	if(exists $hash{$Start1}{$End1}{$Ref1}{$Alt1} && exists $hash{$Start2}{$End2}{$Ref2}{$Alt2}){
		print OUT "$AlleleID\t$hash{$Start1}{$End1}{$Ref1}{$Alt1}->[0]\n";
	}
}
close IN2;
close OUT;

foreach my $key1(keys %hash){
	foreach my $key2(keys %{$hash{$key1}}){
		foreach my $key3(keys %{$hash{$key1}{$key2}}){
			if($key1 eq $key2  && $key2 eq "133634057" && $key3 eq "T"){
				`sed -i s/25094/-/g $ARGV[3].AlleleID`;
			}
		}
	}	
}

if(-s "$ARGV[3].AlleleID"){
	open IN3,"<$ARGV[3].AlleleID";
	my %hashd;
	while(<IN3>){
		chomp;
		my @tmp=split /\t/;
		$hashd{$tmp[0]}=$tmp[1];
	}
	close IN3;
	open IN4,"<$ARGV[2]";
	open OUT2,">$ARGV[3].info";
	#print OUT2 "Gene\tHGVS(c.)\tHGVS(p.)\tGenotype\tPosition\tType\tDisease\tLevel\tCitation\n";
	my %hashs;
	while(<IN4>){
		chomp;
		my @atm=split /\t/;
		my ($Chr,$Start,$End,$Ref,$Alt,$Gene,$Level,$Disease,$Type,$c,$p,$Citation)=($atm[0],$atm[1],$atm[2],$atm[3],$atm[4],$atm[5],$atm[7],$atm[8],$atm[9],$atm[14],$atm[15],$atm[17]);
		if(exists $hashd{$atm[6]}){
			print OUT2 "$Gene\t$c\t$p\t$hashd{$atm[6]}\tchrX:$Start-$End\t$Type\t$Disease\t$Level\t$Citation\n";
		}
	}
	close IN4;
	close OUT2;
}else{
	exit;
}
