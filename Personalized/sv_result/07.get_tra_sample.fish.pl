use strict;
use warnings;
use Cwd;
use Data::Dumper;
use List::MoreUtils qw(uniq);


my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
my %hash;
my @chr1=(1..22,"X","Y");
my @chr2=(1..22,"X","Y");
# my @location;
foreach my $num1(0..$#chr1){
	my $num3=$num1;
	foreach my $num2($num3..$#chr2){
		my $chr1=$chr1[$num1];
		my $chr2=$chr1[$num2];
		# next if $chr1 eq $chr2;
		# print "$chr1\t$chr2\n";
		$hash{"chr$chr1;chr$chr2"}="";
	}
	
}
# print Dumper \%hash;
# exit;
open OUT,">tra_chr_join_samplecount.xls";
my $head=join "\tchr",@chr1;
print OUT "\tchr$head\n";;
my %sample;
foreach my $sample(@sample){
	open F,"SV_result/$sample.somaticSV.xls";
	my %id;
	while(<F>){
		chomp;
		next if $.==1;
		my @line=split/\t/,$_;
		next unless $line[10] eq 'BND';
		my $loc=join ";",@line[0,3];
		my $loc2=join ";",@line[3,0];
		my $chr1=$line[0];
		# $chr1=~s/chr//g;
		my $chr2=$line[3];
		# $chr2=~s/chr//g;	
		# if ($_=~/SVLEN=(\d+)/){
			# next if $1<500;
			# print "$1\n";
		# }
		if(defined $hash{$loc}){
			$hash{$loc}.="$sample\t";
			if(defined $sample{$sample}{$loc}){
				$sample{$sample}{$loc}+=1;
			}else{
				$sample{$sample}{$loc}=1;
			}
		}elsif(defined $hash{$loc2}){
			$hash{$loc}.="$sample\t";
			if(defined $sample{$sample}{$loc2}){
				$sample{$sample}{$loc2}+=1;
			}else{
				$sample{$sample}{$loc}=1;
			}
		}else{
			print "Error $loc\t$loc2\n";
			# print STDERR "$sample\t$loc\n";
		}
		
		
	}
	close F;
}
# print Dumper \%hash;
foreach my $chr1(@chr1){
	print OUT "chr$chr1\t";
	foreach my $chr2(@chr2){
		my $len;
		my @s;
		if(exists $hash{"chr$chr1;chr$chr2"}){
			$hash{"chr$chr1;chr$chr2"}||="";
			@s=uniq (split"\t",$hash{"chr$chr1;chr$chr2"});
			$len=@s;
		}elsif(exists $hash{"chr$chr2;chr$chr1"} ){
			$hash{"chr$chr2;chr$chr1"}||="";
			@s=uniq (split"\t",$hash{"chr$chr2;chr$chr1"});
			$len=@s;
		}
		print OUT "$len\t";
		
	}
	print OUT "\n";
}
close OUT;
#每个样本的易位
open OUT2,">tra_each_sample_chr.xls";
print OUT2 "sample\t";
foreach my $num1(0..$#chr1){
	my $num3=$num1;
	foreach my $num2($num3..$#chr2){
		my $chr1=$chr1[$num1];
		my $chr2=$chr1[$num2];
		print OUT2 "chr$chr1;chr$chr2\t";
	}
	
}
print OUT2 "\n";
foreach my $sample(@sample){
	print OUT2 "$sample\t";
	foreach my $num1(0..$#chr1){
		my $num3=$num1;
		foreach my $num2($num3..$#chr2){
			my $chr1="$chr1[$num1]";
			my $chr2="$chr1[$num2]";
			my $num;
			if(exists $sample{$sample}{"chr$chr1;chr$chr2"}){
				$num=$sample{$sample}{"chr$chr1;chr$chr2"};
			}elsif(exists $sample{$sample}{"chr$chr2;chr$chr1"}){
				$num=$sample{$sample}{"chr$chr2;chr$chr1"};
			}
			$num||=0;
			print OUT2 "$num\t";
			
		}
	}
	print OUT2 "\n";
}
