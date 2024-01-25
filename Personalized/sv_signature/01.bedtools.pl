use strict;
use warnings;
use Cwd;
use Data::Dumper;
my $samplelist="../clinical/sample.list";
my $out_path="./bedtools/";
open IN,$samplelist;
chomp(my @sample=<IN>);
close IN;
my $svpath="../sv_result/SV_result";
foreach my $sample(@sample){
	open OUT,">$out_path/$sample.xls";
	open F,"$svpath/$sample.somaticSV.xls";
	my $num;
	while(<F>){
		chomp;
		next if /^#|^\s+$/;
		my @line=split/\t/,$_;
		next if $line[10] eq "INS";
		if($line[10] eq "BND" || $line[12]=~/BND/){
			# my $out=join "\t",@line[0,1,5];
			#my $end=$line[1]+1000000;
			print OUT "$line[0]\t$line[1]\t$line[2]\tTRA\tTRA\n$line[3]\t$line[4]\t$line[5]\tTRA\tTRA\n";
		}else{
			$line[12]=~s/.*;SVLEN=(.*?);.*/$1/g;
			$line[12]=abs($line[12]);
			my $out1=join "\t",@line[0,1,2,10,12];
			my $out2=join "\t",@line[3,4,5,10,12];
			print OUT "$out1\n$out2\n";
		}
		
		
		
	}
}
