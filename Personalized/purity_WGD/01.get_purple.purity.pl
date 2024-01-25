use strict;
use warnings;
use Cwd;
use Data::Dumper;
my %hash;
my $samplelist="sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
chomp(my $head=`cat head`);
open N,">ALL.purity.txt";
print N "name\t$head\n";
open OUT,">01.purity_WGD.txt";
print OUT "sample\tpurity\tgender(bio)\twholeGenomeDuplication\n";
foreach my $sample(@sample){
		next if $sample=~/^#/;
        	my $path="/path/wgs_analysis"
        	my @file_list=glob("$path/*/cnv/$sample/purple*/$sample.purple.purity.tsv");
        	print "$sample\t@file_list\n";
		#my $file="$sample.purple.purity.tsv";
		if(@file_list){
		
			#`ln -s $file ./`;
			open F,"@file_list";
			while(<F>){
				chomp;
				next if /^purity/;
               			my @line=split/\t/,$_;
                		my $out=join "\t",@line[0,5,16];
                		print OUT "$sample\t$out\n";
                		print N "$sample\t$_\n";
			}
		}else{
			print STDERR "$sample is error\n";
		}
		
		
}













