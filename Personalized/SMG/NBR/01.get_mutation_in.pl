use strict;
use warnings;
use Cwd;
use Data::Dumper;
my %hash;
my $samplelist="../../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
open OUT,">mutations.txt";
print OUT "sampleID\tchr\tpos\tref\tmut\n";
foreach my $sample(@sample){
		next if $sample=~/^#/;
		my $file="../../snv_indel_result/$sample.final.xls";
		open F,$file;
		#open OUT,">/share/work1/hanwj4457/project/kongyan_jiayan_zhiduan+nianmo/acral/noncoding/IN/$sample.snv.txt";
		#print OUT "sampleID\tchr\tpos\tref\tmut\n";
		while(<F>){
			chomp;
			my @line=split/\t/,$_;
			next if /^SampleName/;
			$line[1]=~s/chr//g;
			print OUT "$sample\t$line[1]\t$line[2]\t$line[4]\t$line[5]\n";
		}
		close F;
}













