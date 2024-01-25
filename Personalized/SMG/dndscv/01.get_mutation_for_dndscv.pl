use strict;
use warnings;
use Cwd;
use Data::Dumper;
my %hash;
my $samplelist="../../clinical/sample.list";
my $path="../../snv_indel_result";
open IN,$samplelist;
chomp(my @sample=<IN>);
open OUT,">acral_coding_snv_indel_fordndscv.xls";
print OUT "\#SampleName\tChr\tStart\tRef\tAlt\n";
foreach my $sample(@sample){
    if(-e "$path/$sample.final.xls"){
		open S,"$path/$sample.final.xls";
		while(<S>){
			chomp;
			next if /^SampleName/;
			my @line=split/\t/,$_;
			$line[1]=~s/chr//g;
			# next unless ($line[16]=~/^exonic/|| $line[16]=~/^splicing/);
			my $out=join "\t",@line[0,1,2,4,5];
			print OUT "$out\n";
		}
	}else{
		print "$sample is error\n";
	}
}
