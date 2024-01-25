use strict;
use warnings;
use Cwd;
use Data::Dumper;
my %hash;
readpipe("mkdir snv");
my $snv_path="../snv_indel_result";
my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
foreach my $sample(@sample){
open OUT,">snv/$sample.vcf";
	print OUT "\#CHROM  POS     ID      REF     ALT\n";
        my (@title,@snv);
	my $num=0;
        open F,"$snv_path/$sample.final.xls" or die "$sample file is error \n";
        while(<F>){
                chomp;
                next if /^\s+$/;
		next if /^SampleName/;
		my @line=split/\t/,$_;
		print OUT "$line[1]\t$line[2]\t$line[0]\t$line[4]\t$line[5]\n";
		
        }
}













