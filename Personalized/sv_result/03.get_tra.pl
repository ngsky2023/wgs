use strict;
use warnings;
use Cwd;
use Data::Dumper;

my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
# open C,"../sv/chromothripsis/run/class.xls";
# my %tra;
# while(<C>){
	# chomp;
	# next if $.==1;
	# my @line=split/\t/;
	# next if $line[11]<=34;
	# next if $line[3]=~/NA/;
	# next if $line[0]=~/high1/;
	# my $chr=join "\t",@line[3,4];
	# $tra{$line[1]}{$line[2]}=$chr;
# }

open OUT,">tra.xls";
print OUT "sample\tCHROM_A\tSTART_A\tEND_A\tCHROM_B\tSTART_B\tEND_B\tTYPE\n";
foreach my $sample(@sample){
	open F,"SV_result/$sample.somaticSV.xls" or die "$sample is error\n";
	my %id;
	while(<F>){
		chomp;
		next if /^#/;
		my @line=split/\t/,$_;
		my $out=join "\t",@line[0..5,10];
		my $chr1=$line[0];
		$chr1=~s/chr//g;
		my $chr2=$line[3];
		$chr2=~s/chr//g;
		next if $line[10] eq 'INS';
		next if $.==1;
		
		#if ($_=~/SVLEN=(\d+)/){
		#	next if $1<500;
			# print "$1\n";
		#}
		print OUT "$sample\t$out\n";
		# if(defined $tra{$sample}){
			# if(defined $tra{$sample}{$chr1}){
				# my @location=split/\t/,$tra{$sample}{$chr1};
				# if($line[1] >=$location[0] && $line[2]<=$location[1]){
					# print OUT "$sample\t$out\n";
				# }
				
			# }elsif(defined $tra{$sample}{$chr2}){
				# my @location=split/\t/,$tra{$sample}{$chr2};
				# if($line[4] >=$location[0] && $line[5]<=$location[1]){
					# print OUT "$sample\t$out\n";
				# }
			# }
		# }
		
	}
	close F;
}
