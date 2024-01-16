#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib $Bin;

my ($input , $help);
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %opt;
#open IN , "/share/work1/wangrr/DB/hg19/binCor/RDcfDNACNV.ARM.Q60.zscore.txt";
open IN , "/share/work1/wangrr/DB/hg19/binCor/normal.ARM.Q60.P.txt";
#open IN , "/share/work1/wangrr/DB/hg19/binCor/cfDNA.ARM.Q60.P.txt";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$opt{"$F[0]"} = [@F[1,2]];
}
close IN;

my %arm;
my $total = 0;
open IN , "/share/work1/wangrr/DB/hg19/centromere.bed";
<IN>;
#chrY    10104553        13104553        59373566
#chrM    0       1       16571
while (<IN>){
	chomp;
	my ($chr , $str , $end , $len) = split /\t/ , $_;
	$arm{$chr} = [$str , $end , $len];
}
close IN;

open IN , "$input";
#fixedStep chrom=chr1 start=1 step=1000000 span=1000000
#725970
my %loc;
my ($win , $chr , $ch);
my %count;
while (<IN>){
	chomp;
	if (/fixedStep chrom=(\S+) start=\d+ step=(\d+) span=\d+/){
		($chr , $win) = ($1 , $2);
		$ch = $chr;
		$ch =~ s/chr//;
	}else{
		next unless exists $arm{$chr};
		$total += $_;
		$loc{$chr} += $win;
		if ($loc{$chr} < $arm{$chr}->[0]){
			my $am = $ch;
			if ($am ne 'chrM'){
				$am = $am."p";
			}
			$count{$am} += $_;
		}elsif ($loc{$chr} > $arm{$chr}->[1]){
			my $am = $ch;
			if ($am ne 'M'){
				$am = $am."q";
			}
			$count{$am} += $_;
		}
	}
}
close IN;

#print "ARM\tP\tZscore\tmean\tSD\n";
print "ARM\tZscore\n";
for my $am (sort keys %count){
	next unless exists $opt{$am} and exists $count{$am};
	my ($mean , $sd) = @{$opt{$am}};
	my $P = $count{$am}/$total;
#	print "$am\t$P\t" , ($P-$mean)/$sd , "\t$mean\t$sd\n";
	print "$am\t" , ($P-$mean)/$sd , "\n";
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeWindow2Arm.pl
#
#        USAGE: ./mergeWindow2Arm.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wang Ruiru (wangrr), wangruiru\@berrygenomics.com
# ORGANIZATION: Berry Genomics
#      VERSION: 1.0
#      CREATED: 04/03/18 14:16:43
#     REVISION: ---
#===============================================================================
EOF!
}



