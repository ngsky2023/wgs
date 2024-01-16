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

my ($input , $help , $bed);
GetOptions(
	"i|input=s"	=>	\$input,
	"b|bed=s"	=>	\$bed,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %soft;
open IN , "$Bin/soft.list";
while (<IN>){
	chomp;
	my ($k , $d) = split /\t/ , $_;
	$soft{$k} = $d;
}
close IN;

my %bam;
open IN , "$input";
while (<IN>){
	chomp;
	my ($sample , $lu) = split /\t/ , $_;
	push @{$bam{$sample}} , $lu;
}
close IN;

mkdir 'analysis' unless -d 'analysis';
for my $sample (sort keys %bam){
	mkdir "analysis/$sample" unless -d "analysis/$sample";
	mkdir "analysis/$sample/cut" unless -d "analysis/$sample/cut";
	mkdir "analysis/$sample/2.Realign" unless -d "analysis/$sample/2.Realign";
	my $inm = '';
	for my $bam (@{$bam{$sample}}){
		my $out = basename $bam;
		my $cmd = "$soft{sambamba} view $bam -L $bed -f bam -o analysis/$sample/cut/$out -t 8";
		print "$cmd\n";
		readpipe($cmd);
		$inm .= " analysis/$sample/cut/$out";
	}
	my $cmd = "$soft{sambamba} merge analysis/$sample/2.Realign/$sample.bam $inm -t 8";
	print "$cmd\n";
	readpipe($cmd);
}


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: somaticBam2germ.pl
#
#        USAGE: ./somaticBam2germ.pl  
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
#      CREATED: 11/29/2018 03:19:50 PM
#     REVISION: ---
#===============================================================================
EOF!
}



