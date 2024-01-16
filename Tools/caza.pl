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

my ($input , $help , $normal);
GetOptions(
	"i|input=s"	=>	\$input,
	"n|normal=s"	=>	\$normal,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %opt;
open IN , "$normal";
#1       7000001 8000000 0.479461        0.981314        2.00099443976422        0.0330413041554388      2.04528866523792
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$opt{"$F[0]\t$F[1]\t$F[2]"} = [@F[5,6]];
}
close IN;

open IN , "$input";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	print join("\t" , @F[1..3,5..7]) , "\t" , $opt{"$F[1]\t$F[2]\t$F[3]"}->[0] , "\t" , $opt{"$F[1]\t$F[2]\t$F[3]"}->[1] , "\t";
	if ($F[11] eq 'NA' or $opt{"$F[1]\t$F[2]\t$F[3]"}->[0] eq 'NA'){
		$F[11] = $F[11]*2 if $F[11] ne 'NA';
		print "$F[11]\tNA\n";
	}else{
		$F[11] = $F[11]*2;
		print ("$F[11]\t" , ($F[11]-$opt{"$F[1]\t$F[2]\t$F[3]"}->[0])/$opt{"$F[1]\t$F[2]\t$F[3]"}->[1] , "\n");
	}
}
close IN;


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: caza.pl
#
#        USAGE: ./caza.pl  
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
#      CREATED: 03/08/18 13:56:07
#     REVISION: ---
#===============================================================================
EOF!
}



