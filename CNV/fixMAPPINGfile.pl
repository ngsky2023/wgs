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

my $um;
open IN , "$input";
<IN>;
print "chr\tuniq\tuniqGC%\tuniqGC\tp_uniq\tp_uniqGC%\tp_uniqGC\tq_uniq\tq_uniqGC%\tq_uniqGC\n";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	print join("\t" , @F[0,1,3,2,4,6,5,7,9,8]) , "\n";
	$um += $F[1];
}
close IN;
my $map = int($um * 1.1);
my $total = $um * 3;
print ">no_N_Reads\t$total\n";
print ">totalMapped\t$map\n";
print ">uniqMapped\t$um\n";
printf(">redundancy\t%.6f\n" , ($map-$um)/$map);

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: fixMAPPINGfile.pl
#
#        USAGE: ./fixMAPPINGfile.pl  
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
#      CREATED: 03/21/18 10:41:33
#     REVISION: ---
#===============================================================================
EOF!
}



