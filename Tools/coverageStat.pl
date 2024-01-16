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

my (@input , $help);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"help"	=>	\$help,
);

if ($help or $#input<0){
	&help;
	exit;
}

my @cut = qw/100 150 200 250 300/;
for my $file (@input){
	my ($sample) = ($file =~ /\/([\w-]+)\/depth_frequency.xls/);
	print $sample;
	open IN , "$file";
	my $c = 0;
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		for my $cut (@cut){
			if ($cut == $F[0]){
				printf("\t%.2f" , (1-$c)*100);
				print "%";
			}
		}
		$c += $F[1];
		last if $F[0] >= $cut[-1];
	}
	close IN;
	print "\n";
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: coverageStat.pl
#
#        USAGE: ./coverageStat.pl  
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
#      CREATED: 09/17/2018 01:44:59 PM
#     REVISION: ---
#===============================================================================
EOF!
}



