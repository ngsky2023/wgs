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
my $cutoff = 3;
GetOptions(
	"i|input=s"	=>	\$input,
	"c|cutoff=s"	=>	\$cutoff,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %score;
open IN , "$Bin/CharmScore.txt";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$score{$F[0]} = [@F[1..5]];
}
close IN;

open IN , "$input";
my $h = <IN>;
my @h = split /\t/ , $h;
shift @h;
while (<IN>){
	chomp;
	my ($sample , @F) = split /\t/ , $_;
	my ($del , $dup) = (0 , 0);
	for my $i (0..$#h){
		next unless exists $score{$h[$i]};
		if ($F[$i] > $cutoff){
			$dup += $score{$h[$i]}->[1];
		}elsif ($F[$i] < -1 * $cutoff){
			$del += $score{$h[$i]}->[0];
		}
	}
	print "$sample\t$del\t$dup\n";
}
close IN;


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: CharmScore.pl
#
#        USAGE: ./CharmScore.pl  
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
#      CREATED: 07/25/18 17:08:51
#     REVISION: ---
#===============================================================================
EOF!
}



