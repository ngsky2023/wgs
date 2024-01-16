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

#HWI-D00477:789:H3YVYBCX2:2:1102:11430:9621      99      chr1    2489746 60      5M1I144M        =       2489766 170     C
open IN , "$input";
while (<IN>){
	if (/^@/){
		print $_;
		next;
	}
	chomp;
	my $nm = 0;
	my @F = split /\t/ , $_;
	for my $bb (@F[11..$#F]){
		if ($bb =~ /NM:i:(\d+)/){
			$nm = $1;
			last;
		}
	}
	my ($len , $l , $s) = mode($F[5]);
	$nm -= $l;
	$nm += $s;
	print "$_\n" if $len > 0 and $nm/$len <= 0.04;
}
close IN;

sub mode{
	my $m = $_[0];
	my ($len , $l , $s) = (0 , 0 , 0);
	while ($m =~ /(\d+)([MID])/g){
		if ($2 eq 'M'){
			$len += $1;
		}elsif ($2 eq 'I'){
			$len += $1;
			$l += $1-1;
		}elsif ($2 eq 'D'){
			$l += $1-1;
		}elsif ($2 eq 'S'){
			$s++;
		}
	}
	return ($len , $l , $s);
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: filterMem.pl
#
#        USAGE: ./filterMem.pl  
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
#      CREATED: 01/09/18 14:59:54
#     REVISION: ---
#===============================================================================
EOF!
}



