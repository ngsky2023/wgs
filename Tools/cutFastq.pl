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

my ($input , $out , $help , $r);
my ($cutN) = (20000000);
GetOptions(
	"i|input=s"		=>	\$input,
	"o|output=s"	=>	\$out,
	"r=s"			=>	\$r,
	"help"			=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}
$cutN *= 4;

my $sb = 'aa';
open IN , "/share/work3/wangrr/local/simple/bin/pigz -d -c -p 4 $input |";
my $file = "$out\_$sb\_$r.fq.gz";
open OUT , "| /share/work3/wangrr/local/simple/bin/pigz -p 4 > $file";
my $n = 0;
while (<IN>){
	$n++;
	print OUT $_;
	if ($n % $cutN == 0){
		close OUT;
		$sb++;
		my $file = "$out\_$sb\_$r.fq.gz";
		open OUT , "| /share/work3/wangrr/local/simple/bin/pigz -p 4 > $file";
	}
}
close OUT;
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: cutFastq.pl
#
#        USAGE: ./cutFastq.pl  
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
#      CREATED: 10/16/17 13:35:03
#     REVISION: ---
#===============================================================================
EOF!
}



