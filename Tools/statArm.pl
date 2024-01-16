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

my $cut = 3;
open IN , "$input";
my $h = <IN>;
chomp $h;
my @h = split "\t" , $h;
print "score\tnumber\t$h\n";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my $n = 0;
	my $t = 0;
	for my $i (1..$#F){
		my $d = $F[$i];
		next unless $h[$i] =~ /\d/;
		$n++ if abs($d) > $cut;
		$t += abs($d);
	}
	print "$t\t$n\t$_\n";
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: statArm.pl
#
#        USAGE: ./statArm.pl  
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
#      CREATED: 08/09/2018 02:41:10 PM
#     REVISION: ---
#===============================================================================
EOF!
}



