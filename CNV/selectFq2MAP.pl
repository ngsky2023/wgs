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

my ($input , $help , $dir);
GetOptions(
	"i|input=s"	=>	\$input,
	"d|dir=s"	=>	\$dir,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %fq;
open IN , "$input";
while (<IN>){
	chomp;
	my ($sample , $r1) = split /\t/ , $_;
	my @stat = stat($r1);
	my $size = $stat[7];
	if (not exists $fq{$sample} or $fq{$sample}->[1] < $size){
		$fq{$sample} = [$r1 , $size];
	}
}
close IN;

mkdir $dir unless -d $dir;
for my $f (values %fq){
	my $file = $f->[0];
	readpipe("ln -s $file $dir/");
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: selectFq2MAP.pl
#
#        USAGE: ./selectFq2MAP.pl  
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
#      CREATED: 08/07/2018 05:17:23 PM
#     REVISION: ---
#===============================================================================
EOF!
}



