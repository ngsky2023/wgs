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
my $prx = './out';
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
	"prx=s"	=>	\$prx,
);

if ($help or ! $input){
	&help;
	exit;
}

my @chr = qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM/;
my %fh;
for my $chr (@chr){
	open $fh{$chr} , ">$prx.$chr.vcf";
}

if ($input =~ /\.gz$/){
	open IN , "zcat $input|";
}else{
	open IN , "$input";
}
while (<IN>){
	if (/^#/){
		for my $chr (@chr){
			print {$fh{$chr}}  $_;
		}
		next;
	}
	my ($chr) = (/^(\S+)/) or next;
	print {$fh{$chr}}  $_;
}
close IN;

for my $chr (@chr){
	close $fh{$chr};
}



sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: cutVcf.pl
#
#        USAGE: ./cutVcf.pl  
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
#      CREATED: 10/25/17 14:34:41
#     REVISION: ---
#===============================================================================
EOF!
}



