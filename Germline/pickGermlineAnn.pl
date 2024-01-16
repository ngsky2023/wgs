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

my (@input , $help , $vcf);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"v|vcf=s"	=>	\$vcf,
	"help"	=>	\$help,
);

if ($help){
	&help;
	exit;
}

my %pos;
open IN , "$vcf";
while (<IN>){
	next if /^#/;
	chomp;
	my @F = split /\t/ , $_;
	$pos{"$F[0]\t$F[1]"} = 1;
}
close IN;

my $i = 0;
for my $file (@input){
	open IN , "$file";
	my $head = <IN>;
	print $head if $i == 0;
	while (<IN>){
		my @F = split /\t/ , $_;
		if ($F[5] eq '-'){
			$F[2]--;
		}
		print $_ if exists $pos{"$F[1]\t$F[2]"};
	}
	close IN;
	$i++;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: pickGermlineAnn.pl
#
#        USAGE: ./pickGermlineAnn.pl  
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
#      CREATED: 06/12/18 09:33:11
#     REVISION: ---
#===============================================================================
EOF!
}



