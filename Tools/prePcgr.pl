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

if ($help or $#input < 0){
	&help;
	exit;
}

my $h = 1;
for my $file (@input){
	if ($file =~ /\.gz$/){
		open IN , "zcat $file |";
	}else{
		open IN , "$file";
	}
	while (<IN>){
		if (/^#/){
			print $_ if $h;
			next;
		}
		$h = 0;
		my @F = split /\t/ , $_;
		$F[6] = 'PASS';
		print join("\t" , @F);
	}
	close IN;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: prePcgr.pl
#
#        USAGE: ./prePcgr.pl  
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
#      CREATED: 09/30/17 16:49:32
#     REVISION: ---
#===============================================================================
EOF!
}



