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

open IN , "$input";
my $head = <IN>;
print $head;
while (<IN>){
	my @F = split /\t/ , $_;
	if ($F[8] >= 4 and $F[9] >= 0.01 and ($F[43] eq '.' or $F[43] < 0.05) and $F[84] eq '.'){
		print $_;
	}
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: filterNonsynonymous.pl
#
#        USAGE: ./filterNonsynonymous.pl  
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
#      CREATED: 04/13/18 10:28:31
#     REVISION: ---
#===============================================================================
EOF!
}



