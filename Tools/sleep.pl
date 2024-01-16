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
	"help"			=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

my $ct = 501;
my $tt = 0;
for my $file (@input){
	while (! -e $file){
		sleep $ct;
		$tt += $ct;
		last if $tt > 259200;
	}
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: sleep.pl
#
#        USAGE: ./sleep.pl  
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
#      CREATED: 09/28/17 09:19:14
#     REVISION: ---
#===============================================================================
EOF!
}



