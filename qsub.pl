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

my (@input , $help , $shell , $l , $queue);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"s|shell=s"	=>	\$shell,
	"l|l=s"	=>	\$l,
	"q|queue=s"	=>	\$queue,
	"help"			=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

for my $m (@input){
	`rm -f $m`;
}

my $dir = dirname $shell;
chdir $dir;
`rm -f $dir/ERROR`;
my $qsub = readpipe("qsub -cwd -V -l $l -q $queue $shell") or exit 1;
print $qsub;

my $ct = 501;
my $tt = 0;

sleep 1;
for my $file (@input){
	while (! -e $file){
		sleep $ct;
		$tt += $ct;
		#last if $tt > 259200;
		exit 1 if -e "ERROR";
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



