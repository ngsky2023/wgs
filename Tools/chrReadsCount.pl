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

my %count;
my %total;
my %len;
my @q = (0,1,10,20,30,40,50,60);
open IN , "samtools view -h $input|";
while (<IN>){
	chomp;
	if (/^@/){
		if (/^\@SQ\tSN:(\S+)\tLN:(\d+)/){
			$len{$1} = $2;
			#print "$1\t$2\n";
			$count{$1} = {}
		}
		next;
	}
	my ($chr , $q) = (split /\t/ , $_)[2,4];
	next unless exists $count{$chr};
	for my $min (@q){
		if ($q >= $min){
			$count{$chr}->{$min}++;
			#print "$chr\t$min\t$q\n";
			$total{$min}++;
		}
	}
}
close IN;

print "chr\t" , join("\t" , @q) , "\n";
for my $chr (sort keys %count){
	print $chr;
	for my $min (@q){
		if (exists $count{$chr}->{$min}){
			printf("\t%.3f" , ($count{$chr}->{$min}*1000000*1000000)/($len{$chr}*$total{$min}));
		}else{
			print "\t0";
		}
	}
	print "\n";
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: chrReadsCount.pl
#
#        USAGE: ./chrReadsCount.pl  
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
#      CREATED: 03/07/18 11:17:04
#     REVISION: ---
#===============================================================================
EOF!
}



