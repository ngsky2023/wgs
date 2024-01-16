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

my ($input , $help , $t);
GetOptions(
	"i|input=s"	=>	\$input,
	"t|t=s"	=>	\$t,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}
my $win = 100000;

my %cfs;
my %cfsize;
open IN , "$input";
while (<IN>){
	next if /Sample/;
	chomp;
	my @F = split /\t/ , $_;
	my @bin = bin(@F[2,3]);
	my ($sample) = ($F[0] =~ /(\d+)/);
	for my $i (@bin){
		$cfs{$sample}->{$F[1]}->{$i} = 1;
		$cfsize{$sample} += $win;
	}
}
close IN;

my %tus;
my %tusize;
my %same;
open IN , "$t";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my @bin = bin(@F[2,3]);
	my ($sample) = ($F[0] =~ /(\d+)/);
	for my $i (@bin){
		$tus{$sample}->{$F[1]}->{$i} = 1;
		$same{$sample} += $win if exists $cfs{$sample}->{$F[1]}->{$i};
		$tusize{$sample} += $win;
	}

}
close IN;

for my $sample (sort keys %same){
	my ($cf , $tu , $same) = (0 , 0 , 0);
	if (exists $cfsize{$sample}){
		$cf = $cfsize{$sample};
	}
	if (exists $tusize{$sample}){
		$tu = $tusize{$sample};
	}
	if (exists $same{$sample}){
		$same = $same{$sample};
	}
	print "$sample\t$tu\t$cf\t$same\t" , $same/$tu , "\t" , $same/$cf , "\n";
}


sub bin{
	my ($str , $end) = @_;
	my @bin = ((int($str/$win))..(int($end/$win)));
	return @bin;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: compareSize.pl
#
#        USAGE: ./compareSize.pl  
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
#      CREATED: 06/22/18 11:20:38
#     REVISION: ---
#===============================================================================
EOF!
}



