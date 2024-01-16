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

my ($input , $help , $t , $n);
GetOptions(
	"i|input=s"	=>	\$input,
	"t|t=s"	=>	\$t,
	"n|n=s"	=>	\$n,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

open IN , "sed '1d' $input | paste - - |";
my %data;
my %t;
while (<IN>){
	#chr1    69085   70013   OR4F5   528     96.6078 Ct180527-T      chr1    69085   70013   OR4F5   432     80.7015 Ct180527
	chomp;
	my @F = split /\t/ , $_;
	next if $F[5] < 10 or $F[12] < 10 or $F[5] > 1e+7 or $F[12] > 1e+7;
	$t{$t} += $F[5]*($F[2]-$F[1]);
	$t{$n} += $F[12]*($F[2]-$F[1]);
	$F[1]++;
	push @{$data{$t}} , ["$F[0]:$F[1]-$F[2]" , $F[5]];
	push @{$data{$n}} , ["$F[0]:$F[1]-$F[2]" , $F[12]];
}
close IN;

my @td = @{$data{$t}};
my @nd = @{$data{$n}};
my $dt = $t{$t};
my $dn = $t{$n};
print "Target\t$t\n";
for my $i (0..$#td){
	my $td = $td[$i];
	my $nd = $nd[$i];
	if ($nd->[1] == 0 or $td->[1] == 0){
		print $td->[0] , "\tNaN\n";
	}else{
		print $td->[0] , "\t" , log(($td->[1]/$dt)/($nd->[1]/$dn))/log(2) , "\t" , $td->[1] , "\t" , $nd->[1] , "\t$dt\t$dn\n";
	}
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: LRR_Calculation.pl
#
#        USAGE: ./LRR_Calculation.pl  
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
#      CREATED: 11/06/2018 03:14:55 PM
#     REVISION: ---
#===============================================================================
EOF!
}



