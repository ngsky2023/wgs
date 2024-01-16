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

my ($input , $help , $bed , $syn);
GetOptions(
	"i|input=s"	=>	\$input,
	"b|bed=s"	=>	\$bed,
	"s|syn=s"	=>	\$syn,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my $l = 0;
open IN , "$bed";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$l+= $F[2]-$F[1]+1;
}
close IN;

my %s;
my %h;
open IN , "$input";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$s{$F[0]}++;
	$h{$F[0]}++ if $F[18] eq '1';
}
close IN;

my %y;
open IN , "$syn";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$y{$F[0]}++;
}
close IN;



print "Sample\tTMB/Mb\n";
for my $s (sort keys %s){
	$y{$s} = 0 unless exists $y{$s};
	$h{$s} = 0 unless exists $h{$s};
	printf("$s\t%.2f\t%.2f\t%.2f\t%.2f\n" , $s{$s}*1000000/$l , ($s{$s}+$y{$s})*1000000/$l , ($s{$s}-$h{$s})*1000000/$l , ($s{$s}+$y{$s}-$h{$s})*1000000/$l);
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: tmb.pl
#
#        USAGE: ./tmb.pl  
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
#      CREATED: 11/05/2018 10:07:30 AM
#     REVISION: ---
#===============================================================================
EOF!
}



