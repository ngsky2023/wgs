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

my ($input , $vcf , $help);
GetOptions(
	"i|input=s"	=>	\$input,
	"v|vcf=s"	=>	\$vcf,
	"help"	=>	\$help,
);

if ($help or ! $vcf){
	&help;
	exit;
}

my %sample;
open IN , "$input";
<IN>;
while (<IN>){
	my ($sample , $chr , $str , $end , $ref , $var, $tref , $talt, $zona , $class) = (split /\t/ , $_)[0..5,7,8,11,14];
	next if $talt/($tref + $talt) < 0.05;
	next if ((not $zona=~/exonic/) and (not $zona=~/splicing/));
	next if  $zona eq 'ncRNA_exonic';
	next if (($class eq "synonymous SNV") and (not $zona=~/splicing/));

	$str-- if $var eq '-';
	$sample{"$chr\t$str"} = 1;
}
close IN;

open IN , "$vcf";
while (<IN>){
	if (/^#/){
		print $_;
		next;
	}
	my @F = split /\t/ , $_;
	$F[2] = '.';
	if (exists $sample{"$F[0]\t$F[1]"}){
		$F[0] =~ s/chr//;
		$F[6] = 'PASS';
		print join("\t" , @F);
	}
}
close IN;



sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: pickVcfByReport.pl
#
#        USAGE: ./pickVcfByReport.pl  
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
#      CREATED: 09/25/17 14:20:09
#     REVISION: ---
#===============================================================================
EOF!
}



