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

my %to;
open IN , "$Bin/gene_NM.list";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$to{$F[0]} = $F[1];
}
close IN;

open IN , "$input";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	#11,12,13,14,15;
	#6,7,8,9,10
	unless ($F[11] =~ /exonic|splicing/ and $F[11] !~ /RNA/){
		print "$_\n";
		next;
	}
	my $trans = "";
	my $exon = "";
	my $cds = "";
	my $aa = ".";
	$F[15] = $F[13] if $F[15] eq ".";
	$trans = $to{$F[12]} if exists $to{$F[12]};
	for my $ll (split /,/, $F[15]){
		my @l = split /:/ , $ll;
		shift @l if $F[11] =~ /exonic/;
		$cds = $l[2];
		$exon = $l[1];
		if ($#l >= 3){
			$aa = $l[3];
			if ($trans eq ""){
				$trans = $l[0];
				last;
			}elsif ($l[0] eq $trans){
				last;
			}
		}
	}
	$F[13] = "$trans:$exon:$cds:$aa";
	print join("\t" , @F) , "\n";
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: MainTras.pl
#
#        USAGE: ./MainTras.pl  
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
#      CREATED: 06/01/2020 03:02:35 PM
#     REVISION: ---
#===============================================================================
EOF!
}



