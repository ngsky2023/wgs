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

if ($help or $#ARGV<0){
	&help;
	exit;
}

my $h = <<"EOF!";
##fileformat=VCFv4.2
###FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths (number of reads in each observed allele)">
###FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
###FORMAT=<ID=FT,Number=1,Type=String,Description="Variant filters">
###INFO=<ID=t_ref_count,Number=1,Type=Integer,Description="Tumour ref count">
###INFO=<ID=t_alt_count,Number=1,Type=Integer,Description="Tumour alt count">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TUMOR
EOF!

for my $file (@ARGV){
	my $base = basename $file;
	open IN , "$file";
	open O , ">$base";
	my $sample = $base;
	$sample =~ s/\.vcf//;
	print O $h;
	while (<IN>){
		next if /^#/;
		chomp;
		my @F = split /\t/ , $_;
		$F[0] =~ s/^chr//;
		my @j = split /:/ , $F[9];
		my @n = split /,/ , $j[1];
		$F[7] = "t_alt_count=$n[1];t_ref_count=$n[0]";
		print O join("\t" , @F) , "\n";
	}
	if (exists $d{$sample}){
		print O $d{$sample};
	}
	close IN;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: tranvcf.pl
#
#        USAGE: ./tranvcf.pl  
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
#      CREATED: 06/06/2020 12:47:55 PM
#     REVISION: ---
#===============================================================================
EOF!
}



