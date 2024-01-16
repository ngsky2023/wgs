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

my (@input , $help , $prx);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"p|prx=s"	=>	\$prx,
	"help"	=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

my $head = '';
my $vcf = '';
my ($fn , $totalnum , $dbnum , $snp , $indel , $cox) = (0 , 0 , 0 , 0 , 0 , 0);
for my $file (@input){
	open IN , "$file";
	while (<IN>){
		if (/^#/){
			$head .= $_ if $fn == 0;
		}else{
			$vcf .= $_;
			chomp;
			my @F = split /\t/ , $_;
			my $rl = length($F[3]);
			my $al = length($F[4]);
			if ($F[4] =~ /,/){
				$cox++;
			}elsif ($rl == $al and $al == 1){
				$snp++;
			}elsif (($rl == 1 and $al > 1) or ($rl > 1 and $al == 1)){
				$indel++;
			}else{
				$cox++;
			}
			$dbnum++ if $F[2] =~ /rs/;
			$totalnum++;
		}
	}
	close IN;
	$fn++;
}

open VCF , ">$prx.vcf";
print VCF "$head$vcf";
close VCF;

my $dbb = sprintf("%.2f" , $dbnum*100/$totalnum);
open STAT , ">$prx.stat";
print STAT << "EOF!";
Number of mutation\t$totalnum
Number of SNP\t$snp
Number of INDEL\t$indel
Number of complex\t$cox
EOF!
close STAT;



sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeVcfStat.pl
#
#        USAGE: ./mergeVcfStat.pl  
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
#      CREATED: 03/01/18 14:11:34
#     REVISION: ---
#===============================================================================
EOF!
}



