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

my ($input , $help , $dir);
GetOptions(
	"i|input=s"	=>	\$input,
	"d|dir=s"	=>	\$dir,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}


my $head = << "EOF!";
##fileformat=VCFv4.2
##reference=file:///share/work1/wangrr/DB/hg19/hg19AddVirus.fa
##contig=<ID=chrM,length=16571>
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr4,length=191154276>
##contig=<ID=chr5,length=180915260>
##contig=<ID=chr6,length=171115067>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr8,length=146364022>
##contig=<ID=chr9,length=141213431>
##contig=<ID=chr10,length=135534747>
##contig=<ID=chr11,length=135006516>
##contig=<ID=chr12,length=133851895>
##contig=<ID=chr13,length=115169878>
##contig=<ID=chr14,length=107349540>
##contig=<ID=chr15,length=102531392>
##contig=<ID=chr16,length=90354753>
##contig=<ID=chr17,length=81195210>
##contig=<ID=chr18,length=78077248>
##contig=<ID=chr19,length=59128983>
##contig=<ID=chr20,length=63025520>
##contig=<ID=chr21,length=48129895>
##contig=<ID=chr22,length=51304566>
##contig=<ID=chrX,length=155270560>
##contig=<ID=chrY,length=59373566>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Tandem Duplication">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
EOF!

my $sample = '';
mkdir $dir unless -d $dir;
open IN , "$input";
while (<IN>){
	#ZL201029T       1       724621  726809  2189    0.62    -0.48   1       DLOH    1       0       1       2       0.50
	#1CNVCF328754    chr1    238220000       239580000       1.061   del     1.36    0.12245928332543
	chomp;
	my @F = split /\t/ , $_;
	next if $F[4] == 2;
	if ($sample ne $F[0]){
		close O if $sample;
		$sample = $F[0];
		open O , ">$dir/$sample.cnv.vcf";
		print O $head;
	}
	my $type = '';
	if ($F[4] > 2){
		$type = 'DUP';
	}else{
		$type = 'DEL';
	}
	my $len = $F[3]-$F[2];
	print O "$F[1]\t$F[2]\tCNV\t.\t<$type>\t.\tPASS\tEND=$F[3];SVTYPE=$type;SVLEN=$len;CN=$F[4];\tGT\t0/1\n";
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: TitanSegs2Vcf.pl
#
#        USAGE: ./TitanSegs2Vcf.pl  
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
#      CREATED: 03/13/18 14:48:15
#     REVISION: ---
#===============================================================================
EOF!
}




