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

open IN , "$input";
my %sample;
while (<IN>){
	chomp;
	my ($xm , $fc , $index , $sample) = split /\t/ , $_;
	$sample{$sample}++;
	my $dir = finddir($xm , $fc);
	next unless -e "$dir/report";
	my @fq1 = </$dir/*$index.R1.clean.fastq.gz>;
	my @fq2 = </$dir/*$index.R2.clean.fastq.gz>;
	next if $#fq1 < 0 or $#fq2 < 0;
	my $fq1 = $fq1[0];
	my $fq2 = $fq2[0];
	$fq1 =~ s#/{2,}#/#g;
	$fq2 =~ s#/{2,}#/#g;
	print "$sample\t$fq1\t$fq2\tL$sample{$sample}\n";
}
close IN;


sub finddir{
	my ($xm , $fc) = @_;
	my @raw = ('/share/seq_dir/ngs/SS' , '/share/seq_dir/ngs/OC' , '/share/seq_dir/ngs' , '/share/seq_dir/OC/ngs', '/share/seq_dir/ngs/RD');
	my $datadir = '';
	for my $raw (@raw){
		my @datadir = </$raw/$xm/*$fc>;
		next if $#datadir == -1;
		$datadir = $datadir[0];
		last;
	}
	return $datadir;
}



sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: findData.pl
#
#        USAGE: ./findData.pl  
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
#      CREATED: 02/11/18 14:14:06
#     REVISION: ---
#===============================================================================
EOF!
}



