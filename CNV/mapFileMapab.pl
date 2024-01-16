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

my ($chr) = ('');
my %mapab;
open IN , "/share/work1/wangrr/DB/hg19/binCor/wgEncodeCrgMapabilityAlign36mer.20K.wig";
while (<IN>){
	chomp;
	if (/^fixedStep chrom=(chr\w+)/){
		$chr = $1;
	}else{
		push @{$mapab{$chr}} , $_;
	}
}
close IN;

mkdir $dir unless -d $dir;
my @file  = <$input/*.20K>;

for my $file (@file){
	my $base = basename($file);
	my $total = 0;
	my $out = '';
	open IN , "$file";
	while (<IN>){
		next if /^#/;
		chomp;
		my ($chr , @F) = split /\t/ , $_;
		$out .= $chr;
		for my $i (0..$#F){
			my $d = 0;
			if ($i <= $#{$mapab{$chr}} and $mapab{$chr}->[$i] != 0){
				$d =  int($F[$i]/$mapab{$chr}->[$i]);
			}else{
				$d = $F[$i];	
			}
			$out .= "\t$d";
			$total += $d;
		}
		$out .= "\n";
	}
	close IN;
	$out = "#uniqMapped\t$total\n" . $out;
	open O , ">$dir/$base";
	print O $out;
	close O;

	$out = '';
	open IN , "$file.GC";
	while (<IN>){
		next if /^#/;
		chomp;
		my ($chr , @F) = split /\t/ , $_;
		$out .= $chr;
		for my $i (0..$#F){
			my $d = 0;
			if ($i <= $#{$mapab{$chr}} and $mapab{$chr}->[$i] != 0){
				$d =  int($F[$i]/$mapab{$chr}->[$i]);
			}else{
				$d = $F[$i];	
			}
			$out .= "\t$d";
		}
		$out .= "\n";
	}
	close IN;
	$out = "#uniqMapped\t$total\n" . $out;
	open O , ">$dir/$base.GC";
	print O $out;
	close O;

}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mapFileMapab.pl
#
#        USAGE: ./mapFileMapab.pl  
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
#      CREATED: 08/21/2018 05:35:26 PM
#     REVISION: ---
#===============================================================================
EOF!
}



