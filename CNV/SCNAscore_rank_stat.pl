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
my ($tcut , $ocut) = (300 , 250);
#my ($tcut , $ocut) = (0.761180826 , 0.666515166);
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my (%TSG , %OG);
my (@tsg , @og);
open IN , "$Bin/Probability_TUSON.txt";
<IN>;
while (<IN>){
	chomp;
	my ($gene , $tw , $ow) = (split /\t/ , $_)[0,6,12];
	if ($tw <= $tcut){
		$tw = $tcut - $tw + 1;
		$TSG{$gene} = $tw;
		push @tsg , $tw;
	}
	if ($ow <= $ocut){
		$ow = $ocut - $ow + 1;
		$OG{$gene} = $ow;
		push @og , $ow;
	}
}
close IN;

my %EG;
my $ew = 0;
open IN , "$Bin/Essential.txt";
while (<IN>){
	$ew++;
	chomp;
	$EG{$_} = $ew;
}
close IN;

#69.64838988			49.89833412
#91.19023569			64.0052685			107.2980291
my ($otx , $etx) = (91.19023569/64.0052685 , 91.19023569/107.2980291);
#


my ($totalNumber , $totalts , $totalos , $totaldts , $totaldos , $totalits , $totalios) = (0) x 7;
my %sc;
my $base = basename $input;
$base =~ s/\..*$//;
my %lenstat;
my $cnvn;
my $cnvlen;
my @cutlen = (0.5 , 2 , 5 , 10 , 20 , 30);
open IN , "$input";
while (<IN>){
	chomp;
	#chr16   1600000 31200000        0       -       exonic  ABAT,ABC
	my ($chr , $str , $end , $zone , $genes) = (split /\t/ , $_)[0..2,5,6];
	next if $str eq 'Start';
	next if $chr eq 'chrX' or $chr eq 'chrY';
	$cnvn++;
	$cnvlen += $end - $str;
	for my $i (0..$#cutlen){
		if ($i < $#cutlen){
			if ($end - $str >= $cutlen[$i]*1000000 and $end - $str < $cutlen[$i+1]*1000000){
				$lenstat{"$cutlen[$i]M_$cutlen[$i+1]M"}++;
				last;
			}
		}else{
			if ($end - $str >= $cutlen[$i]*1000000){
				$lenstat{"$cutlen[$i]M+"}++;
			}
		}
	}
	
	next unless $zone eq 'exonic';
	my ($type) = (/SVTYPE=(\w+);/);
	my @gene = split /[,;]/ , $genes;
	my $number = $#gene + 1;
	$totalNumber += $number;
	my ($tws , $ows , $ews) = (0 , 0 , 0);
	my ($tn , $on , $en) = (0 , 0 , 0);
	for my $gene (@gene){
		if (exists $TSG{$gene}){
			$tws += $TSG{$gene};
			$tn++;
		}
		if (exists $OG{$gene}){
			$ows += $OG{$gene};
			$on++;
		}
		if (exists $EG{$gene}){
			$ews += $EG{$gene};
			$en++;
		}
	}
	my ($ts , $os , $es) = ($tws/$number , $ows/$number , $ews/$number);

	$sc{"$chr\t$str\t$end\t$type"} = [$ts , $os , $es , $tn , $on , $en];
}
close IN;

my ($totalScore) = (0);
for my $key (sort keys %sc){
	my ($chr , $str , $end , $type) = split /\t/ , $key;
	my ($ts , $os , $es , $tn , $on , $en) = @{$sc{$key}};
	my $score = $ts - $otx*$os;# - $etx*$es;
#	print "$key\t$ts\t$os\t$es\t$score\t$tn,$on,$en\n";
	if ($type eq 'DEL'){
		$score = 0 if $score < 0;
	}elsif ($type eq 'DUP'){
		$score = 0 if $score > 0;
	}
	$totalScore += abs($score);
	#print(($ts-$os) , "\t" , ($ts - $totalts*$os/$totalos) , "\n");
}
print "$base\t$totalScore\t",$cnvlen/3095693983,"\t$cnvn";
for my $i (0..$#cutlen){
	if ($i < $#cutlen){
		if (exists $lenstat{"$cutlen[$i]M_$cutlen[$i+1]M"}){
			print "\t" , $lenstat{"$cutlen[$i]M_$cutlen[$i+1]M"};
		}else{
			print "\t0";
		}
	}else{
		if (exists $lenstat{"$cutlen[$i]M+"}){
			print "\t" , $lenstat{"$cutlen[$i]M+"};
		}else{
			print "\t0";
		}
	}
}
print "\n";


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: SCNAscore.pl
#
#        USAGE: ./SCNAscore.pl  
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
#      CREATED: 07/17/18 10:33:16
#     REVISION: ---
#===============================================================================
EOF!
}



