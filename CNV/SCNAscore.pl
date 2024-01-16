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
#my ($tcut , $ocut) = (300 , 250);
#my ($tcut , $ocut) = (0.761180826 , 0.666515166);
my ($tcut , $ocut) = (0.761180826 , 0.599787207);
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
	my ($gene , $tw , $ow) = (split /\t/ , $_)[0,5,11];
	if ($tw >= $tcut){
		$TSG{$gene} = $tw;
		push @tsg , $tw;
	}
	if ($ow >= $ocut){
		$OG{$gene} = $ow;
		push @og , $ow;
	}
}
close IN;

@tsg = sort {$a<=>$b} @tsg;
@og = sort {$a<=>$b} @og;
for my $gene (sort keys %TSG){
	my $tw = $TSG{$gene};
	$TSG{$gene} = 300*$tw/$tsg[-1];
}
for my $gene (sort keys %OG){
	my $ow = $OG{$gene};
	$OG{$gene} = 300*$ow/$og[-1];
}

my ($totalNumber , $totalts , $totalos , $totaldts , $totaldos , $totalits , $totalios) = (0) x 7;
my %sc;
open IN , "$input";
while (<IN>){
	chomp;
	#chr16   1600000 31200000        0       -       exonic  ABAT,ABC
	my ($chr , $str , $end , $zone , $genes) = (split /\t/ , $_)[0..2,5,6];
	next if $chr eq 'chrX' or $chr eq 'chrY';
	next unless $zone eq 'exonic' and $end - $str >= 5000000;
	my ($type) = (/SVTYPE=(\w+);/);
	my @gene = split /,/ , $genes;
	my $number = $#gene + 1;
	$totalNumber += $number;
	my ($tws , $ows) = (0 , 0);
	my ($tn , $on) = (0 , 0);
	for my $gene (@gene){
		if (exists $TSG{$gene}){
			$tws += $TSG{$gene};
			$tn++;
		}
		if (exists $OG{$gene}){
			$ows += $OG{$gene};
			$on++;
		}
	}
	my ($ts , $os) = ($tws/$number , $ows/$number);
	$totalts += $ts;
	$totalos += $os;
	if ($type eq 'DEL' ){#and $ts > $os){
		$totaldts += $ts;
		$totaldos += $os;
	}elsif ($type eq 'DUP' ){#and $ts < $os){
		$totalits += $ts;
		$totalios += $os;
	}
	$sc{"$chr\t$str\t$end\t$type"} = [$ts , $os , $tn , $on];
}
close IN;

for my $key (sort keys %sc){
	my ($ts , $os , $tn , $on) = @{$sc{$key}};
	print "$key\t$ts\t$os\t" , ($ts-$os) , "\t$tn,$on\n";
	#print(($ts-$os) , "\t" , ($ts - $totalts*$os/$totalos) , "\n");
}
print $totaldts - $totaldos , "\n";
print $totalits - $totalios , "\n";

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



