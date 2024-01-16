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

my %s;
open IN , "/share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/strelka/pon/CombineVariants.txt";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my $k = join("\t" , @F[0..4]);
	$s{$k} = 1;
}
close IN;

my %indel;
open IN , "/share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/strelka/pon/indelPON.bed";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	for my $p ($F[1]..$F[2]){
		$indel{"$F[0]:$p"} = 1;
	}
}
close IN;

open IN , "/share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/strelka/result/count/m.txt";
while (<IN>){
	next if /^#/;
	chomp;
	my @F = split /\t/ , $_;
	my $k = join("\t" , @F[0..4]);
	#print STDERR "$k\n";
	$s{$k} = 1 if $F[5]>=3;
}
close IN;







open IN , "$input";
my $head = <IN>;
print $head;
my $n = 0;
while (<IN>){
	my @F = split /\t/ , $_;
	next unless $F[1] =~ /chr[\dXY]+/;
	my $k = join("\t" , @F[1..5]);
	if (exists $s{$k}){
		#print STDERR "PON\t$_";
		next;
	}

	my ($chr , $str , $end , $ref , $alt) = @F[1..5];
	if (length($ref)>=3 or length($alt)>=3){
		my $pk = 0;
		for my $p ($str..$end){
			if (exists $indel{"$chr:$p"}){
				$pk = 1;
				last;
			}
		}
		if ($pk){	
			#print STDERR "Indel\t$_";
			next;
		}
	}

	if ($F[8] >= 4 and ($F[7]+$F[8]) >= 8 and ($F[25] eq '.' or $F[25] < 0.05) and $F[26] eq '.'){
		print $_;
	}
		#	}elsif ($F[25] ne '.' and $F[25] >= 0.05){
		#		$n++;
		#	}
		#}
}
close IN;
#print STDERR "$n\n";

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: filterNonsynonymous.pl
#
#        USAGE: ./filterNonsynonymous.pl  
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
#      CREATED: 04/13/18 10:28:31
#     REVISION: ---
#===============================================================================
EOF!
}



