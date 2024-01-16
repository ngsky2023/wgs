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

$/ = '>';
my %fa;
open IN , "/share/public/database/hg19/hg19.fa";
<IN>;
while (<IN>){
	chomp;
	my $seq = $_;
	s/^(\S+).*\n//;
	my $id = $1;
	s/\n//g;
	$fa{$id} = $_;
}
close IN;
$/ = "\n";

my ($bchr , $bloc , $bref , $balt , $bdep , $brate , $btype , $binfo) = ('') x 8;
my $id = '';
my $n = 1;
open IN , "$input";
while (<IN>){
	if (/^#/){
		print $_;
		next;
	}
	chomp;
	my @F = split /\t/ , $_;
	my ($chr , $loc , $ref , $alt) = @F[0,1,3,4];
	my $type;
	if (length($ref) == 1 and length($alt) == 1){
		$type = 'SNP';
	}else{
		$type = 'INDEL';
	}
	my @ti = split /:/ , $F[9];
	my ($dep , $dalt) = split /,/ , $ti[1];
	$dep += $dalt;
	
	my $brefl = length($bref);
	my $bloce = $bloc+$brefl-1 if $bloc;
	my $rate = $dalt/$dep;
	my @rate = sort {$a<=>$b} ($rate , $brate) if $brate;
	if ($bchr ne $F[0]
			or (abs($rate-$brate) > 0.05 and $rate[1]/$rate[0]>1.3)
			or ($btype eq 'SNP' and $type eq 'SNP')
			or ($loc-$bloce > 5 or $loc-$bloce < 0) ){
		if ($bchr){
			print "$bchr\t$bloc\t$id\t$bref\t$balt\t.\t$binfo\n";
		}
		($bchr , $bloc , $id , $bref , $balt , $bdep , $brate , $btype , $binfo) = (@F[0..4] , $dep , $rate , $type , join("\t" , @F[6..$#F]));
		if ($id eq '.'){
			$id = "Mut$n";
			$n++;
		}
	}else{
		my $ib = '';
		if ($loc - $bloce > 1){
			$ib = substr($fa{$bchr} , $bloce , $loc-$bloce-1);
		}
		$bref .= $ib . $ref;
		$balt .= $ib . $alt;
		$btype = 'INDEL';
		$id = "Mutm$n";
		$n++;
	}
}
close IN;
print "$bchr\t$bloc\t$id\t$bref\t$balt\t.\t$binfo\n" if $bchr;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: combinePosInVcf.pl
#
#        USAGE: ./combinePosInVcf.pl  
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
#      CREATED: 07/23/18 15:35:22
#     REVISION: ---
#===============================================================================
EOF!
}



