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
open IN , "/share/work1/wangrr/DB/hg19/hg19AddVirus.fa";
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

my ($bchr , $bref , $balt , $btype , $binfo , $bpid) = ('') x 5;
my ($bloc , $baltDep , $brate) = (0) x 3;
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
		if (length($ref) > length($alt) and length($alt) == 1){
			$type = 'del';
		}elsif (length($ref) < length($alt) and length($ref) == 1 and substr($ref , 0 , 1) eq substr($alt , 0 , 1)){
			$type = 'ins';
		}
	}
	my @ti = split /:/ , $F[9];
	my $pid = '';
	if ($F[8] =~ /:PID:/){
		$pid = $ti[12];
	}else{
		$pid = '';
	}
	my ($dep , $dalt) = split /,/ , $ti[1];
	my $minalt = (sort {$a<=>$b} ($dalt , $baltDep))[0];
	$dep += $dalt;
	
	my $brefl = length($bref);
	my $bloce = $bloc+$brefl-1 if $bloc;
	my $rate = $dalt/$dep;
	my @rate = sort {$a<=>$b} ($rate , $brate) if $brate;
	if ($pid eq '' or $pid ne $bpid){
		if ($bchr){
			print "$bchr\t$bloc\t$id\t$bref\t$balt\t.\t$binfo\n";
		}
		($bchr , $bloc , $id , $bref , $balt , $baltDep , $brate , $btype , $bpid , $binfo) = (@F[0..4] , $dalt , $rate , $type , $pid , join("\t" , @F[6..$#F]));
		if ($id eq '.'){
			$id = $chr."Mut$n";
			$n++;
		}
	}else{
		my $ib = '';
		if ($loc - $bloce > 1){
			$ib = substr($fa{$bchr} , $bloce , $loc-$bloce-1);
		}elsif ($loc - $bloce == 0 and $type eq 'ins'){
			$ref =~ s/^\w//;
			$alt =~ s/^\w//;
		}
		$bref .= $ib . $ref;
		$balt .= $ib . $alt;
		$btype = 'INDEL';
		$id = $chr."Mutm$n";
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



