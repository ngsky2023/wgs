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

my ($input , $help , $prx);
GetOptions(
	"i|input=s"	=>	\$input,
	"p|prx=s"	=>	\$prx,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %geneloc;
open IN , "/share/work1/wangrr/DB/hg19/hg19.gff";
while (<IN>){
	chomp;
	next unless /\tgene\t/;
	my @F = split /\t/ , $_;
	my ($gene) = (/ID=([^;]+);/);
	$geneloc{$F[0]}->{$gene} = [@F[3,4]];
}
close IN;


my %mate;
my %mateo;
open FU , ">$prx.fusion.txt";
open RE , ">$prx.region.txt";
open IN , "$input";
my $head = <IN>;
print RE $head;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	#next unless $F[30] eq 'PASS';
	#next if $F[5] eq 'intergenic';
	if (/SVTYPE=DEL;/ or /SVTYPE=DUP;/){
		my ($len) = (/SVLEN=([\-\d]+);/);
		$len = abs($len);
		next if $len > 500000;
		print RE "$_\n";
	}elsif (/SVTYPE=INV;/){
		my ($gene1 , $gene2);
		for my $gene (split /,/ , $F[6]){
			next unless exists $geneloc{$F[0]}->{$gene};
			my ($gstr , $gend) = @{$geneloc{$F[0]}->{$gene}};
			if ($gstr < $F[1] and $F[2] < $gend){
				next;
			}elsif ($gstr < $F[1] and $F[1] < $gend){
				$gene1 = $gene;
			}elsif ($gstr < $F[2] and $F[2] < $gend){
				$gene2 = $gene;
			}
		}
		if ($gene1 and $gene2){
			print FU "$F[0]\t$F[1]\t$gene1\t$F[0]\t$F[2]\t$gene2\n";
		}else{
			print RE "$_\n";
		}
	}elsif (/SVTYPE=BND;/){
		if (exists $mate{$F[25]}){
			if ($F[5] eq 'intronic' or $F[5] eq 'exonic'){
				print FU join("\t" , (@{$mate{$F[25]}}[1..3] , @F[0,1,6])) , "\n";
			}else{
				print RE $mate{$F[25]}->[0] , "\n";
				print RE "$_\n";
			}
		}elsif (exists $mateo{$F[25]}){
			print RE $mateo{$F[25]} , "\n";
			print RE "$_\n";
		}else{
			my ($md) = (/MATEID=([^;]+);/);
			if ($F[5] eq 'intronic' or $F[5] eq 'exonic'){
				$mate{$md} = [$_ , @F[0,1,6]];
			}else{
				$mateo{$md} = $_;
			}
		}
	}
}
close IN;
close FU;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: inrMantaVcfAnnot.pl
#
#        USAGE: ./inrMantaVcfAnnot.pl  
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
#      CREATED: 01/24/18 10:53:14
#     REVISION: ---
#===============================================================================
EOF!
}



