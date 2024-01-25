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

my ($input ,$TMB,$sig,$samplenumber,$topgene, $outdir, $help);
GetOptions(
	"i|input=s"	=>	\$input,
	"tmb|TMB=s" => \$TMB,
	"sig|sig=s" => \$sig,
	"num|num=s" => \$samplenumber,
	"o|outdir=s" => \$outdir,
        "topgene|topgene=s" => \$topgene,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %tmb;
open IN , "$TMB";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$tmb{$F[0]} = $F[1];
}
close IN;

my %p;
open IN , "$sig";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$p{$F[0]} = $F[1];
}
close IN;

my $sn = $samplenumber;
my %gs;
my %n;
open IN , "$input";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my $sample = $F[0];
	next unless exists $tmb{$sample};
	my $gene = $F[12];
	next unless exists $p{$gene};
	if ($F[14] eq 'nonsynonymous SNV' and (not exists $gs{$gene}->{$sample} or (exists $gs{$gene}->{$sample} and $gs{$gene}->{$sample} < 1))){
		$gs{$gene}->{$sample} = 1;
	}elsif ($F[11] =~ /splicing/ and (not exists $gs{$gene}->{$sample} or (exists $gs{$gene}->{$sample} and $gs{$gene}->{$sample} < 2))){
		$gs{$gene}->{$sample} = 2;
	}elsif ($F[14] eq 'stoploss' and (not exists $gs{$gene}->{$sample} or (exists $gs{$gene}->{$sample} and $gs{$gene}->{$sample} < 3))){
		$gs{$gene}->{$sample} = 3;
	}elsif ($F[14] eq 'stopgain' and (not exists $gs{$gene}->{$sample} or (exists $gs{$gene}->{$sample} and $gs{$gene}->{$sample} < 4))){
		$gs{$gene}->{$sample} = 4;
	}elsif ($F[14] =~ /nonframeshift/ and (not exists $gs{$gene}->{$sample} or (exists $gs{$gene}->{$sample} and $gs{$gene}->{$sample} < 5))){
		$gs{$gene}->{$sample} = 5;
	}elsif ($F[14] =~ /frameshift/ and (not exists $gs{$gene}->{$sample} or (exists $gs{$gene}->{$sample} and $gs{$gene}->{$sample} < 6))){
		$gs{$gene}->{$sample} = 6;
	}else{
		next;
	}
}
close IN;

for my $g (keys  %gs){
	my $v = $gs{$g};
	my @v = keys %{$v};
	$n{$g} = $#v+1;
}

my @gene = (sort {$n{$b}<=>$n{$a}} keys %n)[0..$topgene-1];
my @sample = sort {&sorts($b , $a)} keys %tmb;

open FC , ">$outdir/F3_NumMut_Sorted.txt";
print FC "SampID\tNumMut_MB\n";
for my $sample (@sample){
	print FC "$sample\t" , $tmb{$sample} , "\n";
}
close FC;

open FA , ">$outdir/F1_MutGene_Matrix.txt";
print FA "\t" , join("\t" , @sample) , "\n";
for my $gene (@gene){
	print FA $gene;
	for my $sample (@sample){
		if (exists $gs{$gene}->{$sample}){
			print FA "\t" , $gs{$gene}->{$sample};
		}else{
			print FA "\t0";
		}
	}
	print FA "\n";
}
close FA;

open FB , ">$outdir/F2_MutGene_Stats.txt";
print FB "GeneName\tMutPrev\tMutFreq\tP_Value\tMinusLogP\n";
for my $gene (@gene){
	if ($p{$gene} == 0){
		$p{$gene} = 1e-16;
	}
	print FB "$gene\t" , $n{$gene} , "\t" , $n{$gene}/$sn , "\t" , $p{$gene} , "\t" , -1*log($p{$gene})/log(10) , "\n";
}

close FB;


sub sorts{
	my ($s1 , $s2) = @_;
	for my $gene (@gene){
		unless (exists $gs{$gene}){
			return 0;
		}
		if (! exists $gs{$gene}->{$s1} and exists $gs{$gene}->{$s2}){
			return -1;
		}elsif (exists $gs{$gene}->{$s1} and ! exists $gs{$gene}->{$s2}){
			return 1;
		}elsif (! exists $gs{$gene}->{$s1} and ! exists $gs{$gene}->{$s2}){
			next;
		}elsif ($gs{$gene}->{$s1} < $gs{$gene}->{$s2}){
			return -1;
		}elsif ($gs{$gene}->{$s1} > $gs{$gene}->{$s2}){
			return 1;
		}else{
			next;
		}
	}
	return 0;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: toF1_MutGene_Matrix.pl
#
#        USAGE: ./toF1_MutGene_Matrix.pl  
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
#      CREATED: 07/26/18 11:17:55
#     REVISION: ---
#===============================================================================
EOF!
}



