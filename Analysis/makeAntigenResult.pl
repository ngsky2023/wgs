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

my ($input , $help , $ann);
GetOptions(
	"i|input=s"	=>	\$input,
	"a|ann=s"	=>	\$ann,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %ti;
open(TI,"$Bin/gene_NM.list")||die "$!";
while(my $ti=<TI>){
	chomp $ti;
	my @ti=split(/\t/,$ti);
	$ti{$ti[0]}=$ti[1];
}
close(TI);


my %loc;
open IN , "$ann";
<IN>;
while (<IN>){
	my @b = split /\t/ , $_;
	next if $b[8]/($b[7]+$b[8])<0.03;
	my @transc=split(/,/,$b[15]);
	my $i=0;
	my $control=0;
	my ($ti,$exon,$cdna,$aa);
	while($i<@transc){
		my @info=split(/:/,$transc[$i]);
		if (exists $ti{$info[0]} and ($ti{$info[0]} eq $info[1]) and $ti{$info[0]}){
			$ti=$info[1];
			$exon=$info[2];
			$cdna=$info[3];
			$aa=$info[4];
			$control=1;
		}
		$i++;
	}
	if($control==0 and (not $b[15]=~/^\./)){
		my @info=split(/:/,$transc[0]);
		if($#info==4){
			$ti=$info[1];
			$exon=$info[2];
			$cdna=$info[3];
			$aa=$info[4];
		}elsif($#info==2 and $info[1]=~/^exon/){
			$ti=$info[0];
			$exon=$info[1];
			$cdna=$info[2];
			$aa="\.";
		}elsif($#info==2 and (not $info[1]=~/^exon/)){
			$ti="\.";$exon="\.";$cdna="\.";$aa="\.";
		}
	}elsif($control==0 and ($b[15]=~/^\./)){
		$ti="\.";$exon="\.";$cdna="\.";$aa="\.";
	}
	$loc{"$b[1]\t$b[3]"} = $aa;
}
close IN;

my @out;
open IN , "$input";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	unless ($F[0] =~ /chr/){
		$F[0] = "chr$F[0]";
	}
	next unless exists $loc{"$F[0]\t$F[2]"};
	next unless $F[10] eq 'Positive';
	push @out , [$F[5], $loc{"$F[0]\t$F[2]"}, $F[8], $F[7], $F[6], $F[9]];
}
close IN;

for my $ll (sort {$a->[0] cmp $b->[0]} @out){
	print join("\t" , @$ll) , "\n";
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: makeAntigenResult.pl
#
#        USAGE: ./makeAntigenResult.pl  
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
#      CREATED: 01/26/18 15:54:51
#     REVISION: ---
#===============================================================================
EOF!
}



