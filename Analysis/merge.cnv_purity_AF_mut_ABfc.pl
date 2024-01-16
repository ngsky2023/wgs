#!/usr/bin/env perl
use strict;
use warnings;

my ($file_snv, $file_cnv, $file_purity, $out) = @ARGV;
if (@ARGV < 4) {
	die "\n\tperl $0 /home/liuya/workdir/project/exon/somatic_mutation/IBFC2017199/analysis/08CCF/test/H2_15S06464_5.merge.snp_indel.xls /home/liuya/workdir/project/exon/somatic_mutation/IBFC2017199/analysis/03cnv_control_FREEC/01callcnv/H2_15S06464_5/H2_15S06464_5.bam.pileup_ratio.txt /home/liuya/workdir/project/exon/somatic_mutation/IBFC2017199/analysis/09purity_sequenza/H2_15S06464_5/H2_15S06464_5_confints_CP.txt H2_15S06464_5.AF_CN_purity.xls\n\n";
	}

my ($max, $ploidy) = (0, 0);
open P, $file_purity or die $!;
<P>;
while (<P>) {
	chomp;
	my @arr = split /\t/, $_;
	$max = $arr[3];
	$ploidy = $arr[4];
}
close IN;

my %hash;
open CNV, $file_cnv or die $!;
<CNV>;
while (<CNV>) {
	chomp;
	my @arr = split /\t/, $_;
	my ($chrom, $start, $end , $copy_number , $nMinor, $nMajor) = @arr[0,1,1,4,7,7];
	$end += 50000-1;
	my $ab = $nMinor;
	$nMinor = 0;
	$nMajor = 0;
	next if $nMinor eq 'NA' or $nMajor eq 'NA';

	for my $base (split // , $ab){
		if ($base eq 'A'){
			$nMajor++;
		}elsif ($base eq 'B'){
			$nMinor++;
		}
	}
	($nMinor, $nMajor) = sort {$a<=>$b} ($nMinor, $nMajor);
	$chrom =~ s/"//g;
	my $genotype = '-';
	$genotype = ('A' x $nMajor) . ('B' x $nMinor);
	my $chr = $chrom;
	$chr = "chr$chrom" if $chrom !~ /chr/;
	push @{$hash{$chr}}, [$start, $end, $copy_number, $genotype, $nMajor, $nMinor];
}
close CNV;

open IN, $file_snv or die $!;
my $sample = $file_snv;
$sample =~ s/\..*$//;
<IN>;
open OUT, "> $out" or die $!;
print OUT "#Chr\tPos\tRef\tAlt\tDepth_ref\tDepth_alt\tAllele_Fraction\tCNV_position\tCopyNumber\testimatedBAF\tnMajor\tnMinor\tPurity\tPloidy\n";
while (<IN>) {
	chomp;
	my @arr = split /\t/, $_;
	#my ($sample , $chr, $pos, $ref, $alt, $depth_ref, $depth_alt) = (@arr[0..2] , 'A' , 'A' , @arr[4,5]);
	##ZL126676	chr5	1295228	1295228	G	A	375415	232073	0.618177217212951
	my ($chr , $pos, $ref, $alt, $depth_ref, $depth_alt) = (@arr[1,2,4,5] , $arr[6]-$arr[7] , $arr[7]);
	$depth_ref = 0 if $depth_ref < 0;
	if ($depth_alt) {
		my $ratio = sprintf ("%.2f", $depth_alt / ($depth_ref + $depth_alt));
		$ratio = 0.99 if $ratio > 0.99;
		if (defined $hash{$chr}) {
			my @tmp = @{$hash{$chr}};
			foreach my $i (0..$#tmp) {
				if ($pos >= $tmp[$i][0] && $pos <= $tmp[$i][1]) {
					print OUT "$chr\t$pos\t$ref\t$alt\t$depth_ref\t$depth_alt\t$ratio\t$chr:$tmp[$i][0]-$tmp[$i][1]\t$tmp[$i][2]\t$tmp[$i][3]\t$tmp[$i][4]\t$tmp[$i][5]\t$max\t$ploidy\t$_\n";
				}
			}
		}
	}
}
close IN;
close OUT;

