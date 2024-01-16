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
	$max = $arr[0] if $arr[0] > $max;
	$ploidy = $arr[1];
}
close IN;

my %hash;
open CNV, $file_cnv or die $!;
<CNV>;
while (<CNV>) {
	chomp;
	my @arr = split /\t/, $_;
	my ($chrom, $start, $end , $copy_number , $nMinor, $nMajor);
	if ($arr[0] =~ /"/){
		($chrom, $start, $end , $copy_number , $nMinor, $nMajor) = @arr[0..2,9..11];
	}else{
		($chrom, $start, $end , $copy_number , $nMinor, $nMajor) = @arr[1..3,9..11];
	}
	next if $nMinor eq 'NA' or $nMajor eq 'NA';
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
open OUT, "> $out" or die $!;
my $head = <IN>; chomp $head;
my @f = (2 , 'AB' , 1 , 1);
print OUT "#Chr\tPos\tRef\tAlt\tDepth_ref\tDepth_alt\tAllele_Fraction\tCNV_position\tCopyNumber\testimatedBAF\tnMajor\tnMinor\tPurity\tPloidy\t$head\n";
while (<IN>) {
	chomp;
	my @arr = split /\t/, $_;
	my ($chr, $pos, $end, $ref, $alt, $depth_ref, $depth_alt) = (@arr[1..5,7,8]);
	my $k = 0;
	if ($depth_alt) {
		my $ratio = sprintf ("%.2f", $depth_alt / ($depth_ref + $depth_alt));
		if (defined $hash{$chr}) {
			my @tmp = @{$hash{$chr}};
			foreach my $i (0..$#tmp) {
				if ($pos >= $tmp[$i][0] && $pos <= $tmp[$i][1]) {
					print OUT "$chr\t$pos\t$ref\t$alt\t$depth_ref\t$depth_alt\t$ratio\t$chr:$tmp[$i][0]-$tmp[$i][1]\t$tmp[$i][2]\t$tmp[$i][3]\t$tmp[$i][4]\t$tmp[$i][5]\t$max\t$ploidy\t$_\n";
					@f = @{$tmp[$i]}[2..5];
					$k = 1;
					last;
				}
			}
		}
	#chr10:43287148-106602537        3       AAB     2       1
	print OUT "$chr\t$pos\t$ref\t$alt\t$depth_ref\t$depth_alt\t$ratio\t$chr:NA-NA\t$f[0]\t$f[1]\t$f[2]\t$f[3]\t$max\t$ploidy\t$_\n" if $k == 0;
	}
}
close IN;
close OUT;

