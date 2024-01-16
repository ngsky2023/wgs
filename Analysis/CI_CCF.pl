#!/usr/bin/env perl
use strict;
use warnings;
use lib '/share/work1/wangrr/local/Plib/lib/perl5/';
use Statistics::R;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Math::BigFloat;

my ($file, $out) = @ARGV;
if (@ARGV < 2) {
	die "\n\tperl $0 /home/liuya/workdir/project/exon/somatic_mutation/IBFC2017199/analysis/08CCF/test/H2_15S06464_5.AF_CN_purity.xls H2_15S06464_5.CI_CCF.xls\n\n";
	}

open IN, $file or die $!;
open OUT, "> $out" or die $!;
my $head = <IN>;
chomp $head;
my @title = split /\t/, $head;
my @first = splice (@title, 0, 14);
my $first = join "\t", @first;
my $last = join "\t", @title;
print OUT "$first\tAF_Lo\tAF_Hi\tCCF\tCCF_Lo\tCCF_Hi\tClonalStatus\t$last\n";
while (<IN>) {
	chomp;
	my @arr = split /\t/, $_;
	my ($depth_ref, $depth_alt) = ($arr[4], $arr[5]);
	my ($af, $cn, $purity, $nMajor) = ($arr[6], $arr[8], $arr[12], $arr[10]);
	next if $cn == 0;
	my ($af_lo, $af_hi, %nchr);
	my ($ccf, $ccf_lo, $ccf_hi, $stat) = ('NA', 'NA', 'NA', 'NA');

## R to calculate AF_CI
	my $R = Statistics::R->new();
	$R -> send(qq`ref<-rep(0,$depth_ref)`);
	$R -> send(qq`var<-rep(1,$depth_alt)`);
	$R -> send(qq`cov<-$depth_ref+$depth_alt`);
	$R -> send(qq`total<-c(ref,var)`);
	$R -> send(qq`s<-0`);
	$R -> send(qq`for(n in 1:10000){
		boot<-sample(total,size=cov,replace=TRUE)
		s[n]<-sum(boot)/cov
		}`);
	$R -> send(qq`ci<-quantile(s, c(.025,.975))`);
#	$R -> send(qq`x<-as.numeric(min(ci),max(ci))`);
#	my ($ci_min, $ci_max) = $R->get('x');
	$af_lo = $R->get('min(ci)');
	$af_hi = $R->get('max(ci)');
	$R->stop();

## calculate CCF
	my $n_mut = $af*(1/$purity)*($purity*$cn+(1-$purity)*2);
	my $n_mut_lo = $af_lo*(1/$purity)*($purity*$cn+(1-$purity)*2);
	my $n_mut_hi = $af_hi*(1/$purity)*($purity*$cn+(1-$purity)*2);
	
	if ($cn<=0 or $nMajor eq "NA") {
		($ccf, $ccf_lo, $ccf_hi, $stat) = ('NA', 'NA', 'NA', 'NA');
	}

## Maximum Likelihood
	if ($nMajor ne "NA") {
	for (my $i=0; $i<=$nMajor; $i++) {
		my $f = $i*$purity/($purity*$cn+(1-$purity)*2);
		my $n = $depth_ref+$depth_alt;
		my $k = $depth_alt;
		## input $n, output $x $x=$n
		my $x = Math::BigFloat->new($n);
		## $x->bnok($y);   # x over y (binomial coefficient n over k)
		my $z = $x->bnok($k);
		## Maximum Likelihood
		my $l = $z*($f**$depth_alt)*((1-$f)**$depth_ref);
		$nchr{$l}=$i;
	}
	my @l = keys %nchr;
	## max in perl module: List::Util qw(first max maxstr min minstr reduce shuffle sum);
	my $ml = max @l;
	my $nchr=$nchr{$ml};
	if ($nchr == 0){
		print STDERR  "$_\t$n_mut\t$n_mut_lo\t$n_mut_hi\n";
		next;

	}
	$ccf = $n_mut/$nchr;
	$ccf_lo = $n_mut_lo/$nchr;
	$ccf_hi = $n_mut_hi/$nchr;
	
	## clonal status, $ccf_hi >= 1, 95% CI in area
	if ($ccf_hi >= 1) {
		$stat = "Clonal";
	} else {
		$stat = "Subclonal";
	}
	}

## print result
	my @tmp = splice (@arr, 0, 14);
	my $first_tmp = join "\t", @tmp;
	my $last_tmp = join "\t", @arr;
	print OUT "$first_tmp\t$af_lo\t$af_hi\t$ccf\t$ccf_lo\t$ccf_hi\t$stat\t$last_tmp\n";
}
close IN;
close OUT;

