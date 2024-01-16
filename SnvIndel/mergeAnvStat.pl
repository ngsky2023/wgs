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

my (@input , $help , $prx , $vcf);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"p|prx=s"	=>	\$prx,
	"v|vcf=s"	=>	\$vcf,
	"help"	=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

my %s;
open IN , "$vcf";
while (<IN>){
	next if /^#/;
	chomp;
	my @F = split /\t/ , $_;
	my ($rlen , $alen) = (length($F[3]) , length($F[4]));
	next if $F[4] =~ /,/;
	if ($rlen > $alen){
		$F[1]++;
	}
	$s{$F[0]}->{$F[1]} = 1;
}
close IN;

my $head = '';
my $ann = '';
my ($fn , $totalnum , $snp , $indel , $cox , $pop , $dups) = (0 , 0 , 0 , 0 , 0 , 0 , 0);
my ($dbsnp , $total) = (0 , 0);
for my $file (@input){
	open IN , "$file";
	while (<IN>){
		if (/^SampleName/){
			$head .= $_ if $fn == 0;
		}else{
			chomp;
			my @F = split /\t/ , $_;
			next unless exists $s{$F[1]}->{$F[2]};
			next if $F[14] =~ /unknown/;
			$total++;

			if ($F[43] ne '.' and $F[43] >= 0.05){
				$pop++;
				next;
			}
			if ($F[84] ne '.'){
				$dups++;
				next;
			}
			if ($F[23] =~ /rs/){
				$dbsnp++;
			}
			$ann .= "$_\n";
			my $rl = length($F[4]);
			my $al = length($F[5]);
			if ($F[4] =~ /,/){
				$cox++;
			}elsif ($F[4] eq '-' or $F[5] eq '-'){
				$indel++;
			}elsif ($rl == $al and $al == 1){
				$snp++;
			}else{
				$cox++;
			}
			$totalnum++;
		}
	}
	close IN;
	$fn++;
}

open ANN , ">$prx.xls";
print ANN "$head$ann";
close ANN;

my $dbb = sprintf("%.2f" , $dbsnp*100/$totalnum);
open STAT , ">$prx.stat";
print STAT << "EOF!";
Number of mutation\t$totalnum
Number of SNP\t$snp
Number of INDEL\t$indel
Number of complex\t$cox
Number in dbSNP\t$dbsnp ($dbb%)

Number in PopFreqMax\t$pop
Number in genomicSuperDups\t$dups
EOF!
close STAT;


sub help{
print << "EOF!";
#===============================================================================
#===============================================================================
EOF!
}



