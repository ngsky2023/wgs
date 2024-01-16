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

my ($input , $help , $bed);
#$bed = "/share/work2/yangrutao/workdir/research/databases/55panel_primer_IAD117925_197_384WellPlateDataSheet.bed";
GetOptions(
	"i|input=s"	=>	\$input,
	"b|bed=s"	=>	\$bed,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %mask;
open IN , "$bed";
while (<IN>){
	chomp;
	my ($chr , $str , $end) = split /\t/ , $_;
	for my $p ($str..$end){
		$mask{"$chr:$p"} = 1;
	}
}
close IN;

open IN , "samtools view -h $input |";
while (<IN>){
	if (/^\@/){
		print $_;
		next;
	}
	chomp;
	my @F = split /\t/ , $_;
	my $flag = sprintf("%b" , $F[1]);
	my $chain;
	if ($flag =~ /1\d\d\d\d$/){
		$chain = -1;
	}else{
		$chain = 1;
	}
	my ($chr , $pos) = @F[2,3];
	if ($chain == 1){
		my ($mstr , $mend) = (0 , 0);
		if (exists $mask{"$chr:$pos"}){
			while ($F[5] =~ /(\d+)([MIDSH])/g){
				my $k = 0;
				if ($2 eq 'M'){
					for my $p ($pos..($pos+$1-1)){
						if (exists $mask{"$chr:$p"}){
							$mend++;
						}else{
							$k = 1;
							last;
						}
					}
					$pos += $1;
				}elsif ($2 eq 'I'){
					if (exists $mask{"$chr:$pos"}){
						$mend += $1;
					}else{
						last;
					}
				}elsif ($2 eq 'D'){
					for my $p ($pos..($pos+$1-1)){
						if (exists $mask{"$chr:$p"}){
						}else{
							$k = 1;
							last;
						}
					}
					$pos += $1;
				}elsif ($2 eq 'S'){
					$mend += $1;
				}
				last if $k;
			}
			substr($F[10], $mstr, $mend-$mstr) = '!' x ($mend-$mstr);
#			print STDERR join("\t" , @F) , "\n" if length($F[9]) != length($F[10]);
			print  join("\t" , @F) , "\n";
		}else{
			print "$_\n";
			next;
		}
	}else{
		my ($maplen , $m) = mode($F[5]);
		my $mstr = length($F[9]);
		my $mend = $mstr;
		$pos = $pos + $maplen - 1;
		if (exists $mask{"$chr:$pos"}){
			while ($m =~ /(\d+)([MIDSH])/g){
				my $k = 0;
				if ($2 eq 'M'){
					for my $p (reverse(($pos-$1+1)..$pos)){
						if (exists $mask{"$chr:$p"}){
							$mstr--;
						}else{
							$k = 1;
							last;
						}
					}
					$pos -= $1;
				}elsif ($2 eq 'I'){
					if (exists $mask{"$chr:$pos"}){
						$mstr -= $1;
					}else{
						last;
					}
				}elsif ($2 eq 'D'){
					for my $p (reverse(($pos-$1+1)..$pos)){
						if (exists $mask{"$chr:$p"}){
						}else{
							$k = 1;
							last;
						}
					}
					$pos -= $1;
				}elsif ($2 eq 'S'){
					$mstr -= $1;
				}
				last if $k;
			}
			substr($F[9], $mstr, $mend-$mstr) = 'N' x ($mend-$mstr);
			substr($F[10], $mstr, $mend-$mstr) = '!' x ($mend-$mstr);
#			print STDERR join("\t" , @F) , "\n" if length($F[9]) != length($F[10]);
			print  join("\t" , @F) , "\n";
		}else{
			print "$_\n";
			next;
		}
	}
}
close IN;


sub mode{
	my ($m) = @_;
	my ($len) = (0);
	my $out = '';
	while ($m =~ /(\d+)([MIDSH])/g){
		if ($2 eq 'M'){
			$len += $1;
		}elsif ($2 eq 'D'){
			$len += $1;
		}
		$out = "$1$2$out";
	}
	return ($len , $out);
}



sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: maskPrimerBase.pl
#
#        USAGE: ./maskPrimerBase.pl  
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
#      CREATED: 01/30/18 13:51:49
#     REVISION: ---
#===============================================================================
EOF!
}



