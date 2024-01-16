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

my (@input , $help);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"help"	=>	\$help,
);

if ($help or $#input<0){
	&help;
	exit;
}

my %pos;
for my $n (</share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/result/sv/BXB/T/*/results/variants/tumorSV.vcf.gz>){
	open IN , "zcat $n|";
	while (<IN>){
		next if /^#/;
		chomp;
		my @F = split /\t/ , $_;
		$pos{"$F[0]:$F[1]"}++;
		if (/END=(\d+)/ and $1 != $F[1]){
			$pos{"$F[0]:$1"}++;
		}
	}
	close IN;
}

for my $file (@input){
	my %mate;
	my $out = basename $file;
	$out =~ s/\..*$//;
	$out = "$out.filter.vcf";
	open O , ">$out";
	open IN , "$file";
	while (<IN>){
		if (/^#/){
			print O $_;
			next;
		}
		chomp;
		my @F = split /\t/ , $_;
		my $k1 = "$F[0]:$F[1]";
		if (/SVTYPE=BND;/){
			my $id = $F[2];
			$id =~ s/:\d+$//;
			if (exists $mate{$id}){
				if ($F[0] =~ /_/ or $mate{$id}->[0] =~ /_/){
					next;
				}
				my $k2 = $mate{$id}->[0].":".$mate{$id}->[1];
				my $kn = 0;
				if ((exists $pos{$k1} and $pos{$k1}>1 and $kn=$pos{$k1}) or (exists $pos{$k2} and $pos{$k2}>1 and $kn=$pos{$k2})){
					print STDERR "$out\t$kn\t",join("\t" , @{$mate{$id}}) , "\t$_\n";
					next;
				}else{
					print O join("\t" , @{$mate{$id}}) , "\n$_\n";
				}
			}else{
				$mate{$id} = [@F];
			}
		}else{
			if ($F[0] =~ /_/){
				next;
			}
	
			if (exists $pos{$k1} and $pos{$k1} > 1){
				print STDERR "$out\t$pos{$k1}\t$_\n";
				next;
			}else{
				print O "$_\n";
			}
		}
	}
	close O;
	close IN;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: filterPosN.pl
#
#        USAGE: ./filterPosN.pl  
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
#      CREATED: 02/21/2020 05:24:44 PM
#     REVISION: ---
#===============================================================================
EOF!
}



