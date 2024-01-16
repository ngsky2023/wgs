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
my $chr = '';
my $k = 'N';
GetOptions(
	"i|input=s"	=>	\$input,
	"c|chr=s"	=>	\$chr,
	"k|k=s"	=>	\$k,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my $q = 0;
$q = 10 if $k eq 'T';

open IN , "/share/public/software/samtools/1.4/samtools view -h -q $q -F 1024 $input $chr |" or die;
while (<IN>){
	if ($k eq 'N'){
		print $_;
		next;
	}

	if (/^\@/){
		print $_;
		next;
	}

	chomp;
	my @as=();
	my @xs=();
	my $nm = 0;
	my @F = split /\t/ , $_;
	for my $bb (@F[11..$#F]){
		if ($bb =~ /NM:i:(\d+)/){
			$nm = $1;
		}elsif ($bb =~ /^AS:/){
			@as=split(/:/,$bb);
		}elsif ($bb =~ /^XS:/){
			@xs=split(/:/,$bb);
		}
	}
	my ($len , $l , $s) = mode($F[5]);
	$nm -= $l;
	$nm += $s;
	my $diff=$as[2]-$xs[2];
	print "$_\n" if $len > 0 and $nm/$len <= 0.04 and $diff > 10;
}
close IN;

sub mode{
	my $m = $_[0];
	my ($len , $l , $s) = (0 , 0 , 0);
	while ($m =~ /(\d+)([MID])/g){
		if ($2 eq 'M'){
			$len += $1;
		}elsif ($2 eq 'I'){
			$len += $1;
			$l += $1-1;
		}elsif ($2 eq 'D'){
			$l += $1-1;
		}elsif ($2 eq 'S'){
			$s++;
		}
	}
	return ($len , $l , $s);
}




sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: filterMemBam.pl
#
#        USAGE: ./filterMemBam.pl  
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
#      CREATED: 03/26/18 15:50:28
#     REVISION: ---
#===============================================================================
EOF!
}



