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

my ($input , $help , $meerkat , $manta);
GetOptions(
	"i|input=s"	=>	\$input,
	"meerkat=s"	=>	\$meerkat,
	"manta=s"	=>	\$manta,
	"help"	=>	\$help,
);

if ($help){
	&help;
	exit;
}

my %zone;
my %sv;
my %mk;
for my $file ($manta , $meerkat){
	open IN , "$file";
	while (<IN>){
		if (/^#/){
			next;
		}
		chomp;
		my @F = split /\t/ , $_;
		my $chr = $F[0];
		my $str = $F[1];
		if ($F[4] =~ /DUP|INV|DEL/){
			if ($F[7] =~ /^END=(\d+);/ or $F[7] =~ /;END=(\d+);/){
				my $end = $1;
				($str , $end) = sort {$a<=>$b} ($str , $end);
				my $k = '';
				my @zone = sort keys %zone;
				for my $z (@zone){
					my ($zf , $zr) = sort {$a<=>$b} (split /\t/ , $z);
					next if $zr <= $str or $end <= $zf;
					my ($z1 , $z2 , $z3 , $z4) = sort {$a<=>$b} ($zf , $zr , $str , $end);
					if (($z3-$z2)/($z4-$z1) > 0.8){
						$k = $zone{$z};
						$zone{"$str\t$end"} = $k;
						$sv{$k} .= "\t$_";
						$mk{$k} .= ":$file";
						last;
					}
				}
				unless ($k){
					$zone{"$str\t$end"} = $F[2];
					$sv{$F[2]} = $_;
					$mk{$F[2]} = $file;
				}
			}
		}
	}
}
close IN;

for my $key (sort keys %sv){
	my $ll = $sv{$key};
	my $mk = $mk{$key};
	print "$mk\t$ll\n";
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeSV.pl
#
#        USAGE: ./mergeSV.pl  
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
#      CREATED: 12/11/17 15:29:37
#     REVISION: ---
#===============================================================================
EOF!
}



