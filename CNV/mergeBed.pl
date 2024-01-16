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
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %bed;
open IN , "$input";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	push @{$bed{$F[0]}} , [@F];
}
close IN;

for my $chr (sort keys %bed){
	my @bed = sort {$a->[1] <=> $b->[1]} @{$bed{$chr}};
	my ($c , $f , $r) = @{$bed[0]};
	shift @bed;
	for my $bed (@bed){
		my ($chr , $str , $end) = @$bed;
		if ($r < $str-1){
			print "$c\t$f\t$r\n";
			($c , $f , $r) = ($chr , $str , $end);
		}else{
			$r = $end if $end > $r;
		}
	}
	print "$c\t$f\t$r\n";
}


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeBed.pl
#
#        USAGE: ./mergeBed.pl  
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
#      CREATED: 02/23/18 14:58:16
#     REVISION: ---
#===============================================================================
EOF!
}



