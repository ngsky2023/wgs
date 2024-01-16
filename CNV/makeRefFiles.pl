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

my (@input , $help , $ref , $dir);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"r|ref=s"	=>	\$ref,
	"d|dir=s"	=>	\$dir,
	"help"	=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

my %ref;
open IN , "$ref";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	last if $F[0] > 22;
	my $chr = shift @F;
	$ref{$chr} = [@F];
}
close IN;

mkdir $dir unless -d $dir;
for my $file (@input){
	open IN , "$file";
	my $base = basename $file;
	open O , ">$dir/$base";
	while (<IN>){
		chomp;
		my ($chr , @F) = split /\t/ , $_;
		last if $chr > 22;
		print O $chr;
		my @ref = @{$ref{$chr}};
		for my $i (0..$#F){
			if ($i <= $#ref and $F[$i] ne 'NA'){
				print O "\t" , $F[$i] * $ref[$i];
			}else{
				print O "\t" , $F[$i];
			}
		}
		print O "\n";
	}
	close IN;
	close O;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: makeRefFiles.pl
#
#        USAGE: ./makeRefFiles.pl  
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
#      CREATED: 08/27/2018 01:31:09 PM
#     REVISION: ---
#===============================================================================
EOF!
}



