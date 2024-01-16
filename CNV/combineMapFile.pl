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

my (@input , $help , $outdir);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"o|outdir=s"	=>	\$outdir,
	"help"	=>	\$help,
);

if ($help or $#input<0 or ! $outdir){
	&help;
	exit;
}

my @head;
my %data;
for my $dir (@input){
	$dir = abs_path($dir);
	my @file  = </$dir/*.20K>;
	push @file , </$dir/*.20K.GC>;
	for my $file (@file){
		my $base = basename($file);
		@head = ();
		open IN , "$file";
		my $i = 0;
		while (<IN>){
			chomp;
			my ($chr , @d) = split /\t/ , $_;
			push @head , $chr;
			for my $j (0..$#d){
				$data{$base}->[$i]->[$j] += $d[$j];
			}
			$i++;
		}
		close IN;
	}
}

mkdir $outdir unless -d $outdir;
for my $base (sort keys %data){
	my @h = @head;
	open O , ">$outdir/$base";
	for my $d (@{$data{$base}}){
		my $h = shift @h;
		print O $h;
		for my $n (@{$d}){
			print O "\t$n";
		}
		print O "\n";
	}
	close O;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeMapFile.pl
#
#        USAGE: ./mergeMapFile.pl  
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
#      CREATED: 08/10/2018 12:42:57 PM
#     REVISION: ---
#===============================================================================
EOF!
}



