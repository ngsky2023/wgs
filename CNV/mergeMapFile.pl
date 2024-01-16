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

if ($help or ! $#input<0){
	&help;
	exit;
}

my @data;
my @head;
for my $file (@input){
	@head = ();
	open IN , "$file";
	my $i = 0;
	while (<IN>){
		chomp;
		my ($chr , @d) = split /\t/ , $_;
		push @head , $chr;
		for my $j (0..$#d){
			$data[$i]->[$j] += $d[$j];
		}
		$i++;
	}
	close IN;
}

for my $d (@data){
	my $h = shift @head;
	print $h;
	for my $n (@{$d}){
		print "\t$n";
	}
	print "\n";
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



