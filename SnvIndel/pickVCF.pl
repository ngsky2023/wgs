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

my %pair;
open IN , "/share/Oncology/RD_HCC/IBFC2018028/WgsReSeq/newadd/pair.txt";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$pair{$F[0]} = $F[1];
}
close IN;

my %p;
open IN , "$ARGV[0]";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my $k = join("\t" , @F[1..5]);
	$p{$k} = 1;
}
close IN;

open IN , "$input";
my $sample = basename $ARGV[0];
$sample = (split /\./ , $sample)[0];
my $normal = $pair{$sample};
while (<IN>){
	if (/^##/){
		print $_;
		next;
	}elsif (/^#/){
		s/\tTUMOR\tNORMAL/\t$sample\t$normal/;
		print $_;
		next;
	}
	chomp;

	my @F = split /\t/ , $_;
	my ($chr , $str , $end , $ref , $alt) = @F[0..4];
	$end = length($ref) + $str - 1;
	for (my $i=0; $i<length($ref) and $i<length($alt); $i++){
		if (substr($ref, $i, 1) eq substr($alt, $i, 1)){
			$ref =~ s/^.//;
			$alt =~ s/^.//;
			$str++;
		}else{
			last;
		}
	}
	$ref = "-" if $ref eq "";
	$alt = "-" if $alt eq "";
	if ($ref eq "-"){
		$str--;
	}
	my $k = "$chr\t$str\t$end\t$ref\t$alt";
	print "$_\n" if exists $p{$k};
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: pickVCF.pl
#
#        USAGE: ./pickVCF.pl  
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
#      CREATED: 11/20/2019 12:39:23 PM
#     REVISION: ---
#===============================================================================
EOF!
}



