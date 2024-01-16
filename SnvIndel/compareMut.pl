#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use Data::Dumper;
use lib $Bin;

my ($gatk , $strelka , $help , $prx , $bed);
GetOptions(
	"g|gatk=s"	=>	\$gatk,
	"s|strelka=s"	=>	\$strelka,
	"p|prx=s"	=>	\$prx,
	"b|bed=s"	=>	\$bed,
	"help"	=>	\$help,
);

my %bed;
if ($bed){
	open IN , "$bed";
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		$F[0] = "chr$F[0]";
		push @{$bed{$F[0]}} , [@F[1,2]];
	}
	close IN;
}


my %gatk;
my (%gn , %sn , %same);
open IN , "$gatk";
<IN>;
#print "$h";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	if ($bed){
		my $k = 0;
		for my $ll (@{$bed{$F[1]}}){
			if ($ll->[0] <= $F[2] and $F[2] <= $ll->[1]){
				$k = 1;
				last;
			}
		}
		next if $k == 0;
	}
	my $k = join("\t" , @F[1..5]);
	#print STDERR "$k\n";
	#print "G\t$k\n";
	$gatk{$k} = "$_\n";
}
close IN;

open G , ">$prx.gatk.xls";
open S , ">$prx.strelka.xls";
open IN , "$strelka";
my $h = <IN>;
print $h;
while (<IN>){
	next if /^#/;
	chomp;
	my @F = split /\t/ , $_;
	my ($chr , $str , $end , $ref , $alt) = @F[1..5];

	if ($bed){
		my $k = 0;
		for my $ll (@{$bed{$chr}}){
			if ($ll->[0] <= $str and $end <= $ll->[1]){
				$k = 1;
				last;
			}
		}
		next if $k == 0;
	}

	my $k = "$chr\t$str\t$end\t$ref\t$alt";
	my $type = type($ref , $alt);
	if (exists $gatk{$k}){
		$same{$type}++;
		print "$_\n";
		delete $gatk{$k};
	}else{
		print S "$_\n";
		$sn{$type}++;
	}
}
close IN;

for my $k (keys %gatk){
	my $v = $gatk{$k};
	my ($chr , $str , $end , $ref , $alt) = split /\t/ , $k;
	my $type = type($ref , $alt);
	print G "$v";
	$gn{$type}++;
}

#for my $t (qw/snv indel/){
#	print $t;
#	my ($s , $to) = (0 , 0);
#	if (exists $same{$t}){
#		print "\t" , $same{$t};
#		$s += $same{$t};
#		$to += $same{$t};
#	}else{
#		print "\t0";
#	}
#	if (exists $gn{$t}){
#		print "\t" , $gn{$t};
#		$to += $gn{$t};
#	}else{
#		print "\t0";
#	}
#	if (exists $sn{$t}){
#		print "\t" , $sn{$t};
#		$to += $sn{$t};
#	}else{
#		print "\t0";
#	}
#	print "\t" , $s/$to;
#	print "\n";
#}
#print "$same\t$gn\t$sn\n";


sub type{
	my ($ref , $alt) = @_;
	if (length($ref) > 1 or length($alt) > 1 or $alt eq "-" or $ref eq "-"){
		return "indel";
	}else{
		return "snv";
	}
}



sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: compare.pl
#
#        USAGE: ./compare.pl  
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
#      CREATED: 11/05/2019 02:19:59 PM
#     REVISION: ---
#===============================================================================
EOF!
}



