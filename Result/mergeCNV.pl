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

if ($help){
	&help;
	exit;
}

my %mustcnv;
open IN , "$Bin/CNVgene.txt";
while (<IN>){
        chomp;
        my ($gene , $t) = split /\t/ , $_;
        $mustcnv{$gene} = $t;
}
close IN;

my %cnv;
my %tumor;
for my $file (@input){
	open IN , "$file";
	my $h = <IN>;
	chomp $h;
	my @h = split /\t/ , $h;
	my $sample = $h[4];
	$sample =~ s/_T$/-T/;
	$sample =~ s/_P$/-P/;
	$sample =~ s/_X$/-X/;
	$tumor{$sample} = 1;
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		next if $F[4] eq '0' or $F[0] =~ /^\d+$/;
		my ($type , $num) = split /\|/ , $F[4];
		next unless exists $mustcnv{$F[0]};
		next if 0.5 < $num and $num < 4;
		if ($type == 1){
			$type = "amplification";
		}else{
			$type = "deletion";
		}
		$cnv{$sample}->{$F[0]} = [$sample , $F[0] , $type , $num];
	}
	close IN;

}

print  "sample\tgene\tvariationType\tcopyNumber\n";
for my $sample (sort keys %tumor){
	if (exists $cnv{$sample}){
		for my $key (sort keys %{$cnv{$sample}}){
			print  join("\t" , @{$cnv{$sample}->{$key}}) , "\n";
		}
	}else{
		print "$sample\tNegative\t-\t-\n";
	}
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeCNV.pl
#
#        USAGE: ./mergeCNV.pl  
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
#      CREATED: 11/05/2018 09:51:06 AM
#     REVISION: ---
#===============================================================================
EOF!
}



