#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib "/share/work1/wangrr/local/Plib/lib/perl5/";
use Statistics::Descriptive;
use lib $Bin;

my (@input , $help);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"help"	=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

my %arm;
open IN , "/share/work1/wangrr/DB/hg19/centromere.bed";
<IN>;
#chrY    10104553        13104553        59373566
#chrM    0       1       16571
while (<IN>){
	chomp;
	my ($chr , $str , $end , $len) = split /\t/ , $_;
	$arm{$chr} = [$str , $end , $len];
}
close IN;


my %table;
my @sample = ();

for my $file (@input){
	my $base = basename $file;
	my ($sample) = ($base =~ /^([^\.]+)/);
	push @sample , $sample;

	open IN , "$file";
	#fixedStep chrom=chr1 start=1 step=1000000 span=1000000
	#725970
	my %loc;
	my ($win , $chr , $ch);
	my %count;
	my $total = 0;
	while (<IN>){
		chomp;
		if (/fixedStep chrom=(\S+) start=\d+ step=(\d+) span=\d+/){
			($chr , $win) = ($1 , $2);
			$ch = $chr;
			$ch =~ s/chr//;
		}else{
			next unless exists $arm{$chr};
			$total += $_;
			$loc{$chr} += $win;
			if ($loc{$chr} < $arm{$chr}->[0]){
				my $am = $ch;
				if ($am ne 'chrM'){
					$am = $am."p";
				}
				$count{$am} += $_;
			}elsif ($loc{$chr} > $arm{$chr}->[1]){
				my $am = $ch;
				if ($am ne 'M'){
					$am = $am."q";
				}
				$count{$am} += $_;
			}
		}
	}
	close IN;

	for my $am (sort keys %count){
		my $P = $count{$am}/$total;
		$table{$am}->{$sample} = $P;
	}
}

print "ARM\tmean\tSD\t" , join("\t" , @sample) , "\n";
for my $arm (sort keys %table){
	my @data = values %{$table{$arm}};
	my ($mean , $sd);
	my $stat = Statistics::Descriptive::Full->new();
	@data = sort {$a<=>$b} @data;
	my $md = $data[abs($#data/2)];
	$stat->add_data(\@data);
	$mean = $stat->mean();
	next if $mean == 0 or $md == 0;
	$sd = $stat->standard_deviation();
	print "$arm\t$md\t$sd";
	for my $sample (@sample){
		print "\t" , $table{$arm}->{$sample};
	}
	print "\n";
}




sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeWindow2Arm.pl
#
#        USAGE: ./mergeWindow2Arm.pl  
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
#      CREATED: 04/03/18 14:16:43
#     REVISION: ---
#===============================================================================
EOF!
}



