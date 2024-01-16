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

my (@input , $help , $type);
$type = 'cfDNA';
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"t|type=s"	=>	\$type,
	"help"	=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

my $cnv_fd_cutoff_rank_1 = $type eq 'cfDNA' ? 0.15 : 0.3;                  
my $cnv_fd_cutoff_rank_2 = $type eq 'cfDNA' ? 0.07 : 0.27;                 
my $cnv_fd_cutoff_rank_3 = $type eq 'cfDNA' ? 0.02 : 0.13;                 
my @cutsize = $type eq 'cfDNA' ? (10,15,25) : (10,20,50);
my $xn = $type eq 'cfDNA' ? 1 : 1;
my %xs = (
	chr1 => 1,
	chr2 => 1,
	chr3 => 1,
	chr4 => 1,
	chr5 => 1,
	chr6 => 1,
	chr7 => 1,
	chr8 => 1,
	chr9 => 1,
	chr10 => 1,
	chr11 => 1,
	chr12 => 1,
	chr13 => 1,
	chr14 => 1,
	chr15 => 1,
	chr16 => 2,
	chr17 => 2,
	chr18 => 1,
	chr19 => 5,
	chr20 => 1,
	chr21 => 1,
	chr22 => 2,
	chrX => 2,
	chrY => 2,
);

for my $file (@input){
	open IN , "$file";
	my $head = <IN>;
	my ($sample , $sex) = split /\t/ , $head;
	while (<IN>){
		chomp;
		my (undef , $loc , $copy , $bininfo , $cv) = split /\t/ , $_;
		my ($chr , $str , $end) = split /[:-]/ , $loc;
		next if $chr =~ /[XY]/;
		my $len = $end-$str;
		my $lenm = sprintf("%.2f" , $len/1000000);
		next if $lenm < 2;
		my $class = '';
		my $normal = 2;
		if ($sex eq 'M' and $chr =~ /[XY]/){
			$normal = 1;
		}elsif ($sex eq 'F' and $chr eq 'chrY'){
			next;
		}
		next if $copy < 0;
		$cv =~ s/\(CV\)//;
		my $change = abs($copy-$normal);
		if ($copy > $normal){
			$class = 'dup';
		}elsif ($copy < $normal){
			$class = 'del';
		}else{
			next;
		}
	
		my $nn = 0;
		#258[1625-1882(0 NA]
		if ($bininfo =~ /(\d+)\[(\d+)-(\d+)\((\d+) NA\]/){
			my ($tn , $sn , $en , $rn) = ($1 , $2 , $3 , $4);
			$nn = $tn - $rn;
		}else{
			print STDERR "bb\t$_\n";
			next;
		}
		my $nb = $nn*$change;
		my $out = "$sample\t$chr\t$str\t$end\t$copy\t$class\t$lenm\t$cv\t$nb\t$nn\n";
		my $xs = $xs{$chr};
		next if $nb < 5 * $xs;
		if ($nn > $cutsize[2] and $change > $cnv_fd_cutoff_rank_3*$xn*$xs){
			print $out;
		}elsif ($nn > $cutsize[1] and $nn <= $cutsize[2] and $change > $cnv_fd_cutoff_rank_2*$xn*$xs){
			print $out;
		}elsif ($nn > $cutsize[0] and $nn <= $cutsize[1] and $change > $cnv_fd_cutoff_rank_1*$xn*$xs){
			print $out;
		}else{
			print STDERR "$out";
		}
		
	}
	close IN;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: convertResult.pl
#
#        USAGE: ./convertResult.pl  
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
#      CREATED: 06/11/18 13:26:48
#     REVISION: ---
#===============================================================================
EOF!
}



