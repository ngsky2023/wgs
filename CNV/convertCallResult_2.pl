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
my $cnv_fd_cutoff_rank_3 = $type eq 'cfDNA' ? 0.03 : 0.13;                 
my @cutsize = $type eq 'cfDNA' ? (20,30,50) : (10,20,50);
my $xn = $type eq 'cfDNA' ? 1 : 1;

my %arm;
open IN , "/share/work1/wangrr/DB/hg19/centromere.bed";
<IN>;
while (<IN>){
	chomp;
	my ($chr , $str , $end , $len) = split /\t/ , $_;
	$arm{$chr} = [$str , $end , $len];
}
close IN;

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
	chr16 => 1.5,
	chr17 => 1.5,
	chr18 => 1,
	chr19 => 4,
	chr20 => 1,
	chr21 => 1,
	chr22 => 1.5,
	chrX => 1.5,
	chrY => 1.5,
);

my %axs = (
	'10p' => 1,
	'10q' => 1,
	'11p' => 1,
	'11q' => 1,
	'12p' => 1,
	'12q' => 1,
	'13q' => 1,
	'14q' => 1,
	'15q' => 1,
	'16p' => 2.5,
	'16q' => 1.5,
	'17p' => 2,
	'17q' => 1.5,
	'18p' => 1,
	'18q' => 1,
	'19p' => 4,
	'19q' => 4,
	'1p' => 1,
	'1q' => 1,
	'20p' => 1,
	'20q' => 1,
	'21p' => 1,
	'21q' => 1,
	'22q' => 1.5,
	'2p' => 1,
	'2q' => 1,
	'3p' => 1,
	'3q' => 1,
	'4p' => 1,
	'4q' => 1,
	'5p' => 1.25,
	'5q' => 1,
	'6p' => 1.5,
	'6q' => 1,
	'7p' => 1,
	'7q' => 1,
	'8p' => 1,
	'8q' => 1,
	'9p' => 1,
	'9q' => 1,
	'Xp' => 1.5,
	'Xq' => 1.5,
	'Yp' => 1.5,
	'Yq' => 1.5,
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

		my ($pend , $qstr) = @{$arm{$chr}}[0,1];
		my $cute = ($pend + $qstr)/2;
		my $arm = $chr;
		$arm =~ s/chr//;
		if ($end < $cute){
			$arm .= 'p';
		}elsif ($str > $cute){
			$arm .= 'q';
		}else{
			$arm .= 'c';
		}

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
		my $xs = 1;
		if (exists $axs{$arm}){
			$xs = $axs{$arm};
		}else{
			$xs = $xs{$chr};
		}

		next if $nb <= 5 * $xs;
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




