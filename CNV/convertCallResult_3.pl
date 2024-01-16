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

my %back;
my $fn = 0;
for my $file (<$Bin/back52/XJ*.data.txt>){
	open IN , "$file";
	while (<IN>){
		chomp;
		my ($chr , @F) = split /\t/ , $_;
		$back{"chr$chr"}->[$fn] = [@F];
	}
	close IN;
	$fn++;
}

my %ref;
open IN , "$Bin/reference_cfDNA52_nB5_nS5";
while (<IN>){
	chomp;
	my ($chr , @F) = split /\t/ , $_;
	$ref{"chr$chr"} = [@F];
}
close IN;

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
		my ($tn , $sn , $en , $rn);
		if ($bininfo =~ /(\d+)\[(\d+)-(\d+)\((\d+) NA\]/){
			($tn , $sn , $en , $rn) = ($1 , $2 , $3 , $4);
			$nn = $tn - $rn;
			$sn--;
			$en--;
			next if $sn > $en;
		}else{
			print STDERR "bb\t$_\n";
			next;
		}

		my $xs = 1;
		if (exists $axs{$arm}){
			$xs = $axs{$arm};
		}
		my $nb = $nn*$change;
		next if $nb <= 5 * $xs;

		my @back;
		for my $dl (@{$back{$chr}}){
			my $dt = 0;
			map{$dt += $_ if $_ ne 'NA' and $_ > 0} @{$dl}[$sn..$en];
			push @back , $dt;
		}
		@back = sort {$a<=>$b} @back;
		my @ref = @{$ref{$chr}}[$sn..$en];
		my $mean = 0;
		map{$mean += $_ if $_ ne 'NA' and $_ > 0} @ref;
		my $data = $copy/$normal * $mean;

		my $stat = Statistics::Descriptive::Full->new();
		my @tmp = @back;
		@back = ();
		map{push @back , $_ if $_ ne 'NA' and $_ > 0} @tmp;
		$stat->add_data(\@back);
		#my $mean = $stat->mean();
		my $sd = $stat->standard_deviation();
		print STDERR "$_\n" unless $sd;

		my $bound = 4 * $xs;
		my $lcut = $mean - $sd * $bound;
		my $rcut = $mean + $sd * $bound;
		$cv = $sd/$mean;

		my $out = "$sample\t$chr\t$str\t$end\t$copy\t$class\t$lenm\t$cv\t$data\t$back[0]\t$lcut\t$mean\t$rcut\t$back[-1]\n";

		#next if $nb <= 5 * $xs;
		#next if $nn/$tn < 0.4;
		#if ($nn > $cutsize[2] and $change > $cnv_fd_cutoff_rank_3*$xn*$xs){
		#	print $out;
		#}elsif ($nn > $cutsize[1] and $nn <= $cutsize[2] and $change > $cnv_fd_cutoff_rank_2*$xn*$xs){
		#	print $out;
		#}elsif ($nn > $cutsize[0] and $nn <= $cutsize[1] and $change > $cnv_fd_cutoff_rank_1*$xn*$xs){
		#	print $out;
		#}else{
		#	print STDERR "$out";
		#}
		#$rcut = $rcut > $back[-1] ? $rcut : $back[-1];
		#$lcut = $lcut < $back[0] ? $lcut : $back[0];
		if ((($class eq 'dup' and $data > $rcut) or ($class eq 'del' and $data < $lcut)) and $lenm > 2){
			print $out;
		}else{
			print STDERR $out;
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




