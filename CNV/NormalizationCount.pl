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

my ($input , $help , $prx);
my $window = 200;
GetOptions(
	"i|input=s"	=>	\$input,
	"w|w=s"	=>	\$window,
	"p|prx=s"	=>	\$prx,
	"help"	=>	\$help,
);

if ($help or ! $input){
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
	$chr =~ s/chr//;
	$arm{$chr} = [$str , $end , $len];
}
close IN;

my %exc;
open IN , "/share/work1/wangrr/DB/hg19/hg19Excludable.bed";
while (<IN>){
	chomp;
	my ($chr , $str , $end , @m) = split /\t/ , $_;
	$chr =~ s/chr//;
	push @{$exc{$chr}} , [$str , $end];
}
close IN;

my ($gc , $mapab , $count , $win);
($gc , $win) = readsWig("/share/work1/wangrr/DB/hg19/binCor/hg19.${window}K.gc.wig");
($mapab , $win) = readsWig("/share/work1/wangrr/DB/hg19/binCor/wgEncodeCrgMapabilityAlign100mer.${window}K.map.wig");
($count , $win) = readsWig($input);

my %data;
my %undata;

my $tmpf = "$prx.data.txt";
open D , ">$tmpf";
my $total = 0;
print D "chr\tpos\tlen\tcount\tfilter\tgc\tmapab\n";
for my $chr (sort keys %$count){
	if (exists $mapab->{$chr} and exists $gc->{$chr}){
		my @d = @{$count->{$chr}};
		my $loc = -1 * $win;
		for my $i (0..$#d){
			$loc += $win;

			my $len = $win;
			if ($loc + $win > $arm{$chr}->[2]){
				$len = $arm{$chr}->[2]-$loc;
			}

			my $mk = 0;
			for my $zone (@{$exc{$chr}}){
				my ($str , $end) = @$zone;
				unless ($loc+$win < $str or $end < $loc){
					$mk = 1;
					last;
				}
			}

			unless ($gc->{$chr}->[$i] > 0.31 and $gc->{$chr}->[$i] < 0.55 and $mapab->{$chr}->[$i] > 0.9){
				$mk = 1;
			}
			if ($mk == 0){
				$total += $d[$i];
			}else{
#				next;
			}
			print D "$chr\t$loc\t$len\t$d[$i]\t$mk\t" , $gc->{$chr}->[$i] , "\t" , $mapab->{$chr}->[$i] , "\n";
		}
	}
}
close D;
my $xls = &GCcor($tmpf, $total);
#unlink $tmpf;
my %tmph;
my %out;
open IN , "$xls";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my $loc = "$F[0]\t$F[1]";
	push @{$out{$F[0]}} , $F[-1];
	next if $F[-1] eq 'NA';
#	$tmph{$F[0]}->{t} += $F[-1];
	$tmph{$loc} = $F[-1];
}
close IN;
#unlink $xls;


open O , ">$prx.count.txt";
print O "Chr\tPos\tCount\n";
for my $loc (sort keys %tmph){
	print O sprintf("%s\t%.3f\n" , ($loc , $tmph{$loc}));
}
close O;

my $tn = 12500 / int($window/20);
open O , ">$prx.seqin.txt";
for my $chr (1..24){
	my $n;
	if (exists $out{$chr}){
		print O "$chr\t" , join("\t" , @{$out{$chr}});
		$n = @{$out{$chr}};
	}else{
		if ($chr eq '23'){
			print O "$chr\t" , join("\t" , @{$out{X}});
			$n = @{$out{X}};
		}elsif ($chr eq '24'){
			print O "$chr\t" , join("\t" , @{$out{Y}});
			$n = @{$out{Y}};
		}
	}
	print O "\tNA" x ($tn-$n) , "\n";
}
close O;


sub readsWig{
	open IN , "$_[0]";
	#fixedStep chrom=chr1 start=1 step=1000000 span=1000000
	my %loc;
	my ($win , $chr , $ch);
	my %count;
	my $k = 0;
	my $mi = 0;
	my $am;
	while (<IN>){
		chomp;
		if (/fixedStep chrom=(\S+) start=\d+ step=(\d+) span=\d+/){
			$chr = $1;
			my $tmpwin = $2;
			$win = $tmpwin if $chr =~ /chr/;
			$ch = $chr;
			$ch =~ s/chr//;
			if ($chr =~ /chr[\dMXY]+/ or $chr =~ /HBV/){
				$k = 1;
			}else{
				$k = 0;
			}
		}else{
			next if $k == 0;
			push @{$count{$ch}} , $_;
		}
	}
	close IN;
	return (\%count , $win);
}

sub r{
	my $in = $_[0];
	my $out = ".$input.tmp.cor.xls";
	open RR , "| Rscript -";
	print RR << "EOF!";
d<-read.table('$in',header=TRUE)
Loess <- loess(d\$count ~ d\$gc)
loess.fittedRC <- predict(Loess, as.vector(d\$gc))
MedianRC <- median(d\$count, na.rm=T)
d\$cor<-d\$count + ( MedianRC - loess.fittedRC)
write.table(d, file='$out', quote=F, col.name=T, row.names=F, sep="\\t")
EOF!
	close RR;
	return $out;
}

sub GCcor{
	my ($in , $total) = @_;
	my $out = "$prx.cor.xls";
	print readpipe("Rscript $Bin/correctionCount.R $in $out $total $input ./");
	return $out;
}

sub rr{
	my $in = $_[0];
	my $out = ".$input.tmp.cor.xls";
	open RR , "| Rscript -";
	print RR << "EOF!";
x<-read.table('$in',header=TRUE)
count <- x\$count
gc <- x\$gc
mapab <- x\$mapab
MedianRC <- median(count, na.rm=T)
rough = loess(count ~ gc, span = 0.03)
i <- seq(0, 1, by = 0.001)
final = loess(predict(rough, i) ~ i, span = 0.3)
cor.gc <- count * MedianRC / predict(final, gc)
#MedianRC <- median(cor.gc, na.rm=T)
#final = approxfun(lowess(mapab, cor.gc))
#cor.mapab <- cor.gc * MedianRC / final(mapab)
x\$cor <- cor.gc
write.table(x, file='$out', quote=F, col.name=T, row.names=F, sep="\\t")
EOF!
	close RR;
	return $out;
}





sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: NormalizationCountCNV.pl
#
#        USAGE: ./NormalizationCountCNV.pl  
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
#      CREATED: 04/11/18 11:32:25
#     REVISION: ---
#===============================================================================
EOF!
}



