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
my $window = 1000000;
GetOptions(
	"i|input=s"	=>	\$input,
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
	$arm{$chr} = [$str , $end , $len];
}
close IN;

my $gc = readsWig("/share/work1/wangrr/DB/hg19/binCor/hg19.200K.gc.wig", 0);
my $mapab = readsWig("/share/work1/wangrr/DB/hg19/binCor/wgEncodeCrgMapabilityAlign100mer.200K.map.wig", 0);
my ($count , $total) = readsWig("$input", 1);

my %data;
my %undata;

my $tmpf = ".tmp.data.txt";
open D , ">$tmpf";
print D "arm\tcount\tgc\tmapab\n";
for my $arm (sort keys %$count){
	my @d = @{$count->{$arm}};
	@d = map{$_*(10**6)/$total} @d;
	if (exists $mapab->{$arm} and exists $gc->{$arm}){
		my $tn = 0;
		for my $i (0..$#d){
			if ($arm ne 'M'){
				next unless $gc->{$arm}->[$i] > 0.2 and $gc->{$arm}->[$i] < 0.8 and $mapab->{$arm}->[$i] > 0.9;
				$d[$i] = $d[$i]/$mapab->{$arm}->[$i];
				$tn++;
				if ($arm =~ /\d/){
					print D "$arm\t$d[$i]\t" , $gc->{$arm}->[$i] , "\t" , $mapab->{$arm}->[$i] , "\n";
				}else{
					$undata{$arm} += $d[$i];
				}
			}else{
				$undata{$arm} = $d[$i]/$mapab->{$arm}->[$i];
			}
		}
		if ($arm !~ /\d/ and $arm ne 'M'){
			$undata{$arm} = $undata{$arm} / $tn;
		}
	}else{
		my $td = 0;
		map{$td += $_} @d; 
		$undata{$arm} = $td / ($#d+1);
	}
}
close D;
my $xls = &r($tmpf);
my %tmph;
open IN , "$xls";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	next if $F[-1] eq 'NA';
	$tmph{$F[0]}->{t} += $F[-1];
	$tmph{$F[0]}->{n}++;
}
close IN;

print "ARM\tCount\n";
for my $arm (sort keys %undata){
	printf("%s\t%.3f\n" ,  ($arm , $undata{$arm}));
}
for my $arm (sort keys %tmph){
	printf("%s\t%.3f\n" , ($arm , $tmph{$arm}->{t}/$tmph{$arm}->{n}));
}


sub readsWig{
	my $pk = $_[1];
	open IN , "$_[0]";
	#fixedStep chrom=chr1 start=1 step=1000000 span=1000000
	my %loc;
	my $total = 0;
	my ($win , $chr , $ch);
	my %count;
	my $k = 0;
	while (<IN>){
		chomp;
		if (/fixedStep chrom=(\S+) start=\d+ step=(\d+) span=\d+/){
			($chr , $win) = ($1 , $2);
			$ch = $chr;
			$ch =~ s/chr//;
			if ($chr =~ /chr/ or $chr eq 'HBV_gtB' or $chr eq 'HBV_gtC'){
				$k = 1;
			}else{
				$k = 0;
			}
		}else{
			next if $k == 0;
			$total += $_;
			$loc{$chr} += $win;
			my $am = $ch;
			$_ = $_*(10**6)/$win if $pk;
			push @{$count{$am}} , $_;
		}
	}
	close IN;
	if ($pk){
		return (\%count , $total);
	}else{
		return \%count;
	}
}

sub r{
	my $in = $_[0];
	my $out = ".tmp.cor.xls";
	open RR , "| Rscript -";
	print RR << "EOF!";
d<-read.table('$in',header=TRUE)
Loess <- loess(d\$count ~ d\$gc)
loess.fittedRC <- predict(Loess, as.vector(d\$gc))
MedianRC <- median(d\$count, na.rm=T)
d\$cor<-d\$count + ( MedianRC - loess.fittedRC)
write.table(d, file='.tmp.cor.xls', quote=F, col.name=T, row.names=F, sep="\\t")
EOF!
	close RR;
	return $out;
}

sub rr{
	my $in = $_[0];
	my $out = ".tmp.cor.xls";
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
write.table(x, file='.tmp.cor.xls', quote=F, col.name=T, row.names=F, sep="\\t")
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



