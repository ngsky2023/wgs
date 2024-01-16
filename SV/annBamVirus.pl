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

#chr13	81515749	+	chrVirus	2188	+	0+15	high
my %vir;
open IN , "$input";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my ($chr , $pos);
	($chr , $pos) = @F[2,3];
	$vir{$chr}->{$pos} = $_;
}
close IN;


my %gff;
my %z;
my %i;
my %v;
open IN , "/share/work1/wangrr/DB/hg19/hg19.gff";
while (<IN>){
	next if /^#/;
	my @F = split /\t/ , $_;
	my $chr = $F[0];
	next unless exists $vir{$chr};
	for my $pos (sort keys %{$vir{$chr}}){
		if ($F[3] <= $pos and $pos <= $F[4]){
			if (exists $i{"$chr\t$pos"}){
				$i{"$chr\t$pos"} .= $_;
			}else{
				$i{"$chr\t$pos"} = $vir{$chr}->{$pos} . "\n" . $_;
			}
		}elsif ($F[2] eq 'gene'){
			#chr7    unknown gene    55086725        55275031        .       +       .       ID=EGFR;Name=EGFR;Name=EGFR
			next if exists $i{"$chr\t$pos"};
			my ($gene) = (/ID=([^;\s]+);/);
			if ($pos > $F[4]){
				$z{"$chr\t$pos"} = [$gene , $F[4]+1 , $pos-$F[4]-1];
			}else{
				if (exists $z{"$chr\t$pos"}){
					#print $vir{$chr}->{$pos} , "\tintergenic\tgene1=" , $z{"$chr\t$pos"}->[0] , ";gene2=$gene;dist1=" , $z{"$chr\t$pos"}->[2] , ";dist2=" , $F[3]-$pos-1 , ";\n";
					my @vv = split /\t/ , $vir{$chr}->{$pos};
					#HBV_gtC 81      chr1    120623773       0       2       intergenic      gene1=NOTCH2;gene2=FAM72B;dist1=11455;dist2=21523
					if (exists $v{$z{"$chr\t$pos"}->[0].";$gene"} and exists $v{$z{"$chr\t$pos"}->[0].";$gene"}->{$vv[0]}){
						$v{$z{"$chr\t$pos"}->[0].";$gene"}->{$vv[0]}->{cut} += $vv[4];
						$v{$z{"$chr\t$pos"}->[0].";$gene"}->{$vv[0]}->{pair} += $vv[5];
					}else{
						$v{$z{"$chr\t$pos"}->[0].";$gene"}->{$vv[0]}->{cut} = $vv[4];
						$v{$z{"$chr\t$pos"}->[0].";$gene"}->{$vv[0]}->{pair} = $vv[5];
						$v{$z{"$chr\t$pos"}->[0].";$gene"}->{$vv[0]}->{txt} = [@vv[0..3],'intergenic', "gene1=".$z{"$chr\t$pos"}->[0].";gene2=$gene;dist1=".$z{"$chr\t$pos"}->[2].";dist2=".($F[3]-$pos-1).";"];
					}
					delete $z{"$chr\t$pos"};
				}
			}
		}
	}
}
close IN;

for my $key (sort keys %i){
	my ($bk , @an) = split /\n/ , $i{$key};
	my $k = 1;
	my $gene;
	my $zone = '';
	my @vv = split /\t/ , $bk;
	for my $ll (@an){
		if ($ll =~ /\tgene\t/ and $ll =~ /ID=([^;]+);/){
			$gene = $1;
		}elsif ($ll =~ /\tmRNA\t/ or $ll =~ /\tgene\t/ or $ll =~ /\ttranscript\t/){
			next;
		}else{
			$zone = (split /\t/ , $ll)[2];
			#print "$bk\t$zone\t$gene\n";
			$k = 0;
		}
	}
	$zone = 'intron' if $k;
	if (exists $v{"$gene:$zone"} and exists $v{"$gene:$zone"}->{$vv[0]}){
		$v{"$gene:$zone"}->{$vv[0]}->{cut} += $vv[4];
		$v{"$gene:$zone"}->{$vv[0]}->{pair} += $vv[5];
	}else{
		$v{"$gene:$zone"}->{$vv[0]}->{cut} = $vv[4];
		$v{"$gene:$zone"}->{$vv[0]}->{pair} = $vv[5];
		$v{"$gene:$zone"}->{$vv[0]}->{txt} = [@vv[0..3] , $zone , $gene];
	}
}

for my $key (sort keys %v){
	for my $vv (sort keys %{$v{$key}}){
		print join("\t" , (@{$v{$key}->{$vv}->{txt}}[0..3], $v{$key}->{$vv}->{cut}, $v{$key}->{$vv}->{pair}, @{$v{$key}->{$vv}->{txt}}[4,5])) , "\n";
	}
}
#HBV_gtC 2912    chr13   34871364        0       2       intergenic      gene1=RFC3;gene2=NBEA;dist1=330668;dist2=645059;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: annVF.pl
#
#        USAGE: ./annVF.pl  
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
#      CREATED: 11/06/17 10:21:35
#     REVISION: ---
#===============================================================================
EOF!
}



