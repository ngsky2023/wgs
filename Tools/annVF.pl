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
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	my ($chr , $pos);
	if ($F[0] eq 'chrVirus'){
		($chr , $pos) = @F[3,4];
	}else{
		($chr , $pos) = @F[0,1];
	}
	$vir{$chr}->{$pos} = $_;
}
close IN;


my %gff;
my %z;
my %i;
open IN , "/share/work3/wangrr/DB/hg19/hg19.gff";
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
					#print $vir{$chr}->{$pos} , "$chr\t.\tintergenic\t" , $z{"$chr\t$pos"}->[1] , "\t" , $F[3]-1 , "\t.\t.\t.\tgene1=" , $z{"$chr\t$pos"}->[0] , ";gene2=$gene;dist1=" , $z{"$chr\t$pos"}->[2] , ";dist2=" , $F[3]-$pos-1 , ";\n\n";
					print $vir{$chr}->{$pos} , "\tintergenic\tgene1=" , $z{"$chr\t$pos"}->[0] , ";gene2=$gene;dist1=" , $z{"$chr\t$pos"}->[2] , ";dist2=" , $F[3]-$pos-1 , ";\n";
					delete $z{"$chr\t$pos"};
				}
			}
		}
	}
}
close IN;

for my $key (sort keys %i){
	#print $i{$key} , "\n";
	my ($bk , @an) = split /\n/ , $i{$key};
	my $k = 1;
	my $gene;
	my $zone = '';
	for my $ll (@an){
		if ($ll =~ /\tgene\t/ and $ll =~ /ID=([^;]+);/){
			$gene = $1;
		}elsif ($ll =~ /\tmRNA\t/ or $ll =~ /\tgene\t/){
			next;
		}else{
			$zone = (split /\t/ , $ll)[2];
			print "$bk\t$zone\t$gene\n";
			$k = 0;
		}
	}
	print "$bk\tintron\t$gene\n" if $k;
}

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



