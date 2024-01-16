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

my ($input , $help);
my $segm = 400;
my $gff = "/share/work1/wangrr/DB/hg19/hg19.gff";
my ($cutn , $totaln) = (2 , 4);
GetOptions(
	"i|input=s"	=>	\$input,
	"s|seg=s"	=>	\$segm,
	"g|gff=s"	=>	\$gff,
	"c|cutn=s"	=>	\$cutn,
	"t|totaln=s"	=>	\$totaln,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my $virusid = " HBV_gtA HBV_gtB HBV_gtC HBV_gtD HBV_gtE HBV_gtF HBV_gtG HBV_gtH";

#A00262:57:H37FJDSXX:3:2304:19081:28588  1201    HBV_gtB 2054    10      85S65M  HBV_gtC 2382    0       CCGCAGCCACCCTGCCG
my %v;
open IN , "samtools view  -F 1024 $input $virusid | grep -P 'chr\\w+' |";
while (<IN>){
	chomp;
	my ($id,$flag,$hbv,$hbvloc,$q,$mode,$chro , @bb) = split /\t/ , $_;
	next if $mode =~ /\d+[SH].*\d+[SH]/;
	next if $mode !~ /\d+[SH]/ and $chro !~ /^chr\w+$/;
	$flag = sprintf("%b" , $flag);
	my $R;
	if ($flag =~ /1\d{6}$/){
		$R = 1;
	}else{
		$R = 2;
	}
	my $vcut = $hbvloc;
	my $ccut;
	my $chr;
	my $maplen = &maplen($mode);
	my $k = 0;
	if ($mode =~ /\d+[SH]/){
		if ($mode =~ /\d+[SH]$/){
			$vcut += $maplen-1;
		}
		for my $bb (@bb){
			#SA:Z:chr13,35020442,+,46M104S,60,0;HBV_gtC,1002,-,36M114S,0,1;
			if ($bb =~ /^SA:Z:/){
				$bb =~ s/^SA:Z://;
				$bb =~ s/;$//;
				for my $cc (split /;/ , $bb){
					if ($cc =~ /^chr/){
						my @chrbb = split /,/ , $cc;
						$ccut = $chrbb[1];
						$chr = $chrbb[0];
						if ($chrbb[3] =~ /\d+[SH]$/){
							my $ml = &maplen($chrbb[3]);
							$ccut += $ml-1;
						}
						print STDERR "$id\t$cc\t$ccut\n";
						$k = 1;
						last;
					}
				}
				last;
			}
		}
	}
	if ($k == 1){
		#		$v{"$hbv\t$vcut\t$chr\t$ccut"}->{cut}++;
		$v{"$hbv\t$chr"}->{cut}->{$ccut}->{$vcut}++;
	}else{
		$chr = $chro;
		$ccut = $bb[0];
		if ($R == 1 and $flag =~ /0\d{4}$/){
			$vcut += $maplen;
		}
		#		$v{"$hbv\t$vcut\t$chr\t$ccut"}->{pair}++;
		$v{"$hbv\t$chr"}->{pair}->{$ccut}->{$vcut}++;
	}
}
close IN;

my $bp = 2;
my %vir;
for my $key (sort keys %v){
	next unless exists $v{$key}->{cut};
	my ($hbv , $chr) = split /\t/ , $key;
	my @ccut = sort {$a<=>$b} keys %{$v{$key}->{cut}};
	my @pos;
	my $i = 0;
	push @{$pos[$i]} , $ccut[0];
	for my $j (1..$#ccut){
		if ($ccut[$j] <= $pos[$i]->[-1] + $bp){
			push @{$pos[$i]} , $ccut[$j];
		}else{
			$i++;
			push @{$pos[$i]} , $ccut[$j];
		}
	}
	for my $pos (@pos){
		my @cut;
		my $cn = 0;
		for my $p (@$pos){
			for my $v (sort keys %{$v{$key}->{cut}->{$p}}){
				my $count = $v{$key}->{cut}->{$p}->{$v};
				push @cut , [$p , $v , $count];
				$cn += $count;
			}
		}
		@cut = sort {$b->[2]<=>$a->[2]} @cut;
		my ($ccut , $vcut , $count) = @{$cut[0]};
		my $pn = 0;
		for my $p (sort keys %{$v{$key}->{pair}}) {
			if ($p - $segm < $ccut and $ccut < $p + $segm){
				for my $pv (sort keys %{$v{$key}->{pair}->{$p}}){
					if ($pv - $segm < $vcut and $vcut < $pv + $segm){
						$pn += $v{$key}->{pair}->{$p}->{$pv};
					}
				}
			}
		}
		print STDERR "$hbv\t$vcut\t$chr\t$ccut\t$cn\t$pn\n";
		next if $cn < $cutn  or $cn+$pn < $totaln;
		push @{$vir{$chr}->{$ccut}} , [$hbv,$vcut,$chr,$ccut,$cn,$pn];
	}
}

my %gff;
my %z;
my %i;
undef %v;
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
				$i{"$chr\t$pos"} = $_;
			}
		}elsif ($F[2] eq 'gene'){
			#chr7    unknown gene    55086725        55275031        .       +       .       ID=EGFR;Name=EGFR;Name=EGFR
			next if exists $i{"$chr\t$pos"};
			my ($gene) = (/ID=([^;\s]+);/);
			if ($pos > $F[4]){
				$z{"$chr\t$pos"} = [$gene , $F[4]+1 , $pos-$F[4]-1];
			}else{
				if (exists $z{"$chr\t$pos"}){
					#HBV_gtC 81      chr1    120623773       0       2       intergenic      gene1=NOTCH2;gene2=FAM72B;dist1=11455;dist2=21523

					$v{"$chr\t$pos"} = [$z{"$chr\t$pos"}->[0] , $z{"$chr\t$pos"}->[2] , $gene , $F[3]-$pos-1];
					delete $z{"$chr\t$pos"};
				}
			}
		}
	}
}
close IN;

for my $key (sort keys %i){
	my @an = split /\n/ , $i{$key};
	my $k = 1;
	my $gene;
	my $zone = '';
	for my $ll (@an){
		if ($ll =~ /\tgene\t/ and $ll =~ /ID=([^;]+);/){
			$gene = $1;
		}elsif ($ll =~ /\tmRNA\t/ or $ll =~ /\tgene\t/ or $ll =~ /\ttranscript\t/){
			next;
		}else{
			$zone = (split /\t/ , $ll)[2];
			$k = 0;
		}
	}
	$zone = 'intron' if $k;
	$v{$key} = [$zone , $gene];
}

print "HBV\tHpos\tChr\tPos\tSplitReads\tPairReads\tDepth\tVAF\tRegion\tGene\n";
for my $chr (sort keys %vir){
	for my $pos (sort keys %{$vir{$chr}}){
		my @txt = @{$vir{$chr}->{$pos}};
		my $key = "$chr\t$pos";
		my @ann = @{$v{$key}};
		my $depth = depth($input , $chr , $pos);
		if ($#ann > 1){
			for my $txt (@txt){
				my $tn = $txt->[-1]+$txt->[-2];
				my $vaf = $tn/$depth;
				unless ($ann[0] =~ /TERT|CTNND2|CCNA2|CCNE1|KMT2B|MLL4|PTPRD|UNC5D|NRG3|AHRR|MKK4|TP53|CLPTM1L/ or $ann[2] =~ /TERT|CTNND2|CCNA2|CCNE1|KMT2B|MLL4|PTPRD|UNC5D|NRG3|AHRR|MKK4|TP53|CLPTM1L/ or ($txt->[-2] >= $cutn*2 and $tn >= $totaln*2)){
					next;
				}
				print join("\t" , @$txt) , "\t$depth\t$vaf\tintergenic\tgene1=$ann[0];gene2=$ann[2];dist1=$ann[1];dist2=$ann[3];\n";
			}
		}else{
			for my $txt (@txt){
				my $tn = $txt->[-1]+$txt->[-2];
				my $vaf = $tn/$depth;
				unless ($ann[1] =~ /TERT|CTNND2|CCNA2|CCNE1|KMT2B|MLL4|PTPRD|UNC5D|NRG3|AHRR|MKK4|TP53|CLPTM1L/ or ($txt->[-2] >= $cutn*2 and $tn >= $totaln*2)){
					next;
				}
				print join("\t" , @$txt) , "\t$depth\t$vaf\t$ann[0]\t$ann[1]\n";
			}
		}
	}
}


sub depth{
	my ($bam , $chr , $pos) = @_;
	my $d = readpipe("samtools depth $bam -r $chr:$pos-$pos");
	my ($depth) = ($d =~ /(\d+)$/);
	return $depth;
}



sub maplen{
	my ($mode) = @_;
	my $len = 0;
	while ($mode =~ /(\d+)([MSHID])/g){
		if ($2 eq 'M' or $2 eq 'D'){
			$len += $1;
		}
	}
	return $len;
}

sub htt{
	my $helpt =<< "EOF!";
#===============================================================================
#
#         FILE: bamVirus.pl
#
#        USAGE: ./bamVirus.pl  
#
#  DESCRIPTION: $0 -i bamfile > HBV.txt
#
#      OPTIONS: -c splitReadsNumber -t totalReadsNumber -g gff
#      VERSION: 1.0
#      CREATED: 03/14/18 15:28:16
#     REVISION: ---
#===============================================================================
EOF!
	return $helpt;
}
sub help{
	print &htt;
}


=head1 NAME
bamVirus

=head1 OPTIONS 
	-c	<splitReadsNumber>
	-t	<totalReadsNumber>
	-g	<gff>
	-s	<insert>

=cut
