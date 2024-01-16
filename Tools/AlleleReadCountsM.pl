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

my ($input , $help , $bam , $bed);
GetOptions(
	"i|input=s"	=>	\$input,
	"b|bam=s"	=>	\$bam,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

print "chr\tposition\tref\trefCount\tNref\tNrefCount\n";
open IN , "$input";
while (<IN>){
	my $ll = $_;
	chomp;
	next if /^#/;
	my ($chr , $pos , $ref , $alt , $tinfo) = (split /\t/ , $_)[0,1,3,4,9];
	if ($alt =~ /,/){
		next;
	}

	if ($chr !~ /chr/){
		$chr = "chr$chr";
	}
	my $loc = $pos;
	my $type = '';

	my $rl = length($ref);
	my $al = length($alt);

	if ($rl == $al and $rl == 1){
		$type = 'snp';
	}else{
		next;
	}

	my %dd;
	my $depth = 0;
	my $altn = 0;
	my $refnn = 0;
	for my $file (($bam)){
		open BAM , "-|" , "samtools view -q 10 -F 1024 $file $chr\:$pos-$pos";
		while (my $sam = <BAM>){
			my @sp = split /\t/ , $sam;

			my ($readn , $refn) = (-1 , $sp[3]-1);

			$depth++;
			my $mn = 0;
			while ($sp[5] =~ /(\d+)([IDMS])/g){
				$mn++;
				my ($d , $w) = ($1 , $2);
				if ($w eq 'S'){
					$readn += $d;
				}elsif ($w eq 'I'){
					$readn += $d;
				}elsif ($w eq 'D'){
					$refn += $d;
				}elsif ($w eq 'M'){
					if ($refn+1<=$loc and $loc<=$refn+$d){
						my $base = substr($sp[9] , $loc-$refn+$readn , 1);
						if ($type eq 'snp' and $base eq $alt){
							$altn++;
						}elsif ($type eq 'snp' and $base eq $ref){
							$refnn++;
						}
						last if $type eq 'snp';
					}
					$refn += $d;
					$readn += $d;
				}else{
					print STDERR "$d$w\n";
				}
			}
		}
		close BAM;
	}
	$chr =~ s/chr//;
	next if $refnn + $altn == 0;
	print "$chr\t$pos\t$ref\t$refnn\t$alt\t$altn\n";
}
close IN;

sub maplen{
	my ($m) = @_;
	my ($len) = (0);
	while ($m =~ /(\d+)([MIDSH])/g){
		if ($2 eq 'M'){
			$len += $1;
		}elsif ($2 eq 'D'){
			$len += $1;
		}
	}
	return $len;
}



sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: FilterAmpliconMutByBam.pl
#
#        USAGE: ./FilterAmpliconMutByBam.pl  
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
#      CREATED: 07/09/18 13:31:25
#     REVISION: ---
#===============================================================================
EOF!
}



