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

my $virusid = " HBV_gtA HBV_gtB HBV_gtC HBV_gtD HBV_gtE HBV_gtF HBV_gtG HBV_gtH";

#A00262:57:H37FJDSXX:3:2304:19081:28588  1201    HBV_gtB 2054    10      85S65M  HBV_gtC 2382    0       CCGCAGCCACCCTGCCG
my %v;
open IN , "samtools view -q 10 -F 1024 $input $virusid | grep -P 'chr\\w+' |";
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
							my $ml = &maplen($mode);
							$ccut += $ml-1;
						}
						$k = 1;
						last;
					}
				}
				last;
			}
		}
	}
	if ($k == 1){
		$v{"$hbv\t$vcut\t$chr\t$ccut"}->{cut}++;
	}else{
		$chr = $chro;
		$ccut = $bb[0];
		if ($R == 1 and $flag =~ /0\d{4}$/){
			$vcut += $maplen;
		}
		$v{"$hbv\t$vcut\t$chr\t$ccut"}->{pair}++;
	}
}
close IN;

for my $key (sort keys %v){
	print "$key";
	for my $type ('cut' , 'pair'){
		if (exists $v{$key}->{$type}){
			print "\t" , $v{$key}->{$type};
		}else{
			print "\t0";
		}
	}
	print "\n";
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

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: bamVirus.pl
#
#        USAGE: ./bamVirus.pl  
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
#      CREATED: 03/14/18 15:28:16
#     REVISION: ---
#===============================================================================
EOF!
}



