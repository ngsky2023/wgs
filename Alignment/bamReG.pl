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

my ($input , $help , $check);
GetOptions(
	"i|input=s"	=>	\$input,
	"c|check=s"	=>	\$check,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

#A00262:57:H37FJDSXX:3:2304:19081:28588  1201    HBV_gtB 2054    10      85S65M  HBV_gtC 2382    0       CCGCAGCCACCCTGCCG
my %len;
open IN , "samtools view -H $input | ";
#@SQ     SN:HBV_gtD      LN:3182
while (<IN>){
	chomp;
	next unless /^\@SQ/;
	my ($id , $len) = (/SN:(\S+)\tLN:(\d+)/);
	$len{$id} = $len;
}
close IN;

my %v;
my @zone;
open IN , "$check";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	if ($#F > 1){
		push @zone , [@F[0..2]];
	}elsif ($#F == 0){
		push @zone , [$F[0] , 1 , $len{$F[0]}];
	}
}
close IN;

for my $ze (@zone){
	my @ze = @$ze;
	open IN , "samtools view -q 10 -F 1024 $input $ze[0]:$ze[1]-$ze[2] | grep -P 'chr\\w+' |";
	while (<IN>){
		chomp;
		my ($id,$flag,$chr,$loc,$q,$mode,$chro , @bb) = split /\t/ , $_;
		next if $mode =~ /\d+[SH].*\d+[SH]/;
		next if $mode !~ /\d+[SH]/ and $chro eq '=';
		next unless $chro =~ /chr/ or $chro eq '=';
		next if $chr !~ /chr/ and $chro !~ /chr/;
		$flag = sprintf("%b" , $flag);
		my $R;
		if ($flag =~ /1\d{6}$/){
			$R = 1;
		}else{
			$R = 2;
		}
		my $vcut = $loc;
		my $ccut;
		my $chrs;
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
							$chrs = $chrbb[0];
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
			$v{"$chr\t$vcut\t$chrs\t$ccut"}->{cut}++;
		}else{
			$chrs = $chro;
			$ccut = $bb[0];
			if ($R == 1 and $flag =~ /0\d{4}$/){
				$vcut += $maplen;
			}
			$v{"$chr\t$vcut\t$chrs\t$ccut"}->{pair}++;
		}
	}
	close IN;
}

for my $key (sort keys %v){
	print "$key";
	for my $type ('cut'){
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



