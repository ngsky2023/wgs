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

my (@input , $help);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"help"	=>	\$help,
);

if ($help){
	&help;
	exit;
}

my %out;
my @sex;
my $i = 0;
my ($utotal , $num) = (0 , 0);
for my $file (@input){
	my %cn;
	open IN , "$file";
	<IN>;
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		if ($F[0] eq 'X' or $F[0] eq 'Y'){
			push @{$sex[$i]->{$F[0]}} , $F[-1];
		}elsif ($F[0] =~ /\d/){
			if (exists $cn{$F[0]}){
				$cn{$F[0]}++;
			}else{
				$cn{$F[0]} = 0;
			}
			push @{$out{$F[0]}->[$cn{$F[0]}]} , $F[-1];
			if ($F[-1] ne 'NA'){
				$utotal += $F[-1];
				$num++;
			}
		}
	}
	close IN;
	$i++;
}
my $umean = $utotal/$num;

for my $sex (@sex){
	my @X = @{$sex->{X}};
	my @Y = @{$sex->{Y}};
	my $x = mean(@X);
	#my $y = mean(@Y);
	if ($umean/$x > 1.5){
		for my $i (0..$#X){
			push @{$out{23}->[$i]} , $X[$i];
		}
		for my $i (0..$#Y){
			push @{$out{24}->[$i]} , $Y[$i];
		}
	}else{
		for my $i (0..$#X){
			push @{$out{25}->[$i]} , $X[$i];
		}
		for my $i (0..$#Y){
			push @{$out{26}->[$i]} , $Y[$i];
		}
	}
}

my $tn = 1250;
my $CV = '';
for my $chr (1..26){
	my $n;
	if (exists $out{$chr}){
		$CV .= ($chr+26) . cv($out{$chr});
		print  $chr , merge($out{$chr});
		$n = @{$out{$chr}};
	}
	my $tail = "\tNA" x ($tn-$n) . "\n";
	print $tail;
	$CV .= $tail;
}
print $CV;

sub cv{
	my @data = @{$_[0]};
	my $out = '';
	for my $data (@data){
		my @d;
		map{push @d , $_ if $_ ne 'NA'} @$data;
		my ($mean , $sd);
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(\@d);
		$mean = $stat->mean();
		$sd = $stat->standard_deviation();
		if ($mean and $mean =~ /\d/ and $mean != 0){
			$out .= "\t" . $sd/$mean;
		}else{
			$out .= "\tNA";
		}
	}
	return $out;
}

sub merge{
	my @data = @{$_[0]};
	my $out = '';
	for my $data (@data){
		if ($data->[0] eq 'NA'){
			$out .= "\tNA";
		}else{
#			my ($he , $i) = (0 , 0);
#			map{
#				if ($_ ne 'NA'){
#					$he += $_;
#					$i++;
#				}
#			} @$data;
#			$out .= "\t" . ($he/$i);
			my @dat;
			map{
				if ($_ ne 'NA'){
					push @dat , $_;
				}
			} @$data;
			@dat = sort {$a<=>$b} @dat;
			$out .= "\t" . $dat[int($#dat/2)];
		}
	}
	return $out;
}

sub mean{
	my @data = @_;
	my ($he , $i) = (0 , 0);
	for my $d (@data){
		if ($d ne 'NA'){
			$he += $d;
			$i++;
		}
	}
	return $he/$i;
}


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: makeRefCount.pl
#
#        USAGE: ./makeRefCount.pl  
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
#      CREATED: 06/07/18 17:35:45
#     REVISION: ---
#===============================================================================
EOF!
}



