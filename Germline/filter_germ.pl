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

#chr1    900009  .       CG      CGC,C   261     PASS    CIGAR=2M1I,1M1D;RU=C,G;REFREP=1,1;IDREP=2,0;MQ=53       GT:GQ:GQX:DPI:AD:ADF:A
#chr1    900285  .       C       T       964     PASS    SNVHPOL=2;MQ=60 GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL    1/1:229:21:77:1:0,77:0
#chr1    900286  .       A       G       975     PASS    SNVHPOL=2;MQ=60 GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL    1/1:229:21:77:0:0,77:0
#chr1    900298  .       C       G       420     PASS    SNVHPOL=6;MQ=60 GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL    0/1:222:19:80:2:30,50:
#chr1    900717  .       CTTAT   C       1843    PASS    CIGAR=1M4D;RU=TTAT;REFREP=2;IDREP=1;MQ=60       GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL
my $altn = 4;
my $vaf = 0.1;
if ($input =~ /\.gz/){
	open IN , "zcat $input|";
}else{
	open IN , "$input";
}
while (<IN>){
	unless (/^chr/){
		print $_;
		next;
	}
	next unless /\tPASS\t/;
	chomp;
	my @F = split /\t/ , $_;
	my @form = split /:/ , $F[-2];
	my @info = split /:/ , $F[-1];
	my %info;
	for my $i (0..$#form){
		$info{$form[$i]} = $info[$i];
	}
	my $ad = $info{AD};
	my @ad = split /,/ , $ad;
	my $mk = 0;
	my $depth = $ad[0];
	for my $ai (@ad[1..$#ad]){
		$mk = 1 if $ai >= $altn;
		$depth += $ai;
	}
	my $alt = (sort {$b<=>$a} @ad[1..$#ad])[0];
	my $af = $alt/$depth;
	#	my $vk = 0;
	#	for my $ai (@ad[1..$#ad]){
	#		$vk = 1 if $ai/$depth >= $vaf;
	#	}
	
	if ($mk and $af >= $vaf){
		print "$_\n";
	}else{
		print STDERR "$_\n";
	}
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: filter.pl
#
#        USAGE: ./filter.pl  
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
#      CREATED: 08/18/2019 09:50:26 AM
#     REVISION: ---
#===============================================================================
EOF!
}



