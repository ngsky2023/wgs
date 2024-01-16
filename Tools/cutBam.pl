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

my ($input , $help , $outdir);
GetOptions(
	"i|input=s"	=>	\$input,
	"o|outdir=s"	=>	\$outdir,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %soft;
open IN , "$Bin/soft.list";
while (<IN>){
	chomp;
	my ($k , $d) = split /\t/ , $_;
	$soft{$k} = $d;
}
close IN;
my $samtools = $soft{samtools};

#@SQ     SN:chrM
my $base = basename $input;
$base =~ s/\.bam//;
my @chr;
my $un = '';
open IN , "$samtools view -H $input|";
while (<IN>){
	if (/\@SQ\tSN:(chr\w+)/){
		push @chr , $1;
	}elsif (/\@SQ\tSN:(\S+)/){
		$un .= " $1";
	}
}
close IN;

for my $chr (@chr){
	mkdir "$outdir/$chr" unless -d "$outdir/$chr";
	readpipe("$samtools view -b $input $chr > $outdir/$chr/$base.$chr.bam");
}
readpipe("$samtools view -b $input $un > $outdir/chrUn/$base.chrUn.bam");


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: cutBam.pl
#
#        USAGE: ./cutBam.pl  
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
#      CREATED: 05/18/18 09:58:57
#     REVISION: ---
#===============================================================================
EOF!
}



