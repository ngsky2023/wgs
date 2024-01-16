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

my ($input , $help , $sampledir);
GetOptions(
	"i|input=s"	=>	\$input,
	"s|sampledir=s"	=>	\$sampledir,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %sample;
for my $file (<$sampledir/*.20K>){
	my $base = basename $file;
	$base =~ s/.20K//;
	$sample{$base} = '';
}

open IN , "$input";
while (<IN>){
	my @F = split /\t/, $_;
	$sample{$F[0]} = $_;
}
close IN;

for my $sample (sort keys %sample){
	if ($sample{$sample} eq ''){
		print "$sample" , ("\t0" x 8) , "\n";
	}else{
		print $sample{$sample};
	}
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: addZoneSample.pl
#
#        USAGE: ./addZoneSample.pl  
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
#      CREATED: 08/13/2018 11:04:01 AM
#     REVISION: ---
#===============================================================================
EOF!
}



