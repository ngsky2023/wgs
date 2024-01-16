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

my ($input , $help , $dir);
my $affinity = "/share/public/software/PSSMHCpan-master/PSSMHCpan-1.0/PSSMHCpan-1.0.pl";
GetOptions(
	"i|input=s"	=>	\$input,
	"d|dir=s"	=>	\$dir,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}
mkdir $dir unless -d $dir;

my %fh;
open IN , "$input";
<IN>;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$F[11] =~ s/[\*:]//g;
	unless (exists $fh{$F[11]}){
		open $fh{$F[11]} , ">$dir/$F[11].fa";
	}
	print {$fh{$F[11]}} ">" , join("_",@F[0..5,10,16]) , "\n$F[15]\n";
}
close IN;

for my $fh (keys %fh){
	close $fh{$fh};
	print readpipe("$affinity $dir/$fh.fa 9 $fh | grep -w PSSMHCpan >> $dir/affinity.txt");
}

my @out;
my %e;
open IN , "$dir/affinity.txt";
while (<IN>){
	next unless /\t/;
	chomp;
	my @F = split /\t/ , $_;
	push @out , [@F];
}
close IN;

open O , ">$dir/affinity.result.txt";
for my $ll (sort {$a->[4] <=> $b->[4]} @out){
	my @s = split /_/ , $ll->[2];
	next if exists $e{"$s[0]\t$s[1]\t$s[2]"} or $s[2] eq 'NA';
	$e{"$s[0]\t$s[1]\t$s[2]"} = 1;
	$ll->[1] =~ s/(\d\d)(\d\d)$/*$1:$2/;
	print O "$s[0]\t$s[1]\t$s[2]\t$s[3]\t$s[4]\t$s[6]\t" , $ll->[1] , "\t" , $ll->[3],"\t$s[7]\t" , $ll->[4] , "\t" , $ll->[5] , "\n";
}
close O;




sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: pastePvacseq.pl
#
#        USAGE: ./pastePvacseq.pl  
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
#      CREATED: 01/23/18 12:57:58
#     REVISION: ---
#===============================================================================
EOF!
}



