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

my ($input , $help , $prx);
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
	"prx=s"	=>	\$prx,
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

open IN , "$input";
my $txt = <IN>;
close IN;
chomp $txt;
my ($sample , $r1 , $r2 , $lane) = split /\t/ , $txt;
my @r1 = </$r1>;
my @r2 = </$r2>;

my $n = @r1;

open SHELL , ">$prx.cmd.sh";
my $sb = 'aa';
my $number = 0;
my @bam;
for my $i (0..$#r1){
	my $bx = '';
	open BX , "zcat $r1[$i]|";
	my $name = <BX>;
	close BX;
	if ($name =~ /BX:/){
		$bx = '-C';
	}
	print SHELL "$soft{bwa} mem $bx -t 4 -M -R '\@RG\\tID:$sample\\tPL:illumina\\tSM:$sample' $soft{genome} $r1[$i] $r2[$i] | $soft{samtools} view -bS - | $soft{samtools} sort - -m 4G -T $prx.$sb.tmp -o $prx.$sb.sort.bam && rm $r1[$i] $r2[$i] \n";
	push @bam , "$prx.$sb.sort.bam";
	$sb++;
	$number++;
}
close SHELL;

my $dir = dirname($prx);
my $tc = 20;
$tc = $n if $n < $tc;

print readpipe("cd $dir && qsub -cwd -t 1-$n -tc $tc -l vf=10g,p=2 -sync y -N bwa_$sample\_$lane $Bin/runtask.sh $prx.cmd.sh");

for my $file (@bam){
	while (! -e $file){
		sleep 60;
	}
}





sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: forMap.pl
#
#        USAGE: ./forMap.pl  
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
#      CREATED: 10/18/17 16:51:45
#     REVISION: ---
#===============================================================================
EOF!
}



