##!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

sub Usage{
        print STDERR <<USAGE;
=====================================================================
Description: QC target region of BRCA1/2;
Options:
        -i <s>        input bam file;
        -b <S>        bed file;
        -o <s>        Output dir name;
        -h            Help information;
Author: lizhifeng; lizhifeng\@berrygenomics.com
Version: V1.1
Date: 2016.10.17
=====================================================================
USAGE
}

my ($in,$bed,$out,$help);
GetOptions(
        "h|?|help"=>\$help,
        "i=s"=>\$in,
        "b=s"=>\$bed,
        "o=s"=>\$out,
);
if(!defined($in)|| defined($help)){
        &Usage;
        exit 0;
}
$out ||="out";

$bed="/share/work2/staff/RD/Exon_seq/bed/AMPBRCA_target.bed" unless defined($bed);
my $fout=$in.".result";

my ($n,$sum,$line,@temp);
open(IN,"samtools view $in|");
#open(OUT,">$fout");

while($line=<IN>)
{
chomp($line);
@temp=split(/\t/,$line,);
if((($temp[2] eq "chr13" && $temp[3]>32889616 && $temp[3]<32973809 ) || ($temp[2] eq "chr17" && $temp[3]>41196311 && $temp[3]<41277500)) && abs($temp[8])<500){$sum+=abs($temp[8]);$n+=1;}
}
close(IN);
#print "n $n sum $sum average ",$sum/$n,"\n";
print $sum/$n,"\n";



