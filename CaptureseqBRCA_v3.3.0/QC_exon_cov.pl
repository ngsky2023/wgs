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
        -s <s>        samtools;
        -h            Help information;
Author: lizhifeng; lizhifeng\@berrygenomics.com
Version: V1.1
Date: 2016.10.17
=====================================================================
USAGE
}

my ($in,$bed,$out,$samtools,$help);
GetOptions(
        "h|?|help"=>\$help,
        "i=s"=>\$in,
        "b=s"=>\$bed,
        "o=s"=>\$out,
        "s:s"=>\$samtools,
);
if(!defined($in)|| defined($help)){
        &Usage;
        exit 0;
}
$out ||="out";

$bed="/share/work2/staff/RD/Exon_seq/bed/AMPBRCA_target.bed" unless defined($bed);
$samtools="/share/software/software/samtools-1.3_install/bin/samtools" unless defined($samtools);
my $fout=$in.".result";

my ($from,$to,@temp,@tmp,$line,$Initial_bases_on_target,$Effective_sequences_on_target,$Base_covered_on_target,$Average_sequencing_depth_on_target,$Coverage_of_target_region,$min_depth);
open(IN,"$bed");
open(OUT,">$fout");
print OUT "Mininal_depth\tAverage_depth\tBases_Covered_on_target\tBases_on_target_bed\tCoverage\tInfo_of_bed\n";
while($line=<IN>)
{
chomp($line);
@temp=split(/\t/,$line,);
$from=$temp[1]+1;
$to=$temp[2];
`$samtools depth -m 1000000 -r $temp[0]:$from-$to $in > $out/tmp.txt`;
$Initial_bases_on_target=$temp[2]-$temp[1];
$Effective_sequences_on_target=`awk '{total+=\$3};\$3 != 0 {cov+=1};END{print total"\t"cov}' $out/tmp.txt`;
chomp($Effective_sequences_on_target);
($Effective_sequences_on_target,$Base_covered_on_target)=split(/\t/,$Effective_sequences_on_target);
if(! defined($Effective_sequences_on_target)){$Effective_sequences_on_target = 0;}
$Average_sequencing_depth_on_target=$Effective_sequences_on_target/$Initial_bases_on_target;
$Coverage_of_target_region=$Base_covered_on_target/$Initial_bases_on_target;
$min_depth=`sort -k3n $out/tmp.txt |head -n 1`;
chomp($min_depth);
@tmp=split(/\t/,$min_depth);
print OUT $tmp[2],"\t",$Average_sequencing_depth_on_target,"\t",$Base_covered_on_target,"\t",$Initial_bases_on_target,"\t",$Coverage_of_target_region,"\t",$line,"\n";
}



