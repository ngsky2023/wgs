#!/usr/bin/perl -w
use strict;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use lib $Bin;

my $usage =
"
Usage:
   
Options:
	-sample	sample name
	-outdir	Output file direcory, if this is not exists created, default(./).'
	-read Input file record sample name and it's fastq file.
	-samtools	The all path of samtools,default is '/home/leiyoubing/software/samtools-0.1.19/samtools'
	-bedtools	The dir path of bedtools, default is '/home/leiyoubing/software/bedtools2-2.19.1/bin'
	-bed	Bed region file name, for exom and target region capture.
	-help                Help

For example:
perl $0 -sample \$sample -outdir \$outdir -samtools \$samtools -bedtools \$bedtools -bed \$bed -read \$read
 
";


my ($sample,$outdir,$samtools,$bedtools,$bed,$read,$help);
GetOptions(
  "sample=s"=>\$sample,
  "read=s"=>\$read,
  "outdir=s"=>\$outdir,
  "samtools=s"=>\$samtools,
  "bedtools=s"=>\$bedtools,
  "bed=s"=>\$bed,
  "help=s"=>\$help
);

$outdir ||= abs_path($outdir);
$samtools ||= '/home/leiyoubing/software/samtools-0.1.19/samtools';
$bedtools ||= '/home/leiyoubing/software/bedtools2-2.19.1/bin';
if (!$sample || !$bed || $help){
        die "$usage\n";
}

system("mkdir -m 755 -p $outdir/$sample/6.QC") unless (-e "$outdir/$sample/6.QC");
###QC
open  SH, ">$outdir/$sample/6.QC/QC.Mapping.$sample.sh" || die "$!";
print SH "#!/bin/bash\n";
print SH "echo ===start at : `date` ===\n";
print SH "perl $Bin/QC_exome_mem.pl -i $outdir/$sample/1.Map/$sample.sort.rmdup.bam -r $bed -c $read -o $outdir/$sample/6.QC -s $samtools -b $bedtools -t PE -plot\n";
#print SH "perl $Bin/QC_exon_cov.pl -i $outdir/$sample/1.Map/$sample.sort.rmdup.bam -b $bed -o $outdir/$sample/6.QC -s $samtools\n";

print SH "echo ===end at : `date` ===\n";
close SH;

open SH,">$outdir/$sample/6.QC/QC.anno.snp.$sample.sh" || die "$!";
print SH "#!/bin/bash\n";
print SH "echo ===snp annotate start at : `date` ===\n";
print SH "perl $Bin/annovar_statistic.Tp.pl -a $outdir/$sample/4.Annot/$sample\_snp.annovar.hg19_multianno.txt -s $outdir/$sample/6.QC/$sample -t snp\n\n";
print SH "echo \"$sample snp annotate finish\" >$outdir/$sample/6.QC/annotate_snp.$sample.mark\n";
print SH "echo ===snp annotate end at : `date` ===\n";
close SH;

open  SH, ">$outdir/$sample/6.QC/QC.anno.indel.$sample.sh" || die "$!";
print SH "#!/bin/bash\n";
print SH "echo ===indel annotate start at : `date` ===\n";
print SH "perl $Bin/annovar_statistic.Tp.pl -a $outdir/$sample/4.Annot/$sample\_indel.annovar.hg19_multianno.txt -s $outdir/$sample/6.QC/$sample -t indel\n\n";
print SH "echo \"$sample indel annotate finish\" >$outdir/$sample/6.QC/annotate_indel.$sample.mark\n";
print SH "echo ===indel annotate end at : `date` ===\n";
close SH;
