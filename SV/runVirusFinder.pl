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

my ($r1 , $r2 , $outdir , $help);
my ($cpu) = (6);
GetOptions(
	"r1=s"	=>	\$r1,
	"r2=s"	=>	\$r2,
	"outdir=s"	=>	\$outdir,
	"cpu=s"	=>	\$cpu,
	"help"	=>	\$help,
);

if ($help or ! $r1 or ! $r2){
	&help;
	exit;
}

mkpath($outdir) unless -e $outdir;
$outdir = abs_path($outdir);


open CFG , ">$outdir/config.txt";
print CFG << "EOF!";
fastq1        = $r1
fastq2        = $r2

thread_no      = $cpu

detect_integration = yes
detect_mutation    = yes

blastn_bin      = /share/public/software/ncbi-blast-2.2.25+/bin/blastn
bowtie_bin      = /share/public/software/bowtie2-2.3.2/bowtie2
bwa_bin         = /share/public/software/bwa-0.7.12/bwa
trinity_script  = /share/work3/wangrr/local/trinityrnaseq_r20140717/Trinity
SVDetect_dir    = /share/public/software/SVDetect_r0.8b/bin/

virus_database     = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/Virus/viral_hs.fa
bowtie_index_human = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/HG19_BGI_Bowtie2/hg19
blastn_index_human = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/HG19_BGI_Blast/hg19
blastn_index_virus = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/Virus/viral_hs

detection_mode     = sensitive
flank_region_size  = 10000
sensitivity_level  = 3

min_contig_length  = 100
blastn_evalue_thrd = 0.05
similarity_thrd    = 0.8 
chop_read_length   = 25
minIdentity        = 90
EOF!
close CFG;

print "/share/work3/wangrr/local/VirusFinder2.0/VirusFinder.pl -c $outdir/config.txt --output $outdir";
	

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: runVirusFinder.pl
#
#        USAGE: ./runVirusFinder.pl  
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
#      CREATED: 09/28/17 16:50:38
#     REVISION: ---
#===============================================================================
EOF!
}

my $cfg = "
fastq1        = r1
fastq2        = r2

thread_no      = 6

detect_integration = yes
detect_mutation    = yes

blastn_bin      = /share/public/software/ncbi-blast-2.2.25+/bin/blastn
bowtie_bin      = /share/public/software/bowtie2-2.3.2/bowtie2
bwa_bin         = /share/public/software/bwa-0.7.12/bwa
trinity_script  = /share/work3/wangrr/local/trinityrnaseq_r20140717/Trinity
SVDetect_dir    = /share/public/software/SVDetect_r0.8b/bin/

virus_database     = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/Virus/viral_hs.fa
bowtie_index_human = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/HG19_BGI_Bowtie2/hg19
blastn_index_human = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/HG19_BGI_Blast/hg19
blastn_index_virus = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/Virus/viral_hs

detection_mode     = sensitive
flank_region_size  = 10000
sensitivity_level  = 3

min_contig_length  = 100
blastn_evalue_thrd = 0.05
similarity_thrd    = 0.8 
chop_read_length   = 25
minIdentity        = 90
";
