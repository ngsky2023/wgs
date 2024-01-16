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

my ($bam , $input , $help , $outdir);
my $virus = "/share/work3/wangrr/DB/virus/Hepatitis_B_virus/GCF_000861825.2_ViralProj15428_genomic.fna";
GetOptions(
	"i|input=s"	=>	\$input,
	"b|bam=s"	=>	\$bam,
	"v|virus=s"	=>	\$virus,
	"o|outdir=s"	=>	\$outdir,
	"help"	=>	\$help,
);

if ($help){
	&help;
	exit;
}
mkpath($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

my %fq;
if ($input){
	open IN , "$input";
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		if (exists $fq{$F[0]}){
			$fq{$F[0]}->[0] .= "," . $F[1];
			$fq{$F[0]}->[1] .= "," . $F[2];
		}else{
			$fq{$F[0]}->[0] = $F[1];
			$fq{$F[0]}->[1] = $F[2];
		}
	}
	close IN;
}
my %bam;
if ($bam){
	open IN , "$bam";
	while (<IN>){
		chomp;
		my ($sample , $aln) = split /\t/ , $_;
		$bam{$sample} = $aln;
	}
	close IN;
}

my $opt = "thread_no      = 30
detect_integration = yes
detect_mutation    = yes
blastn_bin      = /share/public/software/ncbi-blast-2.2.25+/bin/blastn
bowtie_bin      = /share/public/software/bowtie2-2.3.2/bowtie2
bwa_bin         = /share/public/software/bwa-0.7.12/bwa
trinity_script  = /share/work3/wangrr/local/trinityrnaseq_r20140717/Trinity
SVDetect_dir    = /share/public/software/SVDetect_r0.8b/
virus_database     = $virus
bowtie_index_human = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/HG19_BGI_Bowtie2/hg19
blastn_index_human = /share/work2/baijian/HumanRefs_Chaoyang/HumanRefs/HG19_BGI_Blast/hg19
blastn_index_virus = $virus
detection_mode     = sensitive
flank_region_size  = 10000
sensitivity_level  = 3
min_contig_length  = 100
blastn_evalue_thrd = 0.05
similarity_thrd    = 0.8 
chop_read_length   = 25
minIdentity        = 90";

open CMD , ">cmd.sh";
for my $sample (sort keys %fq){
	open O , ">$outdir/$sample.cfg";
	mkpath("$outdir/$sample");
	my ($fq1 , $fq2) = @{$fq{$sample}};
	print O << "EOF!";
fastq1        = $fq1
fastq2        = $fq2
$opt
EOF!
	close O;
	print CMD "/share/work3/wangrr/local/VirusFinder2.0/VirusFinder.pl -c $outdir/$sample.cfg --output $outdir/$sample\n";
}
for my $sample (sort keys %bam){
	open O , ">$outdir/$sample.cfg";
	mkpath("$outdir/$sample");
	my $bamf = $bam{$sample};
	print O << "EOF!";
alignment_file = $bamf
$opt
EOF!
	close O;
	print CMD "/share/work3/wangrr/local/VirusFinder2.0/VirusFinder.pl -c $outdir/$sample.cfg --output $outdir/$sample\n";
}
close CMD;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: VF.pl
#
#        USAGE: ./VF.pl  
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
#      CREATED: 11/21/17 09:16:17
#     REVISION: ---
#===============================================================================
EOF!
}



