#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';
use FindBin '$Bin';
use lib $Bin;

my $usage =
"
Options:
	-infile         input bam file (necessary)
	-outdir         output directory, usually the Sample Directory (default: ./ )
	-ref            Reference file of genome (default:/home/zhanghk/software/CancerPipe/hg19/hg19.fa)
	-gatk           The all path of GATK (default: /home/leiyoubing/software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar)
	-samtools       The all path of samtools (default: /home/leiyoubing/software/samtools-0.1.19/samtools)
	-dbsnp          The file of dbsnp (default: /share/work1/staff/guofang/human/dbsnp141.vcf)

Example:
	perl $0 -infile bam -outdir sample_dir

";
my ($info,$Outdir,$ref,$GATK,$snp_vcf,$samtools) = ('','','','','','');
GetOptions(
        "infile:s"=>\$info,
        "outdir:s"=>\$Outdir,
	"ref:s"=>\$ref,
	"gatk:s"=>\$GATK,
	"samtools:s"=>\$samtools,
	"dbsnp:s"=>\$snp_vcf,
);
die "$usage\n" if(!$info);

$Outdir ||= '.';
$Outdir = abs_path($Outdir);
$ref ||= "/home/zhanghk/software/CancerPipe/hg19/hg19.fa";
$GATK ||= "/home/leiyoubing/software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar";
$snp_vcf ||= "/share/work1/staff/guofang/human/dbsnp141.vcf";
$samtools ||= "/home/leiyoubing/software/samtools-0.1.19/samtools";

my $sample_name = basename $Outdir;
mkdir ("$Outdir/variation",0755) unless(-d "$Outdir/variation");
mkdir ("$Outdir/variation/realn", 0755) unless(-d "$Outdir/variation/realn");

open RECALL,">$Outdir/variation/realn/GATK_recall_$sample_name.sh";
print RECALL "echo ========== GATK Realigner and BaseRecall start at : \`date\`==========\njava -Xmx3G -jar $GATK -T RealignerTargetCreator -R $ref -I $info -o $Outdir/variation/realn/$sample_name.intervals -known $snp_vcf\n";
print RECALL "java -jar $GATK -T IndelRealigner -R $ref -targetIntervals $Outdir/variation/realn/$sample_name.intervals -I $info -o $Outdir/variation/realn/$sample_name.realn.bam -known $snp_vcf\n";
print RECALL "java -jar $GATK -T BaseRecalibrator -R $ref -I $Outdir/variation/realn/$sample_name.realn.bam -knownSites $snp_vcf -o $Outdir/variation/realn/$sample_name.grp\n";
print RECALL "java -jar $GATK -T PrintReads -R $ref -I $Outdir/variation/realn/$sample_name.realn.bam -BQSR $Outdir/variation/realn/$sample_name.grp -o $Outdir/variation/realn/$sample_name.bam\n";
print RECALL "$samtools index $Outdir/variation/realn/$sample_name.bam\n";
print RECALL "rm -rf $Outdir/variation/realn/$sample_name.realn.bam\necho GATK Realigner and BaseRecall ==========end at : \`date\` ==========\n";
print RECALL "echo \"$sample_name recall finish\" >recall_$sample_name.mark\n";
close RECALL;

