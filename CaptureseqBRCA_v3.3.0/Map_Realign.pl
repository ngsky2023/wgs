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
	-sample  sample name	
	-read1   read1
	-read2   read2
	-input	Input file record sample name and it's fastq file.
	-outdir	Output file direcory, if this is not exists created, default(./).
	-ref	Reference for BWA, default is '/home/zhanghk/software/CancerPipe/hg19/hg19.fa'.
	-bwa	The all path of BWA,default is '/home/zhanghk/software/bwa-0.7.10/bwa'
	-BwaParam	The parameter of BWA,default is '-t 4 -M'
	-samtools	The all path of samtools,default is '/home/leiyoubing/software/samtools-0.1.19/samtools'
	-picard	The all path of picard,default is '/home/leiyoubing/software/picard-tools-1.57'
	-bedtools	The dir path of bedtools, default is '/home/leiyoubing/software/bedtools2-2.19.1/bin'
	-bed	Bed region file name, for exom and target region capture.
	-infile         input bam file (necessary)
	-gatk           The all path of GATK (default: /home/leiyoubing/software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar)
	-dbsnp          The file of dbsnp (default: /share/work1/staff/guofang/human/dbsnp141.vcf)
	-primerbed      primer bed
	-help                Help

For example:
perl $0 -sample \$sample -read1 \$read1 -read2 \$read2 -outdir \$outdir -ref \$ref -bwa \$bwa -samtools \$samtools -picard \$picard -bedtools \$bedtools -bed \$bed -BwaParam \$bwa_param_mem -gatk \$gatk  -dbsnp \$dbsnp
 
";


my ($bam,$sample,$read1,$read2,$outdir,$ref,$ref1,$ref2,$bwa,$BwaParam,$samtools,$picard,$interpretation,$bedtools,$bed,$GATK,$snp_vcf,$primerbed,$help);
GetOptions(
  "sample=s"=>\$sample,
  "read1=s"=>\$read1,
  "read2=s"=>\$read2,
  "outdir=s"=>\$outdir,
  "ref=s"=>\$ref,
  "bwa=s"=>\$bwa,
  "BwaParam=s"=>\$BwaParam,
  "samtools=s"=>\$samtools,
  "picard=s"=>\$picard,
  "interpretation=s"=>\$interpretation,
  "bedtools=s"=>\$bedtools,
  "bed=s"=>\$bed,
  "gatk=s"=>\$GATK,
  "dbsnp=s"=>\$snp_vcf,
  "primerbed=s"=>\$primerbed,
  "help=s"=>\$help
);
$outdir ||= abs_path($outdir);
$ref ||= '/home/zhanghk/software/CancerPipe/hg19/hg19.fa';
$bwa ||= '/home/zhanghk/software/bwa-0.7.10/bwa';
$BwaParam ||= '-t 4 -M';
$samtools ||= '/home/leiyoubing/software/samtools-0.1.19/samtools';
$picard ||= '/home/leiyoubing/software/picard-tools-1.57';
$bedtools ||= '/home/leiyoubing/software/bedtools2-2.19.1/bin';
$GATK ||= "/home/leiyoubing/software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar";
$snp_vcf ||= "/share/work1/staff/guofang/human/dbsnp141.vcf";
if (!$sample || !$bed || $help){
        die "$usage\n";
}
if($interpretation =~ /fast/){
	$ref1="$ref.panel.fa";
	$ref2="$ref.genome.fa";
}
system("mkdir -m 755 $outdir/$sample") unless (-d "$outdir/$sample");
system("mkdir -m 755 -p $outdir/$sample/1.Map") unless (-e "$outdir/$sample/1.Map");
&mapping($sample,$read1,$read2);
$bam="$outdir/$sample/1.Map/$sample.sort.rmdup.bam";
system("mkdir -m 755 -p $outdir/$sample/2.Realign") unless(-e "$outdir/$sample/2.Realign");
&realign($sample,$bam);


##########mapping and remove duplicated reads##################
sub mapping{
	my ($sample,$r1,$r2) = @_;
	open DUP,">$outdir/$sample/1.Map/bam_rmdup_$sample.sh";
	print DUP "#!/bin/bash\n";
	print DUP "echo ===start at : `date` ===\n";
	if($interpretation=~/fast/){
		print DUP "$bwa mem $BwaParam -R \"\@RG\\tID:$sample\\tPL:illumina\\tSM:$sample\\tCN:BERRY\" $ref1 $r1 $r2 | perl -ne 'my \@ln=split(/\\t/,\$_);if (/^\\\@/) {print \$_;}elsif (!(\$ln[1] & 0x4)) {print \$_;}' > $sample.sam\n";
		print DUP "perl $Bin/change_coordination.pl $sample $sample.sam  /home/sunfl/workdir/database/cache.seq > $sample.sam1\n";
		print DUP "mv $sample.sam1 $sample.sam\n";
		print DUP "$samtools view -bS -q 60 $sample.sam > $sample.bam\n";
	}else{
		#print DUP "$bwa mem $BwaParam -R \"\@RG\\tID:$sample\\tPL:illumina\\tSM:$sample\\tCN:BERRY\" $ref $r1 $r2 | perl -ne 'my \@ln=split(/\\t/,\$_);if (/^\\\@/) {print \$_;}elsif (!(\$ln[1] & 0x4)) {print \$_;}' > $sample.sam\n";
		if ($primerbed){
		print DUP "$bwa mem $BwaParam -R \"\@RG\\tID:$sample\\tPL:illumina\\tSM:$sample\\tCN:BERRY\" $ref $r1 $r2 | perl $Bin/maskPrimerBase.pl -b $primerbed -i - | perl -ne 'my \@ln=split(/\\t/,\$_);if (/^\\\@/) {print \$_;}elsif (!(\$ln[1] & 0x4)) {print \$_;}' > $sample.sam\n";
		}
		else{
		print DUP "$bwa mem $BwaParam -R \"\@RG\\tID:$sample\\tPL:illumina\\tSM:$sample\\tCN:BERRY\" $ref $r1 $r2 | perl -ne 'my \@ln=split(/\\t/,\$_);if (/^\\\@/) {print \$_;}elsif (!(\$ln[1] & 0x4)) {print \$_;}' > $sample.sam\n";
		}
		#print DUP "$bwa mem $BwaParam -R \"\@RG\\tID:$sample\\tPL:illumina\\tSM:$sample\\tCN=BERRY\" $ref $r1 $r2 | perl $Bin/filterMem.pl -i - | perl $Bin/maskPrimerBase.pl -i - | perl -ne 'my \@ln=split(/\\t/,\$_);if (/^\\\@/) {print \$_;}elsif (!(\$ln[1] & 0x4)) {print \$_;}' > $sample.sam\n";
		print DUP "$samtools view -bS -q 20 $sample.sam > $sample.bam\n";
	}
	print DUP "rm -rf $sample.sam\n";
	print DUP "$samtools sort -o $sample.sort.bam $sample.bam\n";	
	print DUP "$samtools flagstat $sample.sort.bam 1>$sample.sort.flagstat\n";
        print DUP "$samtools index $sample.sort.bam\n";
	print DUP "java -Xmx6g -jar $picard/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT INPUT=$sample.sort.bam OUTPUT=$sample.sort.rmdup.bam METRICS_FILE=$sample.rmdup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000\n";
	print DUP "$samtools index $sample.sort.rmdup.bam\n";
        #print DUP "$samtools index $sample.sort.Header.bam $sample.sort.Header.bai\n";
	print DUP "echo \"$sample bwa & rmdup finish\" >rmdup_$sample.mark\n";
	print DUP "echo ===end at : `date` ===\n";
	close DUP;
}


########## realignment ##################
sub realign{
	my ($sample,$info)=@_;
	open RECALL,">$outdir/$sample/2.Realign/GATK_recall_$sample.sh";
	print RECALL "echo ========== GATK Realigner and BaseRecall start at : \`date\`==========\n";
	if($interpretation=~/fast/){
		print RECALL "java -Xmx5G -jar $GATK -T RealignerTargetCreator -nt 3 -R $ref2 -I $info -o $outdir/$sample/2.Realign/$sample.intervals -known $snp_vcf\n";
		print RECALL "java -jar $GATK -T IndelRealigner  -R $ref2 -targetIntervals $outdir/$sample/2.Realign/$sample.intervals -I $info -o $outdir/$sample/2.Realign/$sample.realn.bam -known $snp_vcf\n";
		print RECALL "java -jar $GATK -T BaseRecalibrator  -R $ref2 -I $outdir/$sample/2.Realign/$sample.realn.bam -knownSites $snp_vcf -o $outdir/$sample/2.Realign/$sample.grp\n";
		print RECALL "java -jar $GATK -T PrintReads  -R $ref2 -I $outdir/$sample/2.Realign/$sample.realn.bam -BQSR $outdir/$sample/2.Realign/$sample.grp -o $outdir/$sample/2.Realign/$sample.bam\n";
	}else{
		print RECALL "java -Xmx5G -jar $GATK -T RealignerTargetCreator -nt 3 -R $ref -I $info -o $outdir/$sample/2.Realign/$sample.intervals -known $snp_vcf\n";
                print RECALL "java -jar $GATK -T IndelRealigner  -R $ref -targetIntervals $outdir/$sample/2.Realign/$sample.intervals -I $info -o $outdir/$sample/2.Realign/$sample.realn.bam -known $snp_vcf\n";
                print RECALL "java -jar $GATK -T BaseRecalibrator  -R $ref -I $outdir/$sample/2.Realign/$sample.realn.bam -knownSites $snp_vcf -o $outdir/$sample/2.Realign/$sample.grp\n";
                print RECALL "java -jar $GATK -T PrintReads  -R $ref -I $outdir/$sample/2.Realign/$sample.realn.bam -BQSR $outdir/$sample/2.Realign/$sample.grp -o $outdir/$sample/2.Realign/$sample.bam\n";
	}
	print RECALL "$samtools index $outdir/$sample/2.Realign/$sample.bam $outdir/$sample/2.Realign/$sample.bai\n";
	print RECALL "rm -rf $outdir/$sample/2.Realign/$sample.realn.bam\n";
	print RECALL "rm -rf $outdir/$sample/2.Realign/$sample.realn.bai\n";
	print RECALL "echo GATK Realigner and BaseRecall ==========end at : \`date\` ==========\n";
	print RECALL "echo \"$sample recall finish\" >recall_$sample.mark\n";
	close RECALL;
}
