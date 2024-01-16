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
	-infile          <FILE>  input bam file (necessary)
	-outdir          <DIR>   output directory, usually the Sample Directory (default: ./ )
	-min             <NUM>   Mindepth(default:4)
	-max             <NUM>   Maxdepth(default:500)
	-ref             <FILE>  Reference file of genome (default:/home/zhanghk/software/CancerPipe/hg19/hg19.fa)
	-bed             <FILE>  Bed region file name, for exom and target region capture.
	-gatk            <FILE>  The all path of GATK (default: /home/leiyoubing/software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar)
	-samtools        <FILE>  The all path of samtools (default: /home/leiyoubing/software/samtools-0.1.19/samtools)
	-bcftools        <DIR>   The dir of bcftools (default: /share/software/software/samtools-0.1.18/bcftools)
	-ruby		 <FILE>	 ruby language (default: /home/sunfl/workdir/software/ruby-2.2.2/ruby)
	-ASNP            <FILE>  ATLAS2_SNP (default: /home/sunfl/workdir/software/Atlas2_v1.4.3/Atlas-SNP2/Atlas-SNP2.rb)
	-AINDEL          <FILE>  ATLAS2_SINDEL (default: /home/sunfl/workdir/software/Atlas2_v1.4.3/Atlas-Indel2/Atlas-Indel2.rb)
	-Platypus        <FILE>  PlatyPus_INDEL (default: /home/sunfl/workdir/software/Platypus_0.8.1/Platypus.py)

Example:

perl $0 -infile \$sample.bam -outdir \$outdir/\$sample \$call_variant -ref \$ref -bed \$bed -gatk \$gatk -samtools \$samtools -bcftools \$bcftools -ruby \$ruby -ASNP \$Atlas_SNP -AINDEL \$Atlas_Indel -Platypus  \$Platypus

";

my ($info,$software,$Outdir,$Mindepth,$Maxdepth,$ref,$bed,$GATK,$samtools,$bcftools,$ruby,$Atlas_SNP,$Atlas_Indel,$Platypus,$MustGiven);
GetOptions(
	"infile:s"=>\$info,
	"outdir:s"=>\$Outdir,
	"min:s"=>\$Mindepth,
	"max:s"=>\$Maxdepth,
	"ref:s"=>\$ref,
	"bed:s"=>\$bed,
	"gatk:s"=>\$GATK,
	"samtools:s"=>\$samtools,
	"bcftools:s"=>\$bcftools,
	"ruby:s"=>\$ruby,
	"ASNP:s"=>\$Atlas_SNP,
	"AINDEL:s"=>\$Atlas_Indel,
	"Platypus:s"=>\$Platypus,
	"MustGiven=s"=>\$MustGiven,
);die "$usage\n" if(!$info || !$bed);

$Outdir ||= '.';
$Outdir = abs_path($Outdir);
$software ||= 'GATK';
$Mindepth ||=4;
$Maxdepth ||=1000;
my $clinvar_P = "/share/public/database/Gynecological_cancer_backup/IPDB/ver4/DB_clinvar20180326_GRCh37_P.txt";
$ref ||= "/home/zhanghk/software/CancerPipe/hg19/hg19.fa";
$GATK ||= "/home/leiyoubing/software/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar";
$samtools ||= "/home/leiyoubing/software/samtools-0.1.19/samtools";
$bcftools ||= "/share/software/software/samtools-0.1.18/bcftools";
$ruby ||= "/home/sunfl/workdir/software/ruby-2.2.2/ruby";
$Atlas_SNP ||= "/home/sunfl/workdir/software/Atlas2_v1.4.3/Atlas-SNP2/Atlas-SNP2.rb";
$Atlas_Indel ||= "/home/sunfl/workdir/software/Atlas2_v1.4.3/Atlas-Indel2/Atlas-Indel2.rb";
$Platypus ||="/home/sunfl/workdir/software/Platypus_0.8.1/Platypus.py";

my $bam_pindel = "/share/work2/staff/sunchy/apps/pindel/bam2pindel.pl";
my $pindel2vcf = "/share/software/data/pindel_0.2.4p/pindel2vcf";
my $pindel = "$Bin/get_pindel.pl";

my $sample = basename $Outdir;
system("mkdir -m 755 -p $Outdir/3.Variants/snp" ) unless(-e "$Outdir/3.Variants/snp");
system("mkdir -m 755 -p $Outdir/3.Variants/indel") unless(-e "$Outdir/3.Variants/indel");
#system("mkdir -p $outdir/$tmp[0]/2.Realign") unless(-e "$outdir/$tmp[0]/2.Realign");

###SNP Calling results ######
###SNP//UnifiedGenotyper
open SNP,">$Outdir/3.Variants/snp/UnifiedGenotyper_$sample.sh";
print SNP "echo ========== GATK UnifiedGenotyper call snp start at : \`date\`==========\n";
print SNP "java -Xmx6G -jar $GATK -T UnifiedGenotyper  -R $ref -I $info -L $bed -o $Outdir/3.Variants/snp/$sample\_UnifiedGenotyper.snp.raw.vcf -stand_call_conf 50 -stand_emit_conf 10.0 -A RMSMappingQuality -A BaseCounts -baq CALCULATE_AS_NECESSARY -glm SNP\n";
print SNP "java -jar $GATK -T VariantFiltration -R $ref  -o $Outdir/3.Variants/snp/$sample\_UnifiedGenotyper.snp.raw.vcf.tmp -V $Outdir/3.Variants/snp/$sample\_UnifiedGenotyper.snp.raw.vcf --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP < $Mindepth || DP >$Maxdepth\" --filterName \"LOW_READ_DEPTH\"\n";
print SNP "head -100 $Outdir/3.Variants/snp/$sample\_UnifiedGenotyper.snp.raw.vcf.tmp|awk '\$0~/^#/' >$Outdir/3.Variants/snp/$sample.UnifiedGenotyper.snp.vcf\n";
print SNP "awk '\$7~/PASS/' $Outdir/3.Variants/snp/$sample\_UnifiedGenotyper.snp.raw.vcf.tmp >>$Outdir/3.Variants/snp/$sample.UnifiedGenotyper.snp.vcf\n";
print SNP "echo ========== UnifiedGenotyper call snp end at : \`date\`==========\n";
print SNP "echo \"GATK UnifiedGenotyper call snp finish\" >UnifiedGenotyper_$sample.mark\n";
close SNP;

###SNP//HaplotypeCaller
open SNP,">$Outdir/3.Variants/snp/HaplotypeCaller_$sample.sh";
print SNP "echo ========== GATK HaplotypeCaller call snp start at : \`date\`==========\n";
print SNP "java -Xmx6G -jar $GATK -T HaplotypeCaller -R $ref -I $info -L $bed -o $Outdir/3.Variants/snp/$sample\_HaplotypeCaller.snp.raw.vcf -stand_call_conf 50 -stand_emit_conf 10.0  -A RMSMappingQuality -A BaseCounts\n";
print SNP "java -jar $GATK -T VariantFiltration  -R $ref  -o $Outdir/3.Variants/snp/$sample\_HaplotypeCaller.snp.raw.vcf.tmp -V $Outdir/3.Variants/snp/$sample\_HaplotypeCaller.snp.raw.vcf --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"DP < $Mindepth || DP >$Maxdepth\" --filterName \"LOW_READ_DEPTH\"\n";
print SNP "head -100 $Outdir/3.Variants/snp/$sample\_HaplotypeCaller.snp.raw.vcf.tmp|awk '\$0~/^#/' >$Outdir/3.Variants/snp/HaplotypeCaller_$sample.vcf\n";
print SNP "awk '\$7~/PASS/' $Outdir/3.Variants/snp/$sample\_HaplotypeCaller.snp.raw.vcf.tmp >>$Outdir/3.Variants/snp/HaplotypeCaller_$sample.vcf\n";
print SNP "perl $Bin/split_snp_indel.pl $Outdir/3.Variants/snp/HaplotypeCaller_$sample.vcf $Outdir/3.Variants/snp/$sample.HaplotypeCaller\n";
print SNP "cp $Outdir/3.Variants/snp/$sample.HaplotypeCaller.indel.vcf $Outdir/3.Variants/indel/\n";
print SNP "echo ========== HaplotypeCaller call snp end at : \`date\`==========\n";
print SNP "echo \"$sample call snp finish\" >HaplotypeCaller_$sample.mark\n";
close SNP;

###SNP//atlas2
open SNP,">$Outdir/3.Variants/snp/atlas2_$sample.sh";
print SNP "echo ========== atlas call snp start at : \`date\`==========\n";
print SNP "$ruby $Atlas_SNP -i $info -r $ref -o $Outdir/3.Variants/snp/$sample.atlas.snp.raw -n $sample --Illumina -y $Mindepth -f $Maxdepth\n\n";
print SNP "perl $Bin/restrict2bed.pl $bed $Outdir/3.Variants/snp/$sample.atlas.snp.raw.vcf > $Outdir/3.Variants/snp/$sample.atlas.snp.vcf\n\n";
print SNP "echo ========== atlas call snp end at : \`date\`==========\n";
print SNP "echo \"atlas call snp finish\" >atlas_$sample.mark\n";
close SNP;

###INDEL Calling results ######
###InDel//Atlas2
open INDEL,">$Outdir/3.Variants/indel/indel_atlas_$sample.sh";
print INDEL "echo ========== atlas call indel start at : \`date\`==========\n";
print INDEL "$ruby $Atlas_Indel -b $info -r $ref -o $Outdir/3.Variants/indel/$sample\_atlas.indel.raw.vcf -I\n\n";
print INDEL "perl $Bin/restrict2bed.pl $bed $Outdir/3.Variants/indel/$sample\_atlas.indel.raw.vcf > $Outdir/3.Variants/indel/$sample.atlas.indel.vcf\n";
print INDEL "echo \"$sample atlas call indel finish\" >atlas_indel_$sample.mark\n";
print INDEL "echo ========== atlas call indel end at : \`date\`==========\n";
close INDEL;

###InDel//PlatyPus
open INDEL,">$Outdir/3.Variants/indel/indel_platypus_$sample.sh";
print INDEL "echo ========== platypus call indel start at : \`date\`==========\n";
print INDEL "python $Platypus callVariants --bamFiles=$info --output=$Outdir/3.Variants/indel/$sample\_PlatyPus.indel.raw.vcf --refFile=$ref --assemble=1 --nCPU=5\n\n";
print INDEL "perl $Bin/split_snp_indel.pl $Outdir/3.Variants/indel/$sample\_PlatyPus.indel.raw.vcf $Outdir/3.Variants/indel/$sample.PlatyPus\n";
print INDEL "perl $Bin/restrict2bed.pl $bed $Outdir/3.Variants/indel/$sample.PlatyPus.indel.vcf > $Outdir/3.Variants/indel/$sample.PlatyPus.indel.vcf1\n";
print INDEL "mv $Outdir/3.Variants/indel/$sample.PlatyPus.indel.vcf1 $Outdir/3.Variants/indel/$sample.PlatyPus.indel.vcf\n";
print INDEL "echo \"$sample platypus call indel finish\" >platypus_indel_$sample.mark\n";
print INDEL "echo ========== platypus call indel end at : \`date\`==========\n";
close INDEL;

###merging overlapped SNP and INDEL
open IN,">$Outdir/3.Variants/snp/overlap_snp.sh";
print IN "perl $Bin/overlap_var.pl $Outdir/3.Variants/snp/$sample.HaplotypeCaller.snp.vcf  $Outdir/3.Variants/snp/$sample.UnifiedGenotyper.snp.vcf  $Outdir/3.Variants/snp/$sample.atlas.snp.vcf >$Outdir/3.Variants/snp/$sample.snp.vcf\n";
print IN "echo \"$sample call snp finish\" >overlap_var_$sample.mark\n";
close IN;

open IN,">$Outdir/3.Variants/indel/overlap_indel.sh";
print IN "perl $Bin/overlap_var.pl $Outdir/3.Variants/indel/$sample.HaplotypeCaller.indel.vcf $Outdir/3.Variants/indel/$sample.atlas.indel.vcf $Outdir/3.Variants/indel/$sample.PlatyPus.indel.vcf >$Outdir/3.Variants/indel/$sample.indel.vcf\n";
print IN "perl $Bin/overlap_var1.pl $Outdir/3.Variants/indel/$sample.HaplotypeCaller.indel.vcf  $Outdir/3.Variants/indel/$sample.atlas.indel.vcf  $Outdir/3.Variants/indel/$sample.PlatyPus.indel.vcf > $Outdir/3.Variants/indel/$sample.indel1.vcf\n";
print IN "perl $Bin/SearchIndel.additional.pl  $Outdir/3.Variants/indel/$sample.indel1.vcf $clinvar_P $Outdir/3.Variants/indel/$sample.indel.vcf\n";
print IN "echo \"$sample call indel finish\" >overlap_var_$sample.mark\n";
close IN;

