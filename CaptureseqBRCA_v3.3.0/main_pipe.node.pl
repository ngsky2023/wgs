#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use Cwd;

my $pwd = getcwd;
my $usage =
"
Usage:

  Options:
  -i|input      <FILE>   Input file record sample name and flowcell name.
  -o|outdir     <DIR>    Output file direcory, if this is not exists created, default(./).
  -c|config     <FILE>   Config file.
  -s|step       <num>    which step to run from,default is 0 [option: 0,1,2,3,4]
                          0: Preprocess (trim adapter, merge etc.)
		          1: Mapping and remove duplication reads (bwa mem)
                          2: Realignment
                          3: Variants calling
                          4: Annotation
                          5: Interpretation
			  6: Qulity control( running in step1 and step4)
  -cut          <num>    which step to cut-off end running,default is 4 [option: 1,2,3,4,5]
  -h|help                Help

For example:
        perl $0 -i sample.info -o ./analysis -c config.txt -s 0 -cut 5
";

my ($input,$outdir,$config,$step,$cut,$node,$help);
GetOptions(
  "i|input=s"=>\$input,
  "o|outdir=s"=>\$outdir,
  "c|config=s"=>\$config,
  "s|step=i"=>\$step,
  "cut=i"=>\$cut,
  "node=s"=>\$node,
  "h|help"=>\$help
);
$outdir ||= $pwd;
$step ||= 0;
$cut ||= 5;
$outdir = abs_path($outdir);

if (!$input || !$config || $help){
        die "$usage\n";
}
system("mkdir $outdir") unless (-e "$outdir");

###read config#######
my ($ref,$bwa,$samtools,$picard,$bwa_param_mem,$gatk,$annovar,$dbsnp,$bed,$bcftools,$bedtools,$call_variant,$software,$gene_info,$ruby,$ASNP,$AINDEL,$MustGiven,$Xindel,$hgvs,$Platypus,$flexbar,$rd,$adapter,$Interpretation,$random);
open CONFIG,"$config" || die "$!";
while(my $line = <CONFIG>){
	chomp($line);
	next if ($line =~ /^\s*$/ || $line =~ /^#/);
	if ($line =~ /ref\s*=\s*(.*)$/i){$ref = (split(/\;/,$1,2))[0];next;}
	if ($line =~ /bwa\s*=\s*(.*)$/i){$bwa = $1;next;}
	if ($line =~ /samtools\s*=\s*(.*)$/i){$samtools = $1;next;}
	if ($line =~ /picard\s*=\s*(.*)$/i){$picard = $1;next;}
        if ($line =~ /bwa_param_mem\s*=\s*(.*)$/i){$bwa_param_mem = $1;next;}
	if ($line =~ /gatk\s*=\s*(.*)$/i){$gatk = $1;next;}
	if ($line =~ /annovar\s*=\s*(.*)$/i){$annovar = $1;next;}
	if ($line =~ /bcftools\s*=\s*(.*)$/i){$bcftools = $1;next;}
	if ($line =~ /bedtools\s*=\s*(.*)$/i){$bedtools = $1;next;}
	if ($line =~ /dbsnp\s*=\s*(.*)$/i){$dbsnp = (split(/\;/,$1,2))[0];next;}
	if ($line =~ /bed\s*=\s*(.*)$/i){$bed = $1;next;}
	if ($line =~ /gene_info\s*=\s*(.*)$/i){$gene_info = $1;next;}
	if ($line =~ /call_variant\s*=\s*(.*)$/i){$call_variant = $1;next;}
	if ($line =~ /ruby\s*=\s*(.*)$/i){$ruby= $1;next}
	if ($line =~ /ASNP\s*=\s*(.*)$/i){$ASNP = $1;next;}
	if ($line =~ /AINDEL\s*=\s*(.*)$/i){$AINDEL = $1;next;}
	if ($line =~ /Platypus\s*=\s*(.*)$/i){$Platypus = $1;next;}
	if ($line =~ /flexbar\s*=\s*(.*)$/i){$flexbar =$1;next}
	if ($line =~ /rd\s*=\s*(.*)$/i){$rd=$1;next}
	if ($line =~ /adapter\s*=\s*(.*)$/i){$adapter=(split(/\;/,$1,2))[0];next}
	if ($line =~ /random\s*=\s*(.*)$/i){$random=$1;next}
	if ($line =~ /MustGiven\s*=\s*(.*)$/i){$MustGiven=(split(/\;/,$1,2))[0];next}
	if ($line =~ /Xindel\s*=\s*(.*)$/i){$Xindel=$1;next}
	if ($line =~ /hgvs\s*=\s*(.*)$/i){$hgvs=$1;next}
	if ($line =~ /Interpretation\s*=\s*(.*)$/i){$Interpretation=(split(/\;/,$1,2))[0];next}
        if ($line =~ /random\s*=\s*(.*)$/i){$random=$1;next}
	if ($line =~ /software\s*=\s*(.*)$/i){
		$software = $1;
		$software = $1 if ($software =~ /^(\w+)/);
		next;
	}
}
close CONFIG;
### Step0:Processing ### 
system("perl $Bin/Preprocess.pl -i $input -o $outdir -f $flexbar -r $rd -a '$adapter' -d $random");

###read sample_info#####
my @samples = ();
my %Info = ();
my $name=$input;
$name=~s/.*\///g;
open SAMPLE,"$outdir/Input_fq/sample.$name.txt" || die "$!";
while(my $line = <SAMPLE>){
	chomp($line);
	next if ($line=~ /^\s*$/ || $line=~ /^#/);
	$line =~ s/\s*$//;
	my @tmp = split /\s+/,$line;
	push(@samples,$tmp[0]);
	$Info{$tmp[0]}=[$tmp[1],$tmp[2],$tmp[3]];
}
close SAMPLE;
###make shell########


foreach my $sample(@samples){
	### StepI and II: mapping removing duplication reads and realignment ###
	system("perl $Bin/Map_Realign.pl -sample $sample  -read1  $Info{$sample}[1] -read2  $Info{$sample}[2] -outdir $outdir -ref $ref -bwa $bwa -samtools $samtools -interpretation $Interpretation  -picard $picard -bedtools $bedtools -bed $bed -BwaParam '$bwa_param_mem' -gatk $gatk  -dbsnp $dbsnp");
	### step III variants calling 
	if($Interpretation =~ /fast/){
		my $ref2="$ref.genome.fa";
		system("perl $Bin/Variants_calling.pl -infile  $outdir/$sample/2.Realign/$sample.bam  -outdir $outdir/$sample  $call_variant -ref $ref2 -bed $bed -gatk $gatk -samtools $samtools -bcftools $bcftools -MustGiven $MustGiven -ruby $ruby -ASNP $ASNP -AINDEL $AINDEL -Platypus  $Platypus");
	}else{
		system("perl $Bin/Variants_calling.pl -infile  $outdir/$sample/2.Realign/$sample.bam  -outdir $outdir/$sample  $call_variant -ref $ref -bed $bed -gatk $gatk -samtools $samtools -bcftools $bcftools -MustGiven $MustGiven  -ruby $ruby -ASNP	$ASNP -AINDEL $AINDEL -Platypus  $Platypus");
	}
	### step IV annotation ###
	system("perl $Bin/Annotation.pl -i $outdir/$sample/3.Variants/snp/$sample.snp.vcf -g $gene_info -interpretation $Interpretation  -bed $bed -s $sample -t snp -w $software -anno $annovar -o $outdir -v $hgvs");
	system("perl $Bin/Annotation.pl -i $outdir/$sample/3.Variants/indel/$sample.indel.vcf -g $gene_info -interpretation $Interpretation -bed $bed -s $sample -t indel -w $software -anno $annovar -o $outdir -v $hgvs");
	### step V ###
	 system("perl $Bin/Interpretation.pl -i $outdir/$sample/4.Annot -s $sample  -o $outdir -t $Interpretation -d $Xindel -k $MustGiven");
	### step VI QC ##
	system("perl $Bin/QC.pl -sample $sample -outdir $outdir -samtools $samtools -bedtools $bedtools -bed $bed -read  $Info{$sample}[1]"); 
}
###  qsub shell
if (!$node){
	foreach my $sample(@samples){
		system("nohup perl $Bin/qsub.pl -inputdir $outdir -sample $sample -step $step -type $Interpretation -cut $cut 1>$sample.sh.o 2>$sample.sh.e &");
	}
}
if ($node =~/local/){
	foreach my $sample(@samples){
		system("$Bin/qssh.pl -inputdir $outdir -sample $sample -step $step -type $Interpretation -cut $cut 1>$sample.sh.o 2>$sample.sh.e \n");
	}
}
else{
	foreach my $sample(@samples){
		system("nohup perl $Bin/tqsub.pl -node $node -inputdir $outdir -sample $sample -step $step -type $Interpretation -cut $cut 1>$sample.sh.o 2>$sample.sh.e &");
	}
}
