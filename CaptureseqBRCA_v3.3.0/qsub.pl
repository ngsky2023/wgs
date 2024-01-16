#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
my $usage =
"
Usage:

  Options:
    -inputdir      <dir>         The path of project directory
    -sample        <str>         Sample name
    -help                        Help message
    -step          <num>         which step to run from,default is 1 [option: 1,2,3,4,5,6]
                                  1: bwa aln
                                  2: bwa sampe, rmdup and QC 
                                  3: call snp, call indel
                                  4: annotation
                              
    -cut           <num>         which step to cut-off end running,default is 6 [option: 1,2,3,4,5,6]
  For example:
    perl $0 -inputdir ~/work/DNA/cancer_test2/analysis -sample BJ0717 -step 1 -cut 6
";
my ($shell_dir,$sample,$step,$cut,$help,$type,$finish);
GetOptions(
  "inputdir=s"=>\$shell_dir,
  "sample=s"=>\$sample,
  "step=i"=>\$step,
  "cut=i"=>\$cut,
  "type=s"=>\$type,
  "h|help"=>\$help
);
$step ||= 0;
$cut ||= 4;
my $cmd="";
my $qsub="qsub -cwd ";
#$qsub = "$qsub -q $que " if ($que);
if (!$shell_dir || !$sample || $help){
	die "$usage\n";
}
die "start-step must bigger than the end-step\n" if ($step > $cut);
if ($step <=0 && $cut>=0){
        $cmd="cd $shell_dir/Input_fq/1.Trim/ && $qsub -l p=5 -l vf=4g Preprocess_$sample.sh";
        print "$cmd\n";
	system("cd $shell_dir/Input_fq/1.Trim/ && $qsub -l p=5 -l vf=4g Preprocess_$sample.sh");
	$finish = 1;	
	while($finish){
		sleep 60;
		if(-e "$shell_dir/Input_fq/1.Trim/Preprocess_$sample.mark"){
			$finish=0;
		}
	}
}

##qsub rmdup and QC####

if ($step <=1 && $cut>=1){
	foreach my $rmdup(glob "$shell_dir/$sample/1.Map/bam_rmdup_*.sh"){
		my $dirname = dirname $rmdup;
                $cmd="cd $dirname && $qsub -l p=5 -l vf=8g $rmdup";
                print "$cmd\n";
		system("cd $dirname && $qsub -l p=5 -l vf=8g $rmdup");
	}
	$finish = 1;
	while($finish){
		sleep 60;
		if (-e "$shell_dir/$sample/1.Map/rmdup_$sample.mark"){
			die "$sample.sort.rmdup.bam not exists!\n" unless (-e "$shell_dir/$sample/1.Map/$sample.sort.rmdup.bam");	
			$finish = 0;
		}
	}
        $cmd="cd $shell_dir/$sample/6.QC && $qsub -l vf=3g QC.Mapping.$sample.sh";
        print "$cmd\n";
	system("cd $shell_dir/$sample/6.QC && $qsub -l vf=3g QC.Mapping.$sample.sh");
}

##qsub Realignment ####
if ($step <=2 && $cut>=2){
	if (-e "$shell_dir/$sample/2.Realign"){
                 $cmd="cd $shell_dir/$sample/2.Realign && $qsub -l p=1 -l vf=12g  GATK_recall_$sample.sh";
                 print "$cmd\n";
		system("cd $shell_dir/$sample/2.Realign && $qsub -l p=1 -l vf=12g  GATK_recall_$sample.sh");
		$finish = 1;
		while($finish){
			sleep 60;
			if (-e "$shell_dir/$sample/2.Realign/recall_$sample.mark"){
				$finish = 0;
				die "$sample.bam not exist!\n" unless(-e "$shell_dir/$sample/2.Realign/$sample.bam");
			}else{
				$finish = 1;
			}
		}
	}
}
##qsub Variants calling ####
if ($step<=3 && $cut>=3){
        $cmd="cd $shell_dir/$sample/3.Variants/snp && $qsub -l p=1 -l vf=10g   HaplotypeCaller_$sample.sh";
        print "$cmd\n";
        $cmd="cd $shell_dir/$sample/3.Variants/snp && $qsub -l p=1 -l vf=10g  UnifiedGenotyper_$sample.sh";
        print "$cmd\n";
        $cmd="cd $shell_dir/$sample/3.Variants/snp && $qsub -l p=1 -l vf=10g atlas2_$sample.sh";
        print "$cmd\n";
        $cmd="cd $shell_dir/$sample/3.Variants/indel && $qsub -l p=1 -l vf=10g indel_atlas_$sample.sh";
        print "$cmd\n";
        $cmd="cd $shell_dir/$sample/3.Variants/indel && $qsub -l p=4 -l vf=10g  indel_platypus_$sample.sh";
        print "$cmd\n";
        system("cd $shell_dir/$sample/3.Variants/snp && $qsub -l p=1 -l vf=10g   HaplotypeCaller_$sample.sh");
        system("cd $shell_dir/$sample/3.Variants/snp && $qsub -l p=1 -l vf=10g  UnifiedGenotyper_$sample.sh");
        system("cd $shell_dir/$sample/3.Variants/snp && $qsub -l p=1 -l vf=10g atlas2_$sample.sh");
        system("cd $shell_dir/$sample/3.Variants/indel && $qsub -l p=1 -l vf=10g indel_atlas_$sample.sh");
        system("cd $shell_dir/$sample/3.Variants/indel && $qsub -l p=5 -l vf=10g  indel_platypus_$sample.sh");
        $finish = 1;
        while($finish){
                sleep 60;
                if (-e "$shell_dir/$sample/3.Variants/snp/UnifiedGenotyper_$sample.mark" && -e "$shell_dir/$sample/3.Variants/snp/HaplotypeCaller_$sample.mark" && -e "$shell_dir/$sample/3.Variants/snp/atlas_$sample.mark" && -e "$shell_dir/$sample/3.Variants/indel/atlas_indel_$sample.mark" && -e "$shell_dir/$sample/3.Variants/indel/platypus_indel_$sample.mark"){
                        $finish = 0;
                }else{
                        $finish = 1;
                }
        }
          $cmd="cd $shell_dir/$sample/3.Variants/snp   && $qsub -l p=1 -l vf=1g overlap_snp.sh";
         print "$cmd\n";
          $cmd="cd $shell_dir/$sample/3.Variants/indel && $qsub -l p=1 -l vf=1g overlap_indel.sh";
         print "$cmd\n";         system("cd $shell_dir/$sample/3.Variants/snp   && $qsub -l p=1 -l vf=1g overlap_snp.sh");
        system("cd $shell_dir/$sample/3.Variants/indel && $qsub -l p=1 -l vf=1g overlap_indel.sh");
        $finish = 1;
        while($finish){
                sleep 60;
                if (-e "$shell_dir/$sample/3.Variants/snp/overlap_var_$sample.mark" && -e "$shell_dir/$sample/3.Variants/indel/overlap_var_$sample.mark"){
                        $finish = 0;
                }else{
                        $finish = 1;
                }
        }	
}

###qsub annotation#########
if($type !~ /Gynecological/){
	if ($step <=4 && $cut>=4){
		$cmd="cd $shell_dir/$sample/4.Annot && $qsub -l p=1 -l vf=12g annotate.indel.$sample.sh";
                print "$cmd\n";
		$cmd="cd $shell_dir/$sample/4.Annot && $qsub -l p=1 -l vf=12g annotate.snp.$sample.sh";
                print "$cmd\n";
		system("cd $shell_dir/$sample/4.Annot && $qsub -l p=1 -l vf=12g annotate.indel.$sample.sh");
		system("cd $shell_dir/$sample/4.Annot && $qsub -l p=1 -l vf=12g annotate.snp.$sample.sh");
		$finish = 1;
		while($finish){
			sleep 60;
			if (-e "$shell_dir/$sample/4.Annot/annotate_snp.$sample.mark" && -e "$shell_dir/$sample/4.Annot/annotate_indel.$sample.mark"){
				die "$sample\_snp.annovar.xls not exists!\n" unless (-e "$shell_dir/$sample/4.Annot/$sample\_snp.annovar.xls");
				die "$sample\_indel.annovar.xls not exists!\n" unless ( -e "$shell_dir/$sample/4.Annot/$sample\_indel.annovar.xls");
				$finish = 0;
			}else{
				$finish = 1;
			}
		}
		$cmd="cd $shell_dir/$sample/6.QC && $qsub -l p=1 -l vf=3g QC.anno.indel.$sample.sh";
		$cmd="cd $shell_dir/$sample/6.QC && $qsub -l p=1 -l vf=3g QC.anno.snp.$sample.sh";
		system("cd $shell_dir/$sample/6.QC && $qsub -l p=1 -l vf=3g QC.anno.indel.$sample.sh");
		system("cd $shell_dir/$sample/6.QC && $qsub -l p=1 -l vf=3g QC.anno.snp.$sample.sh");
	}
}
if ($step <=5 && $cut>=5){
		$cmd="cd $shell_dir/$sample/5.Interpretation/$type && $qsub -l p=1 -l vf=2g Interpretation_$sample.sh";
                 print "$cmd\n";
                 system("cd $shell_dir/$sample/5.Interpretation/$type && $qsub -l p=1 -l vf=2g Interpretation_$sample.sh");

}
