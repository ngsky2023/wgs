#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use Encode;
my $BEGIN_TIME=time();
my $version="1.0.0";
my @Original_ARGV=@ARGV;
#######################################################################################
# ==============================================================
# Get Options
# ==============================================================

my $fin;  ## input clinvar hgvs list 
my $aachange;  
my $panel_n;
my $fOut;

GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fin,
				"n:i"=>\$panel_n,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fin && $fOut && $panel_n);

#===============================================================
# Default optional value 
#===============================================================

#my $odir = dirname($fOut);
#mkdir $odir unless (-d $odir) ;
my $aa_file="/share/work2/yangrutao/workdir/research/databases/aa.list.csv" ;
#my $finb=$fin.".bak";
#`cp $fin $finb`;
#$fOut=$fin unless defined($fOut);
#my $clinvar_hgvs="/home/lizf/workdir/AmpliseqBRCA/55gene.list";
my $clinvar_hgvs="/share/work2/yangrutao/workdir/research/databases/55gene.list";
#my $pubmedinfo="/home/lizf/workdir/BRCA_panel/4gene.hgvs.pmidout";
my $panel_genes = "/share/work2/yangrutao/workdir/research/databases/" . $panel_n . "panel_genename.cfg"; 
#===============================================================
# Global value
#===============================================================
my %aa = ();  ## All aa three char name 
my %aa123=(); ## turn aa 1 char name to 3 chars name 
my %nt=();
my %hgvsc = ();  ## all hgvsc 
my @gene_names;
#===============================================================
# Get Data
#===============================================================

load_aa_from_file($aa_file, \%aa,\%aa123);
load_hgvsc_from_file($clinvar_hgvs, \%hgvsc);
load_panel_gene_file($panel_genes,@gene_names);
#print "==== aa names====\n";
#print Dumper %aa;
$nt{"A"}="腺嘌呤";
$nt{"C"}="胞嘧啶";
$nt{"G"}="鸟嘌呤";
$nt{"T"}="胸腺嘧啶";

#===============================================================
# Process
#===============================================================

open(OUT,">$fOut");

my $info="";
my $pmid="";

#print OUT "Gene\tTrancsript\tHGVS(c.)\tHGVS(p.)\tType\tLevel\tPMID\tinfo\n";
print OUT "Sample\tResult\tGene\tTranscript\tExon\tHGVS(c.)\tHGVS(p.)\tZygosity\tCliSig\tPMID\tMutinfo\tPop\tReference\tFunction\tLSDB\n";
foreach  my $s (@gene_names) {
open(IN,$fin) or die $!;
my $line1=<IN>;
while($line1=<IN>)
{
     chomp($line1);
     my @temp=split(/\t/,$line1,);
     my $trans=$temp[11];
     my $allele=$temp[11] .":" . $temp[13];
     my $type="";
	my $popfreq;
	my $in_silico;
	my $lsdb;
	my $refer;
	$pmid = $temp[22];
     #print $allele,"\n";
     $info=cexplain($temp[7],$allele);
     $aachange="$temp[14]:" . $temp[15];
	#print "$aachange\n";
     $info.=pexplain($aachange,$temp[17],$temp[7]);
     if($temp[18]){
		if($temp[18]=~/^synonymous/){
			$info.=",为同义突变。"; $type="同义突变";
		}elsif($temp[18]=~/^nonsynonymous/){
			$info.=",为错义突变。";$type="错义突变";
		}elsif($temp[18]=~/^frameshift/){
     	   	#$info.=",为移码突变";
     	   	$type="移码突变";
		}elsif($temp[18]=~/^nonframeshift/){
			$info.=",为非移码突变。";$type="非移码突变";
		}elsif($temp[18]=~/^stop/){
			$info.=",为无义突变。";$type="无义突变";
		}
	}
	
	if ($temp[28] =~ /Popfreq\s+(1000g2015aug_eas)\s+(.*?)\s+\/\s+(esp6500siv2_ea)\s+(.*?)\s+\/\s+(ExAC_EAS)\s+(.*?)\s+\/\s+(PopFreqMax)\s+(.*?)\s+\;/) {
          my $pop_max = $8;
		my $max;
		if  ($8 < $2 || $8 < $4 || $8 < $6 ){
			my @num;
			push (@num,$2);
			push (@num,$4);
			push (@num,$6);
			foreach (@num) {
				$max = $_ if $_ > $max;
			}
			$pop_max = $max;
		}

		$popfreq = "人群频率信息为";
		if ($2 ne ".") {
			$popfreq .= "千人基因组中东亚人群（1000g2015aug_eas）频率$2，";
		}else{
			$popfreq .= "千人基因组中东亚人群（1000g2015aug_eas）频率无记录，";
		}
		if ($4 ne ".") {
			$popfreq .= "esp6500中东亚人群（esp6500siv2_ea）频率$4，";
		}else{
			$popfreq .= "esp6500中东亚人群（esp6500siv2_ea）频率无记录，";
		}
		if ($6 ne ".") {
			$popfreq .= "ExAC中东亚人群（ExAC_EAS）频率$6，";
		}else{
			$popfreq .= "ExAC中东亚人群（ExAC_EAS）频率无记录，";
		}
		if ($pop_max ne "."){
			$popfreq .= "三个数据库最大人群频率为$pop_max\。";
		}else{
			$popfreq .= "三个数据库最大人群频率无记录。";
		} 	
		if ($temp[28] =~/(BZ8K)\s+(.*?)\;/) {
			if ($2 ne "."){
				$popfreq .= "贝瑞和康遗传性肿瘤基因突变数据库(V1.0)频率$2。";
			}else{
				$popfreq .= "贝瑞和康遗传性肿瘤基因突变数据库(V1.0)频率无记录。";
			}
		}
     }
	if ($temp[28] =~ /silico\s+(.*?);/) {
          if ($1 >= 4) {
               $in_silico = "多种信息预测软件综合判断为有害突变。";
          }elsif ($1 <= -4){
			$in_silico = "多种信息预测软件综合判断为无害变异。";
		}else{
			$in_silico = "多种信息预测软件综合判断为意义未明变异。";
		}      
     }
	if ($pmid ne "." &&  $pmid ne "-"){
		$refer = "该变异有文献报道($pmid)。";
	}else{
		$refer = "该变异无相关文献报道。";
	}
	if ($temp[28] =~/clinvar\s+(\d+)\;/) {
           $lsdb= "ClinVar数据库中收录了该条记录(AlleleID:$1)。";
     }else{
		 $lsdb = "ClinVar数据库中无该条记录。";
	}
	
	#if ($temp[20] =~ /pathogenic/i  && $temp[8]+$temp[9] >50 &&  $temp[10] >0.2){
	if ($temp[20] =~ /pathogenic/i && $temp[7] eq $s && $temp[8]+$temp[9] >50 &&  $temp[10] >0.2){
	#if ($temp[20] =~ /pathogenic/i && $temp[7] eq $s ){
			
			print OUT $temp[0],"\t",$temp[1],"\t",$temp[7],"\t",$trans,"\t",$temp[12],"\t",$temp[13],"\t",$temp[15],"\t",$temp[16],"\t",$temp[20],"\t",$pmid,"\t",$info,"\t",$popfreq,"\t",$refer,"\t",$in_silico,"\t",$lsdb,"\n";
			#Sample\tResult\tGene\tTranscript\tExon\tHGVS(c.)\tHGVS(p.)\tZygosity\tCliSig\tPMID\tMutinfo\tPop\tReference\tFunction\tLSDB\n
			
		}	
	}
}     
close(IN);

foreach  my $s (@gene_names) {
open(IN,$fin) or die $!;
my $line2=<IN>;
while($line2=<IN>)
{
     chomp($line2);
     my @temp=split(/\t/,$line2);
     my $trans=$temp[11];
     my $allele=$temp[11] .":" . $temp[13];
	
     my $type="";
	my $popfreq;
	my $in_silico;
	my $lsdb;
	my $refer;
	$pmid = $temp[22];
     #print $allele,"\n";
     $info=cexplain($temp[7],$allele);
     $aachange="$temp[14]:" . $temp[15];
     $info.=pexplain($aachange,$temp[17],$temp[7]);
     
     if($temp[18]){
		if($temp[18]=~/^synonymous/){
			$info.=",为同义突变。"; $type="同义突变";
		}elsif($temp[18]=~/^nonsynonymous/){
			$info.=",为错义突变。";$type="错义突变";
		}elsif($temp[18]=~/^frameshift/){
     	   	#$info.=",为移码突变";
     	   	$type="移码突变";
		}elsif($temp[18]=~/^nonframeshift/){
			$info.=",为非移码突变。";$type="非移码突变";
		}elsif($temp[18]=~/^stop/){
			$info.=",为无义突变。";$type="无义突变";
		}
	}

	
	if ($temp[28] =~ /Popfreq\s+(1000g2015aug_eas)\s+(.*?)\s+\/\s+(esp6500siv2_ea)\s+(.*?)\s+\/\s+(ExAC_EAS)\s+(.*?)\s+\/\s+(PopFreqMax)\s+(.*?)\s+\;/) {
          my $pop_max = $8;
		my $max = 0;
		my $a2=$2;my $a4=$4;my $a6=$6;
		$a2 = 0 if ($a2 eq ".");
		$a4 = 0 if ($a4 eq ".");
		$a6 = 0 if ($a6 eq ".");
		#print "********$a2\t$a4\t$a6\t$8\n";
		if  ($8 < $a2 || $8 < $a4 || $8 < $a6 ){
			my @num;
			push (@num,$a2);
			push (@num,$a4);
			push (@num,$a6);
			foreach (@num) {
				$max = $_ if ($_ > $max);
			}
			$pop_max = $max;
		}
		if($8 > $a2 && $8 > $a4 && $8 > $a6){
			my @num;
			push (@num,$a2);
			push (@num,$a4);
			push (@num,$a6);
			foreach (@num) {
				$max = $_ if ($_ > $max);
			}
			$pop_max = $max;
			#print "********$a2\t$a4\t$a6\t$8\t$max\n";
		}

		$popfreq = "人群频率信息为";
		if ($2 ne ".") {
			$popfreq .= "千人基因组中东亚人群（1000g2015aug_eas）频率$2，";
		}else{
			$popfreq .= "千人基因组中东亚人群（1000g2015aug_eas）频率无记录，";
		}
		if ($4 ne ".") {
			$popfreq .= "esp6500中东亚人群（esp6500siv2_ea）频率$4，";
		}else{
			$popfreq .= "esp6500中东亚人群（esp6500siv2_ea）频率无记录，";
		}
		if ($6 ne ".") {
			$popfreq .= "ExAC中东亚人群（ExAC_EAS）频率$6，";
		}else{
			$popfreq .= "ExAC中东亚人群（ExAC_EAS）频率无记录，";
		}
		if ($pop_max ne "."){
			$popfreq .= "三个数据库最大人群频率为$pop_max\。";
		}else{
			$popfreq .= "三个数据库最大人群频率无记录。";
		} 	
		if ($temp[28] =~/(BZ8K)\s+(.*?)\;/) {
			if ($2 ne "."){
				$popfreq .= "贝瑞和康遗传性肿瘤基因突变数据库(V1.0)频率$2。";
			}else{
				$popfreq .= "贝瑞和康遗传性肿瘤基因突变数据库(V1.0)频率无记录。";
			}
		}
     }
	if ($temp[28] =~ /silico\s+(.*?);/) {
          if ($1 >= 4) {
               $in_silico = "多种信息预测软件综合判断为有害突变。";
          }elsif ($1 <= -4){
			$in_silico = "多种信息预测软件综合判断为无害变异。";
		}else{
			$in_silico = "多种信息预测软件综合判断为意义未明变异。";
		}      
     }
	if ($pmid ne "." &&  $pmid ne "-"){
		$refer = "该变异有文献报道($pmid)。";
	}else{
		$refer = "该变异无相关文献报道。";
	}
	if ($temp[28] =~/clinvar\s+(\d+)\;/) {
           $lsdb= "ClinVar数据库中收录了该条记录(AlleleID:$1)。";
     }else{
		 $lsdb = "ClinVar数据库中无该条记录。";
	}

	if ($temp[20] =~ /Uncertain significance/ && $temp[7] eq $s && $temp[8]+$temp[9] >50 &&  $temp[10] >0.2){
	#if ($temp[20] =~ /Uncertain significance/ && $temp[7] ~~ @gene_names ){
	#if ($temp[20] =~ /Uncertain significance/ && $temp[7] eq $s ){
		
		print OUT $temp[0],"\t",$temp[1],"\t",$temp[7],"\t",$trans,"\t",$temp[12],"\t",$temp[13],"\t",$temp[15],"\t",$temp[16],"\t",$temp[20],"\t",$pmid,"\t",$info,"\t",$popfreq,"\t",$refer,"\t",$in_silico,"\t",$lsdb,"\n";
	}

	}
}     
close(IN);

close(OUT);

sub cexplain{

     my ($gene, $hgvsc) = @_;
     my $explain="";
     ##### SNP
     if($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)([ACGT])>([ACGT])/){

         my $refNM=$1;
         my $cds=$2;
         my $refnt=$3;
         my $altnt=$4;

         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds."位碱基由".$nt{$refnt}."突变为".$nt{$altnt};
     }
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.([ACGT])(\d+)([ACGT])/){

         my $refNM=$1;
         my $refnt=$2;
         my $cds=$3;
         my $altnt=$4;

         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds."位碱基由".$nt{$refnt}."突变为".$nt{$altnt};
     }
     ### promotor
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.-(\d+)(\S)>(\S)/ ){

     my $refNM=$1;
     my $offset=$2;
     my $refnt=$3;
     my $altnt=$4; 
     $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本编码区上游第".$offset."位碱基由".$nt{$refnt}."突变为".$nt{$altnt};
   }     
     ### 3' of the translation termination codon
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.\*(\d+)([ACGT])>([ACGT])/){

         my $refNM=$1;
         my $cds=$2;
         my $refnt=$3;
         my $altnt=$4;

         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本编码区下游第".$cds."位碱基由".$nt{$refnt}."突变为".$nt{$altnt};
     }  
     ### up boundary SNP
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)(\S)>(\S)/){

     my $refNM=$1;
     my $cds=$2;
     my $offset=$3;
     my $refnt=$4;
     my $altnt=$5;

     $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds."位碱基上游第".$offset."位碱基由".$nt{$refnt}."突变为".$nt{$altnt};
    }  
     ### down boundary
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)\+(\d+)(\S)>(\S)/){

     my $refNM=$1;
     my $cds=$2;
     my $offset=$3;
     my $refnt=$4;
     my $altnt=$5;

 	  $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds."位碱基下游第".$offset."位碱基由".$nt{$refnt}."突变为".$nt{$altnt};
    }
 
     ### del 2&more seq, 
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $del=$4;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."至".$cds2."位碱基".$del."缺失突变";
     }
     ### del 2&more number
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del(\d+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $deln=$4;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."至".$cds2."位缺失".$deln."碱基";
     }  
     ### del 2&more but no number or seq
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."至".$cds2."位碱基缺失";
     }        
     ### del1 
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)del(\S)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $refnt=$3;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."位碱基".$nt{$refnt}."缺失突变";
     }
     ### dup 1
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)dup(\S)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $refnt=$3;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."位碱基".$nt{$refnt}."重复插入";
     }    
     ### dup 2&more
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)dup([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $dup=$4;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."至".$cds2."位碱基".$dup."重复插入";
     }    
    ### ins1+fs
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)ins([ACGT])$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $insnt=$4;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."和".$cds2."位碱基之间插入".$nt{$insnt};
     }    
    ### ins 2&more seq +fs
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)ins([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $insnt=$4;
         $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."和".$cds2."位碱基之间插入".$insnt."序列";  
   }
    ### ins 2&more number +fs 
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)ins(\d+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $insn=$4;
         $explain=$hgvsc."\t".$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."和".$cds2."位碱基之间插入".$insn."个碱基";        	
        }
        
    ### del+ins
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del([ACGT]+)ins([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $del=$4;
         my $ins=$5;
         $explain=$hgvsc."\t".$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."至".$cds2."位碱基缺失".$del."序列并插入".$ins."序列";
        }

    ### up boundary Del
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)del(\S)$/){

     my $refNM=$1;
     my $cds=$2;
     my $offset=$3;
     my $refnt=$4;
     $explain=$hgvsc."\t".$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds."位碱基上游第".$offset."位碱基".$nt{$refnt}."缺失";
    }
    
    ### up boundary Del
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)_(\d+)-(\d+)del(\S+)/){

     my $refNM=$1;
     my $cds1=$2;
     my $offset1=$3;
     my $cds2=$4;
     my $offset2=$5;     
     my $del=$6;
 	  $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."位碱基上游第".$offset1."位碱基至".$cds2."位碱基上游第".$offset2."位碱基之间缺失".$del."序列";
    }
    
     ### up boundary Ins
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)_(\d+)-(\d+)ins(\S+)/){

     my $refNM=$1;
     my $cds1=$2;
     my $offset1=$3;
     my $cds2=$4;
     my $offset2=$5;     
     my $ins=$6;
 	  $explain=$hgvsc."突变是".$gene."基因".$refNM."转录本第".$cds1."位碱基上游第".$offset1."位碱基至".$cds2."位碱基上游第".$offset2."位碱基之间插入".$ins."序列";
    }

return $explain;
}

sub pexplain{
     my ($hgvsp,$ref_func,$pgene) = @_;
     my $explain="";
     ###SNP
     if($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
         my $altaa=$4;
         if($refaa eq $altaa)
         {
           $explain="，其编码蛋白质(".$refNP.")的第".$paa."位氨基酸残基".$aa{$refaa}."不发生改变";
         	 }
         else{
           $explain="，导致其编码蛋白质(".$refNP.")的第".$paa."位氨基酸残基由".$aa{$refaa}."突变为".$aa{$altaa};
         }
     }
     ###SNP and fs
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})fs$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
         my $altaa=$4;
       	$explain="，导致其编码蛋白质（".$refNP."）的第".$paa."位氨基酸残基由".$aa{$refaa}."突变为".$aa{$altaa}."，并发生移码突变。";
     }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)fs$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
       	$explain="，导致其编码蛋白质（".$refNP."）的第".$paa."位氨基酸残基".$aa{$refaa}."突变，并发生移码突变。";
     }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)del$/){
         my $refNP=$1;
         my $refaa1=$2;
         my $paa1=$3;
         my $refaa2=$4;
         my $paa2=$5;       	
        $explain="，导致其编码蛋白质（".$refNP."）的第".$paa1."位氨基酸残基".$aa{$refaa1}."至".$paa2."位氨基酸残基".$aa{$refaa2}."发生缺失突变。";
        }
      elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)delins(.*)$/){
         my $refNP=$1;
         my $refaa1=$2;
         my $paa1=$3;
         my $refaa2=$4;
         my $paa2=$5;          	
         my $aa3=$6;
         my $aaseq="";
         $explain="，导致其编码蛋白质（".$refNP."）的第".$paa1."位氨基酸残基".$aa{$refaa1}."至".$paa2."位氨基酸残基".$aa{$refaa2}."缺失，插入".$aaseq;
        	}
        	
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)=$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;  
        $explain="，通过RNA分析，其编码蛋白质(".$refNP.")的第".$paa."位氨基酸残基".$aa{$refaa}."不发生改变";
      }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)delins([a-zA-Z]{3})$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;  
         my $altaa=$4;
        $explain="，导致其编码蛋白质（".$refNP."）的第".$paa."位氨基酸残基由".$aa{$refaa}."突变为".$aa{$altaa};
      }
	elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)delins([a-zA-Z]{3})([a-zA-Z]{3})$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;  
         my $altaa1=$4;
	    my $altaa2=$5;
        $explain="，导致其编码蛋白质（".$refNP."）的第".$paa."位氨基酸残基由".$aa{$refaa}."突变为".$aa{$altaa1}."和".$aa{$altaa2};
      }
	elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)(Ter)$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;  
         my $altaa1=$4;
        $explain="，导致其编码蛋白质（".$refNP."）的第".$paa."位氨基酸残基由".$aa{$refaa}."突变为终止子";
      }
     else
     {
		if ($ref_func =~/splicing/) {
			$explain = "，该突变会导致$pgene\基因mRNA的异常剪接，可能会影响蛋白质功能";
          }
		$explain .="，蛋白变化未知。";
     	} 
return $explain;
}



sub load_aa_from_file{

my ($faa, $aaname,$aa123) = @_;
open(FAA,$faa);
my $line=<FAA>;
while($line=<FAA>)
{
	my @temp=split(/,/,$line,);
	$aaname->{$temp[0]}=$temp[2];
	$aa123->{$temp[1]}=$temp[0];
}
close(FAA);

}

sub load_hgvsc_from_file{

my ($fhgvs, $hgvsc) = @_;
open(FAA,$fhgvs);
my $line=<FAA>;
while($line=<FAA>)
{
	chomp($line);
	my @temp=split(/\t/,$line,);
	$hgvsc->{$temp[0]}{$temp[1]}{"info"}=$line;
	$hgvsc->{$temp[0]}{$temp[1]}{"pmid"}=$temp[3];
}
close(FAA);

}

sub hgvp_aa1to3{
	my ($hgvp)=@_;
	print $hgvp,"\n";
	$hgvp=~/p.(.*)/;
	my $info=$1;
	my $info3="";
	for (my $i=0;$i<length($info);$i++)
	{
		my $char=substr($info,$i,1);
		if($char=~/[A-Z]/)
		{$info3.=$aa123{$char};}
		else
		{$info3.=$char;}
	}
	$hgvp="p.".$info3;
	return $hgvp;
}

sub load_panel_gene_file{
	my ($cfg_55, $gene_name) = @_;
	open(PAN,$cfg_55);
	while(<PAN>)
	{
		chomp;
		next if ($_ =~/^\s+$/);
		unless ($_ ~~ @gene_names){
			push (@gene_names,$_ );
		}
	}
	close(PAN);
}
			
sub USAGE {
	my $usage=<<"USAGE";
Program: explain.pl
Version: $version
Contact: lizhifeng\@berrygenomics.com
Modified:yangrutao

Description: explain hgvs report site info for clinvar and ACMG result 

Usage:
  Options:
  -i     <int>  input file
  -n     <num>  genes number of panel 
  -o     <str>  output file
  -h		Help

USAGE
	print $usage;
	exit;
}


