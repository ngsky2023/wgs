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
my $aa_file;  ## aa 20 names

my $fOut;

GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fin,
				"a:s"=>\$aa_file,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fin);

#===============================================================
# Default optional value 
#===============================================================

#my $odir = dirname($fOut);
#mkdir $odir unless (-d $odir) ;
$aa_file="/home/lizf/workdir/ExomeSingle/aa.list.csv" unless defined($aa_file);
my $finb=$fin.".bak";
`cp $fin $finb`;
$fOut=$fin unless defined($fOut);
#my $clinvar_hgvs="/home/lizf/workdir/ExomeSingle/4gene.list";
my $clinvar_hgvs="/home/lizf/workdir/AmpliseqBRCA/55gene.list";
#my $pubmedinfo="/home/lizf/workdir/BRCA_panel/4gene.hgvs.pmidout";
#===============================================================
# Global value
#===============================================================
my %aa = ();  ## All aa three char name 
my %aa123=(); ## turn aa 1 char name to 3 chars name 
my %nt=();
my %hgvsc = ();  ## all hgvsc 

#===============================================================
# Get Data
#===============================================================

load_aa_from_file($aa_file, \%aa,\%aa123);
load_hgvsc_from_file($clinvar_hgvs, \%hgvsc);
#print "==== aa names====\n";
#print Dumper %aa;
$nt{"A"}="腺嘌呤";
$nt{"C"}="胞嘧啶";
$nt{"G"}="鸟嘌呤";
$nt{"T"}="胸腺嘧啶";

#===============================================================
# Process
#===============================================================
open(IN,$finb);
open(OUT,">$fOut");
my $line=<IN>;
chomp($line);
my $info="";
my $pmid="";
print OUT $line,"\tPMID\tinfo\n";
while($line=<IN>)
{
     chomp($line);
     my @temp=split(/\t/,$line,);
     if(exists($hgvsc{$temp[0]}{$temp[1]}{"info"}))
     {
     if($hgvsc{$temp[0]}{$temp[1]}{"pmid"}){
     $pmid=$hgvsc{$temp[0]}{$temp[1]}{"pmid"};
   }
     else
     {$pmid="24234437"; }
     	}
     else
     {
     $pmid="25741868";	
     }
     $info=cexplain($temp[0],$temp[1]);
     if($temp[2]){
     	if($temp[2]=~/p.[a-zA-Z]{1}\d/)
     	{
     	$temp[2]=hgvp_aa1to3($temp[2]);
     	print $temp[2],"\n";
      } 
     	$info.=pexplain($temp[2]);
     	}
     print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$pmid,"\t",$info,"。\n";

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
     ### up boundary Ins
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)_(\d+)-(\d+)ins(\S)/){

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
     my ($hgvsp) = @_;
     my $explain="";
     ###SNP
     if($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
         my $altaa=$4;
        $explain="，导致其编码蛋白质(".$refNP.")的第".$paa."位氨基酸残基由".$aa{$refaa}."突变为".$aa{$altaa};
     }
     ###SNP and fs
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})fs$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
         my $altaa=$4;
       	$explain="，导致其编码蛋白质（".$refNP."）的第".$paa."位氨基酸残基由".$aa{$refaa}."突变为".$aa{$altaa}."，并发生移码突变";
     }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)fs$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
       	$explain="，导致其编码蛋白质（".$refNP."）的第".$paa."位氨基酸残基".$aa{$refaa}."突变，并发生移码突变";
     }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)del$/){
         my $refNP=$1;
         my $refaa1=$2;
         my $paa1=$3;
         my $refaa2=$4;
         my $paa2=$5;       	
        $explain="，导致其编码蛋白质（".$refNP."）的第".$paa1."位氨基酸残基".$aa{$refaa1}."至".$paa2."位氨基酸残基".$aa{$refaa2}."发生缺失突变";
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
     else
     {
     	   $explain="，蛋白变化未知";
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
$hgvp=~/(NP_.*):p.(.*)/;
my $id=$1;
my $info=$2;
my $info3="";
for (my $i=0;$i<length($info);$i++)
{
my $char=substr($info,$i,1);
if($char=~/[A-Z]/)
{$info3.=$aa123{$char};}
else
{$info3.=$char;}
}

$hgvp=$id.":p.".$info3;
return $hgvp;

}

sub USAGE {
	my $usage=<<"USAGE";
Program: explain.pl
Version: $version
Contact: lizhifeng\@berrygenomics.com
Description: explain hgvs report site info for clinvar and ACMG result 

Usage:
  Options:
  -i     <int>  input file
  -a     <str>  aa names list file
  -o     <str>  output file
  -h		Help

USAGE
	print $usage;
	exit;
}


