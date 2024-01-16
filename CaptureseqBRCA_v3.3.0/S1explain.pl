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
$nt{"A"}="������";
$nt{"C"}="�����";
$nt{"G"}="������";
$nt{"T"}="�������";

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
     print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$pmid,"\t",$info,"��\n";

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

         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds."λ�����".$nt{$refnt}."ͻ��Ϊ".$nt{$altnt};
     }
     ### promotor
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.-(\d+)(\S)>(\S)/ ){

     my $refNM=$1;
     my $offset=$2;
     my $refnt=$3;
     my $altnt=$4; 
     $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼�����������ε�".$offset."λ�����".$nt{$refnt}."ͻ��Ϊ".$nt{$altnt};
   }     
     ### 3' of the translation termination codon
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.\*(\d+)([ACGT])>([ACGT])/){

         my $refNM=$1;
         my $cds=$2;
         my $refnt=$3;
         my $altnt=$4;

         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼�����������ε�".$cds."λ�����".$nt{$refnt}."ͻ��Ϊ".$nt{$altnt};
     }  
     ### up boundary SNP
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)(\S)>(\S)/){

     my $refNM=$1;
     my $cds=$2;
     my $offset=$3;
     my $refnt=$4;
     my $altnt=$5;

     $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds."λ������ε�".$offset."λ�����".$nt{$refnt}."ͻ��Ϊ".$nt{$altnt};
    }  
     ### down boundary
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)\+(\d+)(\S)>(\S)/){

     my $refNM=$1;
     my $cds=$2;
     my $offset=$3;
     my $refnt=$4;
     my $altnt=$5;

 	  $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds."λ������ε�".$offset."λ�����".$nt{$refnt}."ͻ��Ϊ".$nt{$altnt};
    }
 
     ### del 2&more seq, 
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $del=$4;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λ���".$del."ȱʧͻ��";
     }
     ### del 2&more number
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del(\d+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $deln=$4;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λȱʧ".$deln."���";
     }  
     ### del 2&more but no number or seq
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λ���ȱʧ";
     }        
     ### del1 
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)del(\S)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $refnt=$3;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."λ���".$nt{$refnt}."ȱʧͻ��";
     }
     ### dup 1
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)dup(\S)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $refnt=$3;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."λ���".$nt{$refnt}."�ظ�����";
     }    
     ### dup 2&more
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)dup([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $dup=$4;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λ���".$dup."�ظ�����";
     }    
    ### ins1+fs
     elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)ins([ACGT])$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $insnt=$4;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λ���֮�����".$nt{$insnt};
     }    
    ### ins 2&more seq +fs
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)ins([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $insnt=$4;
         $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λ���֮�����".$insnt."����";  
   }
    ### ins 2&more number +fs 
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)ins(\d+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $insn=$4;
         $explain=$hgvsc."\t".$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λ���֮�����".$insn."�����";        	
        }
        
    ### del+ins
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)_(\d+)del([ACGT]+)ins([ACGT]+)$/){

         my $refNM=$1;
         my $cds1=$2;
         my $cds2=$3;
         my $del=$4;
         my $ins=$5;
         $explain=$hgvsc."\t".$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."��".$cds2."λ���ȱʧ".$del."���в�����".$ins."����";
        }

    ### up boundary Del
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)del(\S)$/){

     my $refNM=$1;
     my $cds=$2;
     my $offset=$3;
     my $refnt=$4;
     $explain=$hgvsc."\t".$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds."λ������ε�".$offset."λ���".$nt{$refnt}."ȱʧ";
    }
     ### up boundary Ins
    elsif($hgvsc=~/(NM_\d+\.\d+):c\.(\d+)-(\d+)_(\d+)-(\d+)ins(\S)/){

     my $refNM=$1;
     my $cds1=$2;
     my $offset1=$3;
     my $cds2=$4;
     my $offset2=$5;     
     my $ins=$6;
 	  $explain=$hgvsc."ͻ����".$gene."����".$refNM."ת¼����".$cds1."λ������ε�".$offset1."λ�����".$cds2."λ������ε�".$offset2."λ���֮�����".$ins."����";
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
        $explain="����������뵰����(".$refNP.")�ĵ�".$paa."λ������л���".$aa{$refaa}."ͻ��Ϊ".$aa{$altaa};
     }
     ###SNP and fs
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)([a-zA-Z]{3})fs$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
         my $altaa=$4;
       	$explain="����������뵰���ʣ�".$refNP."���ĵ�".$paa."λ������л���".$aa{$refaa}."ͻ��Ϊ".$aa{$altaa}."������������ͻ��";
     }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)fs$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;
       	$explain="����������뵰���ʣ�".$refNP."���ĵ�".$paa."λ������л�".$aa{$refaa}."ͻ�䣬����������ͻ��";
     }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)del$/){
         my $refNP=$1;
         my $refaa1=$2;
         my $paa1=$3;
         my $refaa2=$4;
         my $paa2=$5;       	
        $explain="����������뵰���ʣ�".$refNP."���ĵ�".$paa1."λ������л�".$aa{$refaa1}."��".$paa2."λ������л�".$aa{$refaa2}."����ȱʧͻ��";
        }
      elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)delins(.*)$/){
         my $refNP=$1;
         my $refaa1=$2;
         my $paa1=$3;
         my $refaa2=$4;
         my $paa2=$5;          	
         my $aa3=$6;
         my $aaseq="";
         $explain="����������뵰���ʣ�".$refNP."���ĵ�".$paa1."λ������л�".$aa{$refaa1}."��".$paa2."λ������л�".$aa{$refaa2}."ȱʧ������".$aaseq;
        	}
        	
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)=$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;  
        $explain="��ͨ��RNA����������뵰����(".$refNP.")�ĵ�".$paa."λ������л�".$aa{$refaa}."�������ı�";
      }
     elsif($hgvsp=~/(NP_\d+\.\d+):p\.([a-zA-Z]{3})(\d+)delins([a-zA-Z]{3})$/){
         my $refNP=$1;
         my $refaa=$2;
         my $paa=$3;  
         my $altaa=$4;
        $explain="����������뵰���ʣ�".$refNP."���ĵ�".$paa."λ������л���".$aa{$refaa}."ͻ��Ϊ".$aa{$altaa};
      }      
     else
     {
     	   $explain="�����ױ仯δ֪";
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


