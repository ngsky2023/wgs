use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Data::Dumper;

sub Usage{
        print STDERR <<USAGE;
=====================================================================
Description: Rating the variants after annotating;
Options:
        -i <s>        input file;
        -f1 [f]       MAF cutoff;
        -f2 [f]       disorder cutoff;
        -o [s]        Output file name;
        -h            Help information;
Author: zhangyu; zhangyu001\@berrygenomics.com
Version: V1.1
Date: 2016.03.30
Author: lizhifeng; lizhifeng\@berrygenomics.com
Version: V2.1
Date: 2016.09.10
Version: V2.2
Date: 2016.12.26
=====================================================================
USAGE
}

my ($in,$f1,$f2,$out,$help);
GetOptions(
        "h|?|help"=>\$help,
        "in=s"=>\$in,
        "f1=f"=>\$f1,
        "f2=f"=>\$f2,
        "o=s"=>\$out,
);
if(!defined($in)|| defined($help)){
        &Usage;
        exit 0;
}
$out ||="out";
$f1 ||="0.01";
$f2 ||="0.5";

#=====================Configure===================================
my $gedi="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/Gene_Disorder.xls";
my $gediva="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/Gene_Disorder_AA2_expand_20180326.txt";  
my $ge="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/Gene.xls";
my $clinvar= "/share/public/database/Gynecological_cancer_backup/IPDB/ver4/variant_summary_cita_GRCh37_NM_20180326.txt"; 
my $assay="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/BRCA_AssaySummary.txt";
my $NMdb="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/variant_summary.GRCh37.gene_MainNM_ClinVar20171031.txt";
my $ACMGpre="/share/public/database/Gynecological_cancer_backup/IPDB/ver4/ACMG_pre_20180629.txt";
my $BZ8K="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/BZ8K.forGD.txt";
my $BRCA1400K="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/BRCA1.400K.vcf.annovout";
my $BRCA2400K="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/BRCA2.400K.vcf.annovout";
my $dup="/share/public/database/Gynecological_cancer_backup/IPDB/ver2/RepeatRegion.position";
#=================================================================

open IN,"<$in";
open OUT,">$out";

open DB,"<$clinvar";
my (%hash,%hashv,%hashA);

while(<DB>){
	if($_=~/^#/){next;}      #yangrt20170906
	chomp;
	my @tmp=split(/\t/,$_);
	my($chr,$start,$end,$ref,$alt)=($tmp[18],$tmp[19],$tmp[20],$tmp[21],$tmp[22]);  #yangrt  20170901
	my $clin=join("&",($tmp[6],$tmp[12],$tmp[30]));                                 
	$chr="chr".$tmp[18];                                                           
	my $string="";
	if($tmp[1]=~/insertion/){
                $string=join("\t",($chr,$tmp[19],$tmp[19],$tmp[21],$tmp[22]));           #yangrt  20170901
	}
	else{
                $string=join("\t",($chr,$tmp[19],$tmp[20],$tmp[21],$tmp[22]));                       #yangrt  20170901
        }
	$hash{$string}=$clin; 
        $hashA{$string}=$tmp[0];
        my ($AlleleID,$ClinicalSignificance,$reviewstatus,$Citation_PubMed)=($tmp[0],$tmp[6],$tmp[24],$tmp[30]);     #yangrt  20170901
        #$hashv{$AlleleID}=[$ClinicalSignificance,$reviewstatus,$Citation_PubMed];
        $hashv{$AlleleID}->[0].=$ClinicalSignificance;
        $hashv{$AlleleID}->[1].=$reviewstatus;
        $hashv{$AlleleID}->[2].=$Citation_PubMed;
   #     if($AlleleID== 24364 || $AlleleID== 97222)
   #     {print "$AlleleID ClinicalSignificance $hashv{$AlleleID}->[0] reviewstatus $hashv{$AlleleID}->[1] pubmed $hashv{$AlleleID}->[2] \n"; }
}
close DB;
my (%hashAY);
open(ASY,"$assay");
while(my $line=<ASY>)
{
        if($line!~/^BRCA/){next;}
        else
        {
                chomp($line);
                my @tmp=split(/\t/,$line);
                if($line=~/Deleterious/){
                        $hashAY{$tmp[0]}{$tmp[1]}="D";
                        if($tmp[2]=~/c./){$hashAY{$tmp[0]}{$tmp[2]}="D";}
                }elsif($line=~/Neutral/){
                        $hashAY{$tmp[0]}{$tmp[1]}="N";
                        if($tmp[2]=~/c./){$hashAY{$tmp[0]}{$tmp[2]}="N";}
                }
        }
}
close(ASY);

my %hashg;
open GE, "<$ge";
while(<GE>){
	chomp;
	$hashg{$_}=0;
}
close GE;

my (%hashp,%hashc,%hasha,%hasht,%hashr,%hashap,%conf_clinsig,%hashp_clinsig,%hashp_revsta,%hashc_clinsig,%hashc_revsta,%hashtalleleid,%hashtp,%hasht_CliSigSim,);      #yangrt  20170901
open GV,"<$gediva";
while(<GV>){
	chomp;
	my @tmp=split /\t/;   
	#$tmp[0]=="gene";$tmp[1]=="disname/alleleID";$tmp[2]=="p.AAchange";$tmp[3]=="c.hgvs";{$tmp[4]}{$tmp[5]}{$tmp[6]}=="AAref,pos,AAalt"
        my $nm;
	if($tmp[3]){
                if($tmp[3]=~/(NM_\d+).*(c.*)/)
		{
                        $hashc{$tmp[0]}{$1}{$2}=$tmp[1];
                        $hashc_clinsig{$tmp[0]}{$1}{$2}=$tmp[7];        #yangrt  20170901
                        $hashc_revsta{$tmp[0]}{$1}{$2}=$tmp[9];         #yangrt  20170901
                        $nm=$1;
                }
	}
        $hashp{$tmp[0]}{$nm}{$tmp[2]}=$tmp[1];
        $hashp_clinsig{$tmp[0]}{$nm}{$tmp[2]}=$tmp[7];      #yangrt  2017091
        $hashp_revsta{$tmp[0]}{$nm}{$tmp[2]}=$tmp[9];       #yangrt  2017091
	if($tmp[5] ne "-" ){                                            #yangrt  20170901              
	$hasha{$tmp[0]}{$nm}{$tmp[4]}{$tmp[5]}=$tmp[6];
        if($tmp[1]=~/^\d+$/)
        {
          if(exists($hashv{$tmp[1]})){
                $hasht{$tmp[0]}{$nm}{$tmp[4]}{$tmp[5]}=$tmp[7];           #yangrt  20170901
                $hasht_CliSigSim{$tmp[0]}{$nm}{$tmp[4]}{$tmp[5]}=$tmp[8];
                $hashtp{$tmp[0]}{$nm}{$tmp[4]}{$tmp[5]}=$tmp[13]; 
                $hashtalleleid{$tmp[0]}{$nm}{$tmp[4]}{$tmp[5]}=$tmp[1];
                #$hashr{$tmp[0]}{$tmp[4]}{$tmp[5]}=$tmp[9];          #yangrt   20170210
                #$hasht{$tmp[0]}{$tmp[4]}{$tmp[5]}.=$hashv{$tmp[1]}->[0]."; ";
                #$hashr{$tmp[0]}{$tmp[4]}{$tmp[5]}.=$hashv{$tmp[1]}->[1]."; ";
         }
    #     else {print "AlleleID $tmp[1] has no reviewstatus\n";}
          }
	}
#         if($tmp[1]=~/^\d+$/ && ($tmp[1] == 152550 || $tmp[1] == 186542 ||$tmp[1] == 46142  ))
#        {print "$tmp[0] $tmp[4] $tmp[5] $tmp[2] $hashp{$tmp[0]}{$tmp[2]}  $hashv{$tmp[1]}->[0] $hashr{$tmp[0]}{$tmp[4]}{$tmp[5]} \n"; }
        if ($tmp[10] ne "-"){                                       #yangrt   20170210
                $conf_clinsig{$tmp[0]}{$tmp[1]}=$tmp[10];              
        }
}
close GV;

my %hashnm;
my %hashnmid;
open NM,"<$NMdb";
while(<NM>){
	chomp;
	my @tmp=split /\t/;
        if ($tmp[0] eq "MRE11"){
                $tmp[0] = "MRE11A";
        }
	$hashnm{$tmp[0]}{$tmp[1]}+=1;
        $hashnmid{$tmp[0]}=$tmp[1];
        #if($hashnm{$tmp[0]}{$tmp[1]}>1){print $tmp[0],"\t",$tmp[1],"\n";}
}
close NM;

my (%pvs1_LOF,%ref_NM,%pm1_domain,%pp2_mistype);
open PRE,"<$ACMGpre";
<PRE>;
while(<PRE>)
{
     chomp;
     my @tmp=split /\t/;
     $pvs1_LOF{$tmp[0]}=$tmp[2];
     $pp2_mistype{$tmp[0]}=$tmp[5];
     $ref_NM{$tmp[0]}=$tmp[1];
     if($tmp[4])
     {
     $pm1_domain{$tmp[0]}=$tmp[6];
     }
}
close(PRE);

my %hashf;
open(POP,"<$BZ8K");
<POP>;
while(<POP>)
{
    chomp;
    my @tmp=split /\t/;
    my($chr,$start,$end,$ref,$alt,$freq)=($tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4],$tmp[92]);
    my $string=join("\t",($tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4]));
    $hashf{$string}=$freq;    
}
close(POP);
my %f400K;
open(POP,"<$BRCA1400K");
<POP>;
while(<POP>)
{
    chomp;
    my @tmp=split /\t/;
    my $freq=$tmp[5];
    my $string=join("\t",($tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4]));
    $f400K{$string}=$freq;  
}
close(POP);
open(POP,"<$BRCA2400K");
<POP>;
while(<POP>)
{
    chomp;
    my @tmp=split /\t/;
    my $freq=$tmp[5];
    my $string=join("\t",($tmp[0],$tmp[1],$tmp[2],$tmp[3],$tmp[4]));
    $f400K{$string}=$freq;
}
close(POP);

my %dupreg;
open(DUP,"<$dup");
while(<DUP>)
{
chomp;
    my @tmp=split /\t/;
    my $index=int($tmp[1]/1000000);
    $dupreg{$tmp[0]}{$index}{$tmp[1]}=$tmp[2];
}
close(DUP);

my $ngene="BRCA1";
#print $ngene,"\t",$hashg{$ngene},"\t", $pvs1_LOF{$ngene},"\t", $pvs1_LOF{$ngene},"\t",$hashnmid{$ngene},"\n";
$ngene="BRCA2";
#print $ngene,"\t",$hashg{$ngene},"\t", $pvs1_LOF{$ngene},"\t", $pvs1_LOF{$ngene},"\t",$hashnmid{$ngene},"\n";


my $title=<IN>;
chomp $title;
print OUT "Evidence\t$title\tTranscript\tExon\thgvs.c\thgvs.p\tmark\n";

while(<IN>){       
        chomp;
	my @atm=split /\t/;
	my @ACMGrule=();
# 0 pvs1 1 ps1 2 ps2 3 ps3 4 ps4 5 pm1 6 pm2 7 pm3 8 pm4 9 pm5 10 pm6 11 pp1 12 pp2 13 pp3 14 pp4 15 pp5 
#16 ba1 17 bs1 18 bs2 19 bs3 20 bs4 21 bp1 22 bp2 23 bp3 24 bp4 25 bp5 26 bp6 27 bp7  
	my ($ratio,$func,$gene,$genedetail,$exon_func)=($atm[9],$atm[11],$atm[12],$atm[13],$atm[14]);
	my ($chr,$start,$end,$ref,$alt,$AAchange)=($atm[1],$atm[2],$atm[3],$atm[4],$atm[5],$atm[15]);
        my ($intp_dormain,$dpsi_score) = ($atm[111],$atm[124]);
        if($AAchange=~/^\./ && $genedetail=~/NM/){$AAchange=$genedetail;}
	my $flag=0;
	my $cita;
	my $mark="";
	my $string=join("\t",($chr,$start,$end,$ref,$alt));
        if ($gene =~ /^(.*);.*/){     #yangrt   20170901
                $gene = $1;   
        }
#	print "=============start $string ==============\n";
#pvs1 & pm4
	if($func =~ /^exonic/){
		if($exon_func =~ (/stopgain|^frameshift/)){
			$flag=1;
		}elsif($exon_func =~ /nonframeshift insertion|nonframeshift deletion|stoploss/){
                          $ACMGrule[8]=1; # pm4 ; repeat filter
		          my $index=int($start/1000000);
                          foreach my $site (keys %{$dupreg{$chr}{$index}})
                          {
                            if($start>=$site && $end<=$dupreg{$chr}{$index}{$site} )
                            {$ACMGrule[8]=0; $ACMGrule[23]=1;}# bp3 ;  repeat filter
                          }	

		}
		elsif($exon_func=~/SNV/ && exists($pp2_mistype{$gene}))
		{if(!$pp2_mistype{$gene} && $exon_func!~/nonsynonymous/)
			{
				$ACMGrule[21]=1; #bp1
			}
		}
		
	}elsif($func =~ /splicing/){
		$flag=1;
		$mark.="splicing;";
	}
	if($flag == 1 && exists $hashg{$gene} && exists $pvs1_LOF{$gene} && $pvs1_LOF{$gene}=~/YES/ && ($AAchange=~/$hashnmid{$gene}/ || $genedetail=~/$hashnmid{$gene}/)){
		$ACMGrule[0]=1; #pvs1
                if($gene eq "BRCA2" && $start>=32972626){$ACMGrule[0]=0;} 
	}elsif($flag == 1 && exists $hashg{$gene} &&  not  exists $pvs1_LOF{$gene}){
                 $ACMGrule[0]=1; #pvs1
         }

#ps1 & pm5
#        print "genedetail $genedetail AAchange $AAchange\n";
#	$AAchange=~ s/"//g;
	my ($NM,$en,$c,$p)=(".",".",".",".");
	my ($trans,$exon,$hgvsc,$hgvsp)=(".",".",".",".");
        if(exists $hashA{$string}){
        #	print "exists hashA $string $hashv{$hashA{$string}}->[0] $hashv{$hashA{$string}}->[1] \n";
        	$mark.="clinvar $hashA{$string};";
        	if($hashv{$hashA{$string}}->[0] =~ /Pathogenic/i){
        #		print "class $hashv{$hashA{$string}}->[0] status $hashv{$hashA{$string}}->[1] \n";
        		if($hashv{$hashA{$string}}->[1]=~/criteria provided, multiple submitters, no conflicts|criteria provided, single submitter|practice guideline|reviewed by expert panel/)
                        {
                        	$ACMGrule[1]=1; #ps1
                        #	$ACMGrule[4]=1; #ps4
                        	$cita=$hashv{$hashA{$string}}->[1]; 
                        }
                        elsif($hashv{$hashA{$string}}->[1]=~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashv{$hashA{$string}}->[0] =~ /Pathogenic/)
                        {
                        	$ACMGrule[15]=1;#pp5
                        	$cita=$hashv{$hashA{$string}}->[1];
                        }
                
                        #if ($conf_clinsig{$gene}{$hashA{$string}} =~ /Pathogenic/i){        #yangrt   20170120
                        #        $ACMGrule[1]=1; #ps1
                        #        $cita=$hashv{$hashA{$string}}->[1];        
                        #}
                }elsif($hashv{$hashA{$string}}->[0] =~ /Benign/i){
                       #print "class $hashv{$hashA{$string}}->[0] status $hashv{$hashA{$string}}->[1] \n";
                        if($hashv{$hashA{$string}}->[1]=~/criteria provided, multiple submitters, no conflicts|criteria provided, single submitter|practice guideline|reviewed by expert panel/)
                        {
                        	$ACMGrule[18]=1;#bs2
                        	$ACMGrule[26]=1;#bp6
                        	$cita=$hashv{$hashA{$string}}->[1];
                        } 
                        elsif($hashv{$hashA{$string}}->[1]=~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashv{$hashA{$string}}->[0] =~ /Benign/)
                        {                       
                        	$ACMGrule[26]=1;#bp6
                        	$cita=$hashv{$hashA{$string}}->[1];
                        }
                
                        #if ($conf_clinsig{$gene}{$hashA{$string}} =~ /Benign/){        #yangrt   20170120
                        #        $ACMGrule[18]=1;#bs2
                        #        $ACMGrule[26]=1;#bp6
                        #        $cita=$hashv{$hashA{$string}}->[1]; 
                        #}
                        #elsif ($conf_clinsig{$gene}{$hashA{$string}} =~ /Likely benign/){        #yangrt   20170120
                        ##        #$ACMGrule[18]=1;#bs2
                        #        $ACMGrule[26]=1;#bp6
                        #        $cita=$hashv{$hashA{$string}}->[1]; 
                        #}
                }
                
                if ($conf_clinsig{$gene}{$hashA{$string}} =~ /Pathogenic/i){        #yangrt   20170120
                                $ACMGrule[1]=1; #ps1
                                $cita=$hashv{$hashA{$string}}->[1];        
                }
                if ($conf_clinsig{$gene}{$hashA{$string}} =~ /Benign/){        #yangrt   20170120
                                $ACMGrule[18]=1;#bs2
                                $ACMGrule[26]=1;#bp6
                                $cita=$hashv{$hashA{$string}}->[1]; 
                }elsif ($conf_clinsig{$gene}{$hashA{$string}} =~ /Likely benign/){        #yangrt   20170120
                        #        #$ACMGrule[18]=1;#bs2
                        $ACMGrule[26]=1;#bp6
                        $cita=$hashv{$hashA{$string}}->[1]; 
                }
                
                
        }  

	if(exists($ref_NM{$gene}) && $AAchange ne "."){
                my @tmp=split(/,/,$AAchange,);
                for (my $i=0;$i<=$#tmp;$i++){
                	#BRCA1:NM_007298:exon1:c.G71A:p.C24Y	
                        my @exp=split(/:/,$tmp[$i],);
                        #    print "$i $tmp[$i]\n";
                        for (my $l=0;$l<=$#exp;$l++){
                                $NM=$exp[$l] if ($exp[$l]=~ /NM/);
                                $en=$exp[$l] if ($exp[$l]=~ /exon/);
                                $c=$exp[$l] if ($exp[$l]=~ /c\./);
                                $p=$exp[$l] if ($exp[$l]=~ /p\./);
                        }
                        #  print "NM $NM exon $en CDS $c AA $p\n";
  
                        if($NM=~/NM/ && $ref_NM{$gene}=~/$NM/)
                        {
                                #  print "$gene main_tran $ref_NM{$gene} $NM\n";
                                # ps3/bs3
                                if($gene=~/BRCA/){   
                                        if($p=~/p.(\S+)/){
                                                my $paa=$1;
                                                if(exists($hashAY{$gene}{$paa}))
                                                {
                                                        #print "Function Assay ",$hashAY{$gene}{$paa},"\n";
                                                        if($hashAY{$gene}{$paa}=~/D/ && $hashAY{$gene}{$paa}!~/low risk/){
                                                                $ACMGrule[3]=1; $mark.=$hashAY{$gene}{$paa}.";";
                                                        }	 #ps3
                                                        elsif($hashAY{$gene}{$paa}=~/Neutral/){
                                                                $ACMGrule[19]=1;$mark.=$hashAY{$gene}{$paa}.";";
                                                        } #bs3
                                                }
                                        }
                                }
                                #pm1
                                my $pos=();  	
                                if(exists $pm1_domain{$gene}){
                                        if($en=~/exon/ && $pm1_domain{$gene}=~/$en;/){
                                                $ACMGrule[5]=1;
                                        }# pm1
                                        else{
                                                if($p =~ /p\.*(\d+).*/ ){
                                                        $pos=$1;
                                                        my @temp=split(/;/,$pm1_domain{$gene});
                                                        for(my $i=0;$i<=$#temp;$i++){
                                                                if($temp[$i]!~/-/ && $pos eq $temp[$i] ){
                                                                        $ACMGrule[5]=1;;
                                                                } #pm1
                                                                elsif($temp[$i]=~/(\d+)-(\d+)/ ){
                                                                        if($pos>=$1 && $pos<=$2){
                                                                                $ACMGrule[5]=1;
                                                                        }
                                                                }#pm1
                                                        }
                                                } 
                                        } 
                                        #   if($ACMGrule[5]){print $pm1_domain{$gene},"\t",$AAchange,"\n";}

                                }

                                last;	
                        }
                }
        }
        if ($ACMGrule[5]!=1 && $intp_dormain ne "."){   #yangrt 20180508  #pm1
                $ACMGrule[5]=1;
        }   
	if($AAchange ne "."){
                my @tmp=split(/,/,$AAchange,);
                for (my $i=0;$i<=$#tmp;$i++){
                        #BRCA1:NM_007298:exon1:c.G71A:p.C24Y	
                        my @exp=split(/:/,$tmp[$i],);
                        #print "$i $tmp[$i]\n";
                        for (my $l=0;$l<=$#exp;$l++){
                        	$NM=$exp[$l] if ($exp[$l]=~ /NM/);
                        	$en=$exp[$l] if ($exp[$l]=~ /exon/);
                        	$c=$exp[$l] if ($exp[$l]=~ /c\./);
                        	$p=$exp[$l] if ($exp[$l]=~ /p\./);
                        }
                        if(exists $hashnm{$gene}{$NM}){
                                $trans=$NM;
                                $exon=$en;
                                $hgvsc=$c;
                                $hgvsp=$p;	
                        };
                        #}
                        # print "NM $NM exon $en CDS $c AA $p\n";
                        if(exists $hashA{$string}){
                                #print "clinvar recode\n";
                                next;
                        }elsif(exists $hashc{$gene}{$NM}{$c} ){
                                #      print "exists hashc $tmp[$i] $gene $NM $c  ",$hashc{$gene}{$NM}{$c},"\n";
                                if ($mark !~ /clinvar/ ){
                                        $mark = "clinvar $hashc{$gene}{$NM}{$c};" . $mark;
                                }
                                if($hashc{$gene}{$NM}{$c} =~ /^\d+$/){
                                        #if($hashv{$hashc{$gene}{$NM}{$c}}->[0] =~ /Pathogenic/i){
                                        if($hashc_clinsig{$gene}{$NM}{$c} =~ /Pathogenic/i){        #yangrt   20170901
                                                # print "class $hashv{$hashc{$gene}{$NM}{$c}}->[0] status $hashv{$hashc{$gene}{$NM}{$c}}->[1] \n";
                                                #if($hashv{$hashc{$gene}{$NM}{$c}}->[1]=~/criteria provided, multiple submitters, no conflicts|criteria provided, single submitter|practice guideline|reviewed by expert panel/)
                                                #if($hashv{$hashc{$gene}{$NM}{$c}}->[1]=~/criteria provided, multiple submitters|criteria provided, single submitter|practice guideline|reviewed by expert panel/)    {
                                                if($hashc_revsta{$gene}{$NM}{$c} =~ /criteria provided, multiple submitters|criteria provided, single submitter|practice guideline|reviewed by expert panel/)    {      #yangrt    20170901
                                                        $ACMGrule[1]=1; #ps1
                                                        #$cita=$hashv{$hashc{$gene}{$NM}{$c}}->[1];
                                                        $cita=$hashc_revsta{$gene}{$NM}{$c};     #yangrt    20170210
                                                }
                                                #elsif($hashv{$hashc{$gene}{$NM}{$c}}->[1]=~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashv{$hashc{$gene}{$NM}{$c}}->[0] =~ /Pathogenic/)
                                                elsif($hashc_revsta{$gene}{$NM}{$c} =~ /no assertion criteria provided|no assertion for the individual variant|no assertion provided/ ||  $hashc_clinsig{$gene}{$NM}{$c} =~ /Pathogenic/){   #yangrt    20170901
                                                        $ACMGrule[15]=1;#pp5
                                                        #$cita=$hashv{$hashc{$gene}{$NM}{$c}}->[1];
                                                        $cita=$hashc_revsta{$gene}{$NM}{$c};        #yangrt   20170901
                                                }   
                                        }
                                        #elsif($hashv{$hashc{$gene}{$NM}{$c}}->[0] =~ /Benign/i){
                                        elsif($hashc_clinsig{$gene}{$NM}{$c} =~ /Benign/i){     #yangrt    20170901
                                                #if($hashv{$hashc{$gene}{$NM}{$c}}->[1]=~/criteria provided, multiple submitters, no conflicts|criteria provided, single submitter|practice guideline|reviewed by expert panel/){
                                                if($hashc_revsta{$gene}{$NM}{$c} =~ /criteria provided, multiple submitters, no conflicts|criteria provided, single submitter|practice guideline|reviewed by expert panel/){    #yangrt    20170901
                                                        $ACMGrule[18]=1;#bs2
                                                        #$cita=$hashv{$hashc{$gene}{$NM}{$c}}->[1];
                                                        $cita=$hashc_revsta{$gene}{$NM}{$c};        #yangrt   20170901
                                                }
                                                #elsif($hashv{$hashc{$gene}{$NM}{$c}}->[1]=~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashv{$hashc{$gene}{$NM}{$c}}->[0] =~ /Benign/)
                                                elsif($hashc_revsta{$gene}{$NM}{$c} =~ /no assertion criteria provided|no assertion for the individual variant|no assertion provided/  || $hashc_clinsig{$gene}{$NM}{$c} =~ /Benign/){    #yangrt   20170901
                                                        $ACMGrule[26]=1;#bp6
                                                        #$cita=$hashv{$hashc{$gene}{$NM}{$c}}->[1];
                                                        $cita=$hashc_revsta{$gene}{$NM}{$c};    #yangrt    20170901
                                                }
                                        }
                                }
                                else{
                                        $ACMGrule[15]=1;#pp5
                                }
                        }elsif(exists $hashp{$gene}{$NM}{$p}){
                                #print "exists hashp $tmp[$i] $gene $p ",$hashp{$gene}{$p},"\n";
                                
                                if($hashp{$gene}{$NM}{$p} =~ /^\d+$/){
                                        #print "class $hashv{$hashp{$gene}{$p}}->[0] status $hashv{$hashp{$gene}{$p}}->[1] \n";
                                        #if($hashv{$hashp{$gene}{$p}}->[0] =~ /Pathogenic/i){
                                        if($hashp_clinsig{$gene}{$NM}{$p} =~ /Pathogenic/i){        #yangrt  20170901
                                                $mark.="KPDB(clinvar $hashp{$gene}{$NM}{$p});";
                                                #if($hashv{$hashp{$gene}{$p}}->[1]=~/criteria provided, multiple submitters|criteria provided, single submitter|practice guideline|reviewed by expert panel/) {
                                                if($hashp_revsta{$gene}{$NM}{$p} =~/criteria provided, multiple submitters|criteria provided, single submitter|practice guideline|reviewed by expert panel/) {  #yangrt  20170901
                                                        $ACMGrule[1]=1; #ps1
                                                        $cita=$hashv{$hashp{$gene}{$NM}{$p}}->[1];
                                                        $cita=$hashp_revsta{$gene}{$NM}{$p};     #yangrt   20170901
                                                }
                                                #elsif($hashv{$hashp{$gene}{$p}}->[1]=~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashv{$hashp{$gene}{$p}}->[0] =~ /Pathogenic/)
                                                elsif($hashp_revsta{$gene}{$NM}{$p} =~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashp_clinsig{$gene}{$NM}{$p} =~ /Pathogenic/){   #yangrt  20170901
                                                        $ACMGrule[15]=1;#pp5
                                                        #$cita=$hashv{$hashp{$gene}{$p}}->[1];
                                                        $cita=$hashp_revsta{$gene}{$NM}{$p};          #yangrt   20170901
                                 
                                                }
                                        }
                                        #elsif($hashv{$hashp{$gene}{$p}}->[0] =~ /Benign/i){
                                        elsif($hashp_clinsig{$gene}{$NM}{$p} =~ /Benign/i){       #yangrt  20170901
                                                #if($hashv{$hashp{$gene}{$p}}->[1]=~/criteria provided, multiple submitters, no conflicts|criteria provided, single submitter|practice guideline|reviewed by expert panel/){
                                                if($hashp_revsta{$gene}{$NM}{$p} =~/criteria provided, multiple submitters, no conflicts|criteria provided, single submitter|practice guideline|reviewed by expert panel/){     #yangrt   20170901
                                                        $ACMGrule[18]=1;#bs2
                                                        #$cita=$hashv{$hashp{$gene}{$p}}->[1];
                                                        $cita=$hashp_revsta{$gene}{$NM}{$p};      #yangrt   20170901
                                                }
                                                #elsif($hashv{$hashp{$gene}{$p}}->[1]=~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashv{$hashp{$gene}{$p}}->[0] =~ /Benign/)
                                                elsif($hashp_revsta{$gene}{$NM}{$p} =~/no assertion criteria provided|no assertion for the individual variant|no assertion provided/ || $hashp_clinsig{$gene}{$NM}{$p} =~ /Benign/){    #yangrt   20170901
                                                        $ACMGrule[26]=1;#bp6
                                                        #$cita=$hashv{$hashp{$gene}{$p}}->[1];
                                                        $cita=$hashp_revsta{$gene}{$NM}{$p};
                                                }
                                        }
                                }
                                else{
                                        $ACMGrule[15]=1;#pp5
                                        
                                }
                        } 
                        elsif($p =~ /p\.(.)(\d+)(.*)/){ 
                                my $ref=$1; my $pos=$2; my $mut=$3; 
                                #print "$gene $ref $pos $mut\n";
				if(exists $hasha{$gene}{$NM}{$ref}{$pos}){
                                        # print "exists pm5 proof\n";
					if($mut ne $hasha{$gene}{$NM}{$ref}{$pos} && $exon_func eq "nonsynonymous SNV" && exists $hasht{$gene}{$NM}{$ref}{$pos}){
                                                #print "hasht $hasht{$gene}{$ref}{$pos} $hashr{$gene}{$ref}{$pos}\n";
                                                #if($hasht{$gene}{$ref}{$pos}=~/Pathogenic/i && $hasht{$gene}{$ref}{$pos}!~/Benign/i )
                                                if($hasht{$gene}{$NM}{$ref}{$pos}=~/Pathogenic/i && $hasht_CliSigSim{$gene}{$NM}{$ref}{$pos} ne "0"){
                                                        $mark.="KPDB_P(clinvar $hashtalleleid{$gene}{$NM}{$ref}{$pos}); ";
                                                        if ($hashtp{$gene}{$NM}{$ref}{$pos} eq "single nucleotide variant"){
                                                                $ACMGrule[9]=1;
                                                        }
                                                }#pm5
					}
				}
			}
                }
	}
        if($ACMGrule[1]){
                $ACMGrule[15]=0;
        }
#if($ACMGrule[16]){$ACMGrule[18]=0;}
# 0 pvs1 1 ps1 2 ps2 3 ps3 4 ps4 5 pm1 6 pm2 7 pm3 8 pm4 9 pm5 10 pm6 11 pp1 12 pp2 13 pp3 14 pp4 15 pp5
#16 ba1 17 bs1 18 bs2 19 bs3 20 bs4 21 bp1 22 bp2 23 bp3 24 bp4 25 bp5 26 bp6 27 bp7
#pm2 & ba1
#1000g2014oct_all        1000g2014oct_eas        esp6500siv2_all esp6500si_ea    esp6500si_aa    ExAC_ALL        ExAC_AFR        ExAC_AMR        ExAC_EAS        ExAC_FIN        ExAC_NFE    ExAC_OTH ExAC_SAS

#print "Popfreq 29 1000g2015aug_eas $atm[28] 34 esp6500siv2_ea $atm[33] 39 ExAC_EAS $atm[38]  44 PopFreqMax $atm[43] ";
        my $PopFreqMax=$atm[43];
        if($PopFreqMax eq "."){
                if($atm[28] ne "."){
                        $PopFreqMax=$atm[28];
                }
                if($atm[33] ne "."){
                        if($PopFreqMax eq "."){
                                $PopFreqMax=$atm[33];
                        }
                        elsif($PopFreqMax <$atm[33]){
                                $PopFreqMax=$atm[33];
                        }
                }
                if($atm[38] ne "."){
                        if($PopFreqMax eq "."){
                                $PopFreqMax=$atm[38];
                        }
                        elsif($PopFreqMax <$atm[38]){
                                $PopFreqMax=$atm[38];
                        }
                } 
        }
        else{
                if($atm[28] ne "." && $PopFreqMax <$atm[28]){$PopFreqMax=$atm[28];}
                if($atm[33] ne "." && $PopFreqMax <$atm[33]){$PopFreqMax=$atm[33];}
                if($atm[38] ne "." && $PopFreqMax <$atm[38]){$PopFreqMax=$atm[38];}
        }
        $atm[43]=$PopFreqMax;
        $mark.="Popfreq 1000g2015aug_eas $atm[28] / esp6500siv2_ea $atm[33] / ExAC_EAS $atm[38] / PopFreqMax $atm[43] ;";
        #	if($atm[43] eq "." && !exists($f400K{$string})){
        if($atm[43] eq "." || $atm[43]<=0.0001 ){      #yangrutao  20170926
		$ACMGrule[6]=1;#pm2
        }
        elsif($atm[43] ne "." && $atm[43] >0.002 && $atm[43] < 0.05){
                $ACMGrule[17]=1;
        }#bs1
        elsif($atm[43] ne "." && $atm[43] >= 0.05){
                $ACMGrule[16]=1; #ba1
        }
        if(exists($f400K{$string})){
                 $mark.="221e $f400K{$string};";
        #        print "400K $f400K{$string} ";
           #     if($atm[43] eq "."){
         #         print $string," higher frequency in 221e\n";
            #    }
             #  elsif($atm[43]<$f400K{$string}){
          #     print $string," higher frequency in 221e\n";
            #}
        }
        if(not $ACMGrule[6]){ #pm2
                if ($exon_func eq "synonymous SNV" && abs($dpsi_score)<2.0){   # WES  abs($spidex) >= 2.1   pp3
                        $ACMGrule[27]=1; #bp7
                }    
        }elsif($ACMGrule[6] && $ACMGrule[21]){
                if ($exon_func eq "synonymous SNV" && abs($dpsi_score)<2.0){   
                        $ACMGrule[27]=1; #bp7
                }
        }       
        if(exists($hashf{$string})){
                # print "8K $hashf{$string} ";
                $mark.="BZ8K $hashf{$string};";
                my  $freq=$hashf{$string};
                if($freq>0.05 && $freq<0.1){
                        $ACMGrule[17]=1;
                }#bs1
                elsif($freq>=0.1){
                        $ACMGrule[16]=1;
                }#ba1
                if($atm[43] ne "." ){
                        if($atm[43]>10*$hashf{$string}){
                                $ACMGrule[17]=1;
                        }#bs1
                        #elsif($hashf{$string}<0.05 && 10*$atm[43]<$hashf{$string}){$ACMGrule[4]=1;}#ps4
                }
        }
        if($ACMGrule[16]){
                $ACMGrule[18]=0;
        }#ba1 && bs1
# print "\n";
#pp3 & bp4
#score
	#gerp++gt2       SIFT_score      SIFT_pred       Polyphen2_HVAR_pred     LRT_pred        MutationTaster_pred     MutationAssessor_score  MutationAssessor_pred   FATHMM_score    FATHMM_pred     PROVEAN_score   PROVEAN_pred    VEST3_score     CADD_raw        CADD_phred      DANN_score      fathmm-MKL_coding_score fathmm-MKL_coding_pred  MetaSVM_score   MetaSVM_pred    MetaLR_score    MetaLR_pred     integrated_fitCons_score        integrated_confidence_value     GERP++_RS
	my ($sift,$polyphen,$LRT,$FATHMM,$MetaSVM,$MetaLR)=($atm[50],$atm[56],$atm[59],$atm[68],$atm[76],$atm[79]);
	my ($cadd,$gerp,$DANN)=($atm[45],$atm[47],$atm[46]);
	my ($Score_D,$Score_B)=(0,0);
	my $k=0;
	if($sift eq "D"){$k++;}elsif($sift eq "T"){$k--;}
	if($polyphen eq "D" || $polyphen eq "P"){$k++;}elsif($polyphen eq "B"){$k--;}
	if($LRT eq "D"){$k++;}elsif($LRT eq "N"){$k--;}
        if($FATHMM eq "D"){$k++;}elsif($FATHMM eq "T"){$k--;}
	if($MetaSVM eq "D"){$k++;}elsif($MetaSVM eq "T"){$k--;}
	if($MetaLR eq "D"){$k++;}elsif($MetaLR eq "T"){$k--;}
	if($cadd ne "." && $cadd >15){$k++;}
	if($gerp eq "."){$k--;}elsif($gerp >2){$k++;}
	if($DANN ne "." && $DANN >0.9){$k++;}
        #print "silico prediction value $k\n";
        $mark.="silico $k;";
	if($k >= 3){
		$ACMGrule[13]=1;#pp3
	}elsif($k <= -4){
		$ACMGrule[24]=1;#bp4
	}
# 0 pvs1 1 ps1 2 ps2 3 ps3 4 ps4 5 pm1 6 pm2 7 pm3 8 pm4 9 pm5 10 pm6 11 pp1 12 pp2 13 pp3 14 pp4 15 pp5
#16 ba1 17 bs1 18 bs2 19 bs3 20 bs4 21 bp1 22 bp2 23 bp3 24 bp4 25 bp5 26 bp6 27 bp7
 my $PVS=0;
 my $PS=0;
 my $PM=0;
 my $PP=0;
 my $BA=0;     
 my $BS=0;
 my $BP=0;
 my $evidence="";

if($ACMGrule[0]){$evidence="pvs1";$PVS=1;}
my ($i,$j);
for( $i=1;$i<=4;$i++)
{
if($ACMGrule[$i]){$evidence.="ps".$i; $PS+=1;}
}
for( $i=5;$i<=10;$i++)
{
if($ACMGrule[$i]){$j=$i-4;$evidence.="pm".$j; $PM+=1;}
}
for( $i=11;$i<=15;$i++)
{
if($ACMGrule[$i]){$j=$i-10;$evidence.="pp".$j; $PP+=1;}
}
if($ACMGrule[16]){$evidence.="ba1";$BA=1;}
for( $i=17;$i<=20;$i++)
{
if($ACMGrule[$i]){$j=$i-16;$evidence.="bs".$j; $BS+=1;}
}
for( $i=21;$i<=27;$i++)
{
if($ACMGrule[$i]){$j=$i-20;$evidence.="bp".$j; $BP+=1;}
}

print OUT "PVS",$PVS,"PS",$PS,"PM",$PM,"PP",$PP,"BA",$BA,"BS",$BS,"BP",$BP,"_",$evidence,"\t$_\t$trans\t$exon\t$hgvsc\t$hgvsp\t$mark\n";

}
close IN;
close OUT;


