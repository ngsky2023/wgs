use strict;
use warnings;

die "perl $0 <in> <out>" unless @ARGV==2;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";

my $title=<IN>;
print OUT "Level\t$title";
while(my $line=<IN>){
	chomp($line);
	my @atm=split (/\t/,$line,);
	$atm[0]=~/PVS(\d+)PS(\d+)PM(\d+)PP(\d+)BA(\d+)BS(\d+)BP(\d+)_/;
 my $PVS=$1;
 my $PS=$2;
 my $PM=$3;
 my $PP=$4;
 my $BA=$5;
 my $BS=$6;
 my $BP=$7;
 my $P="0";
 my $B="0";
 if($PVS)
 {
 if($PS>=1){$P="D";}	
 elsif($PM>=2){$P="D";}	
 elsif($PM>=1 && $PP>=1){$P="D";}
 elsif($PP>=2){$P="D";}
 elsif($PM==1){$P="LD";}
 }
 elsif($PS>=2){$P="D";}
 elsif($PS==1)
 {
 	if($PM>=3){$P="D";}
 	elsif($PM==2 && $PP>=2){$P="D";}
 	elsif($PM==1 && $PP>=4){$P="D";}
 	elsif($PM>0 && $PM<=2){$P="LD";}
 	elsif($PP>=2){$P="LD";}
 }
 elsif($PM>=3){$P="LD";}
 elsif(($PM==2 && $PP>=2) || ($PM==1 && $PP>=4)){$P="LD";}
 
 if($BA){$B="B";}
 elsif($BS>=2){$B="B";}
 elsif($BS==1 && $BP==1){$B="LB";}
 elsif($BP>=2){$B="LB";}

my $level;
if($P eq "D" && $B eq "0"){$level="Pathogenic";}
elsif($P eq "LD" && $B eq "0"){$level="Likely pathogenic";} 
elsif($P eq "0" && $B eq "B"){$level="Benign";}
elsif($P eq "0" && $B eq "LB"){$level="Likely benign";}
elsif(($P=~/D/ && $B=~/B/)||($P eq "0" && $B eq "0")){$level="Uncertain significance";}
else{print $line,"\n";}
print OUT "$level\t$line\n";
}

close(IN);
close(OUT);
