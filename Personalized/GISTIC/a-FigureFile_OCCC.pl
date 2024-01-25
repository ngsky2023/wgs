#!usr/bin/perl
open(FR,"Markers_ZhongHuaAll_Index.txt")||die "$!";
while(defined($a=<FR>)){
  chomp $a;
  @b=split(/\t/,$a);
  $loc="$b[1]\t$b[2]";
  $c{$loc}=$b[3];
}
close(FR);
$a="";@b=();$loc="";
##############################################################################
open(FR,"result/scores.gistic")||die "$!"; ## 用前先去掉单位数染色体（chr 1~9）编号前的空格
while(defined($a=<FR>)){
  chomp $a;
  @b=split(/\t/,$a);
  $start=$b[2];
  $b[1]=~s/\s//g;
LOCI:while($start<=$b[3]){
    $loc="$b[1]\t$start";
    if(not $c{$loc}){
      $start++;
      next LOCI;
    }
    if($c{$loc} and ($b[0] eq "Amp")){
      $amp{$loc}="$b[4]\t$b[5]\t$b[6]\t$b[7]";
      if($b[7]<=0.02){
        $colamp{$loc}="FB6A4A";
      }elsif($b[7]>0.02 and $b[7]<=0.05){
        $colamp{$loc}="EF3B2C";
      }elsif($b[7]>0.05 and $b[7]<=0.07){
        $colamp{$loc}="CB181D";
      }elsif($b[7]>0.07 and $b[7]<=0.1){
        $colamp{$loc}="A50F15";
      }elsif($b[7]>0.1){
        $colamp{$loc}="67000D";
      }
    }elsif($c{$loc} and ($b[0] eq "Del")){
      $del{$loc}="$b[4]\t-$b[5]\t-$b[6]\t-$b[7]";
      if($b[7]<=0.02){
        $coldel{$loc}="6BAED6";
      }elsif($b[7]>0.02 and $b[7]<=0.05){
        $coldel{$loc}="4292C6";
      }elsif($b[7]>0.05 and $b[7]<=0.07){
        $coldel{$loc}="2171B5";
      }elsif($b[7]>0.07 and $b[7]<=0.1){
        $coldel{$loc}="08519C";
      }elsif($b[7]>0.1){
        $coldel{$loc}="08306B";
      }
    }
    $start++;
  }
}
close(FR);
$a="";@b=();%c=();$start=0;
##############################################################################
open(FR,"Markers_ZhongHuaAll_Index.txt")||die "$!";
open(FW,">result/wgs_Chart.txt")||die "$!";
while(defined($a=<FR>)){
  chomp $a;
  @b=split(/\t/,$a);
  $loc="$b[1]\t$b[2]";
  if($b[0] eq "Name"){
    print FW "$a\tLogQ_Amp\tG_Amp\tLRR_Amp\tFreq_Amp\tColFreq_Amp\tLogQ_Del\tG_Del\tLRR_Del\tFreq_Del\tColFreq_Del\n";
  }elsif($amp{$loc} and $del{$loc}){
    print FW "$a\t$amp{$loc}\t$colamp{$loc}\t$del{$loc}\t$coldel{$loc}\n";
  }
}
close(FR);close(FW);
$a="";@b=();%amp=%del=();%colamp=%coldel=();
##############################################################################
#system("/share/work2/yangar/plot/wgs/GISTIC/Chart.R");
