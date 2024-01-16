#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
die "usage: perl $0 <BED file of captured targets with gene annotation> <LRR Input> <Output: LRR gene summary> <Output: CNV gene summary>\n" if $#ARGV<4;

my $fr1=shift;
my $fr2=shift;
my $fw1=shift;
my $fw2=shift;
my $tp = shift;


my %gene;
open(ANN,"<$fr1")||die "$!";
while(my $ann=<ANN>){
  chomp $ann;
  my @ann=split(/\t/,$ann);
  next if $ann[0] eq "Chr";
  my $start=$ann[1]+1;
  my $pos="$ann[0]\:$start\-$ann[2]";
  $gene{$pos}="$ann[3]";
}
close(ANN);

my ($cut_gain , $cut_loss);
if ($tp eq 'P'){
	$cut_gain=0.35;$cut_loss=-0.5;
}else{
	$cut_gain=0.35;$cut_loss=-0.5;
}
my %samp;
my $sampnum;
my %count;
my %chr;
my %start;
my %end;
my %sumlrr;
my %gain;
my %loss;
open(LRR,"<$fr2")||die "$!";
open(FW,">$fw1")||die "$!";
while(my $lrr=<LRR>){
  chomp $lrr;
  my @lrr=split(/\t/,$lrr);
  my @pos=split(/[:-]/,$lrr[0]);
  my $col=1;
  if($lrr[0] eq "Target"){
    print FW "Gene\tChr\tStart\tEnd";
    while($col<@lrr){
      print FW "\t$lrr[$col]";
      $samp{$col}=$lrr[$col];
      $sampnum++;
      $col++;
    }
    print FW "\n";
  }else{
    next if not exists $gene{$lrr[0]};
    next if $lrr[1] eq "NaN";
	next if $lrr[1] eq 'Inf' or $lrr[1] eq '-Inf';
    my $gene=$gene{$lrr[0]};
    $count{$gene}++;
    $chr{$gene}=$pos[0];
    if($count{$gene}==1){
      $start{$gene}=$pos[1];
    }else{
      $end{$gene}=$pos[1];
    }
    while($col<@lrr){
      my $genesamp="$gene\t$col";
	  
      $sumlrr{$genesamp}+=$lrr[$col];
      if($lrr[$col]>$cut_gain){
        $gain{$genesamp}++;
      }elsif($lrr[$col]<$cut_loss){
        $loss{$genesamp}++;
      }
      $col++;
    }
  }
}
close(LRR);


open(CN,">$fw2")||die "$!";
print CN "Gene\tChr\tStart\tEnd";
my $col=1;
while($col<=($sampnum)){
  print CN "\t$samp{$col}";
  $col++;
}
print CN "\n";

my @gene=sort {$a cmp $b} keys %count;

my $i=0;
my $countCut = 3;
while($i<@gene){
  my $gene = $gene[$i];

  my $count = $count{$gene};
  if($count<=$countCut or $gene =~ /^\d+$/){
    $i++;
    next;
  }
  print FW "$gene\t$chr{$gene}\t$start{$gene}\t$end{$gene}";
  print CN "$gene\t$chr{$gene}\t$start{$gene}\t$end{$gene}";
  $col=1;
  while($col<=($sampnum)){
    my $genesamp="$gene\t$col";
    my $avglrr=$sumlrr{$genesamp}/$count{$gene};
	my $copyNum = sprintf("%.2f" , 2 ** ($avglrr+1));
    if($gain{$genesamp} and $loss{$genesamp}){
      print FW "\t$avglrr|Num=$count{$gene}|Gain=$gain{$genesamp}|Loss=$loss{$genesamp}";
    }elsif($gain{$genesamp} and not $loss{$genesamp}){
      print FW "\t$avglrr|Num=$count{$gene}|Gain=$gain{$genesamp}|Loss=0";
    }elsif(not $gain{$genesamp} and $loss{$genesamp}){
      print FW "\t$avglrr|Num=$count{$gene}|Gain=0|Loss=$loss{$genesamp}";
    }elsif(not $gain{$genesamp} and not $loss{$genesamp}){
      print FW "\t$avglrr|Num=$count{$gene}|Gain=0|Loss=0";
    }
	my $cnv;
	$loss{$genesamp} = 0 unless exists $loss{$genesamp};
	$gain{$genesamp} = 0 unless exists $gain{$genesamp};
    if((($gain{$genesamp}-$loss{$genesamp})/$count{$gene})>=0.7 and $count{$gene}>=$countCut+1){
      $cnv="1|$copyNum";
    }elsif((($gain{$genesamp}-$loss{$genesamp})/$count{$gene})<=-0.7 and $count{$gene}>=$countCut+1){
      $cnv="-1|$copyNum";
    }else{
      $cnv=0;
    }
    print CN "\t$cnv";
    $col++;
  }
  print FW "\n";
  print CN "\n";
  $i++;
}
close(FW);close(CN);
