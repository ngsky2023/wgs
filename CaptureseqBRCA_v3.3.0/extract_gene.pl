#! /usr/bin/perl -w
#perl locate.pl hg19.fa.fai zz hg19.fa >h19.panel.fa
$/=">";
open IN2,$ARGV[1]||die;
while(<IN2>){
	chomp;
	next if($_ eq "");
	my ($chr,$seq)=split(/\s+/,$_,2);
	$seq=~s/\n//g;
	$seq=~s/\s+//g;
	$fa{$chr}=$seq;
}
$/="\n";
open IN1,$ARGV[0]||die;
while(<IN1>){
        chomp;
        my @tt=split/\s+/;
        $tt[8]=~s/\;.*//g;
        $tt[8]=~s/ID\=//g;
	$s=$tt[3]-1;
	$str=substr($fa{$tt[0]},$s,($tt[4]-$tt[3]+1));
	print ">$tt[0]\_$s\_$tt[8]\n$str\n";
		
}
