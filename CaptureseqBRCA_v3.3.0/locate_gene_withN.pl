#! /usr/bin/perl -w
#perl locate.pl hg19.fa.fai zz hg19.fa >h19.panel.fa
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	my @tt=split/\s+/;
	$len{$tt[0]}=$tt[1];
}
open IN1,$ARGV[1]||die;
while(<IN1>){
	chomp;
	my @tt=split/\s+/;
	push @{$ss{$tt[0]}},[$tt[3],$tt[4]];
}
$/=">";
open IN2,$ARGV[2]||die;
while(<IN2>){
	chomp;
	next if($_ eq "");
	my ($chr,$seq)=split(/\s+/,$_,2);
	next if(! defined $ss{$chr});
	$seq=~s/\n//g;
	$seq=~s/\s+//g;
	$fa{$chr}=$seq;
	my @mm=@{$ss{$chr}};
	@mm= sort {$a->[0] <=> $b->[0]} @mm;
	my ($s,$e,$new)=(1,1,"");
	foreach my $i(0..$#mm){
		$e=$mm[$i][0]-1;
		$len1=$e-$s+1;
		$ss= "N" x $len1;
		$s=$mm[$i][1]+1;
		$str=substr($fa{$chr},$e,($s-$e-1));
		$new.=$ss;
		$new.=$str;
	}
	$len1=$len{$chr}-$s+1;
	$ss= "N" x $len1;
	$new.=$ss;
	print ">$chr\n$new\n";
		
}
