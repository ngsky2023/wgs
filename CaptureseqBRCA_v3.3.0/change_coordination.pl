#! /usr/bin/perl -w
print "\@HD\tVN:1.0\tSO:coordinate\n";
print "\@SQ\tSN:chr1\tLN:249250621\n";
print "\@SQ\tSN:chr2\tLN:243199373\n";
print "\@SQ\tSN:chr3\tLN:198022430\n";
print "\@SQ\tSN:chr4\tLN:191154276\n";
print "\@SQ\tSN:chr5\tLN:180915260\n";
print "\@SQ\tSN:chr6\tLN:171115067\n";
print "\@SQ\tSN:chr7\tLN:159138663\n";
print "\@SQ\tSN:chr8\tLN:146364022\n";
print "\@SQ\tSN:chr9\tLN:141213431\n";
print "\@SQ\tSN:chr10\tLN:135534747\n";
print "\@SQ\tSN:chr11\tLN:135006516\n";
print "\@SQ\tSN:chr12\tLN:133851895\n";
print "\@SQ\tSN:chr13\tLN:115169878\n";
print "\@SQ\tSN:chr14\tLN:107349540\n";
print "\@SQ\tSN:chr15\tLN:102531392\n";
print "\@SQ\tSN:chr16\tLN:90354753\n";
print "\@SQ\tSN:chr17\tLN:81195210\n";
print "\@SQ\tSN:chr18\tLN:78077248\n";
print "\@SQ\tSN:chr19\tLN:59128983\n";
#print "\@SQ\tSN:chr20\tLN:63025520\n";
#print "\@SQ\tSN:chr21\tLN:48129895\n";
print "\@SQ\tSN:chr22\tLN:51304566\n";
#print "\@SQ\tSN:chrX\tLN:155270560\n";
#print "\@SQ\tSN:chrY\tLN:59373566\n";
print "\@RG\tID:$ARGV[0]\tPL:illumina\tSM:$ARGV[0]\n";
open IN1,$ARGV[2]||die;
while(<IN1>){
	chomp;
	$info{$_}++;
}
open IN,$ARGV[1]||die;
while(<IN>){
	chomp;
	next if(/^\@/);
	@tt=split/\s+/;
	next if(defined $info{$tt[9]});
	if($tt[2]=~/(chr.*)_(\d+)_.*/){
		$tt[3]=$tt[3]+$2;
		$tt[2]=$1;
		if($tt[6] eq "\="){
			$tt[7]=$tt[7]+$2;
		}elsif($tt[6]=~/(chr.*)_(\d+)_.*/){
			$tt[7]=$tt[7]+$2;
			$tt[6]=$1;
			next if($tt[6] ne $tt[2]);
		}
	}
	foreach my $i (0..$#tt){
		if($tt[$i]=~/^(SA:Z:chr.*)_(\d+)_(.*)\,(\d+)(\,[\+\-].*)/){
			my $coo=$2+$4;
			$tt[$i]=$1.",".$coo.$5;
		}
	}
	$out=join("\t",@tt);
	print "$out\n";
}
