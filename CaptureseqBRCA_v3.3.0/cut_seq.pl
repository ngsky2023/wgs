#! /usr/bin/perl -w
open IN,"gzip -dc $ARGV[0]|"||die;
if($ARGV[3]==0){
	open OUT,"| gzip -f - > $ARGV[2]"||die;
}else{
	open OUT,">$ARGV[2]"||die;
}
while(<IN>){
	$id=$_;
	$seq=<IN>;
	<IN>;
	$qual=<IN>;
	$seq1=substr($seq,0,$ARGV[1]);
	$qual1=substr($qual,0,$ARGV[1]);
	print OUT "$id$seq1\n+\n$qual1\n";
}
