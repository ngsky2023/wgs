#! /usr/bin/perl -w
my ($c,$j)=(0,0);
open IN,"gzip -dc $ARGV[0]|"||die;
while(<IN>){
	chomp;
	$c++;
}
close IN;
$c=$c/4;
if($c<$ARGV[2]){exit}
my %data;
#srand();
my $i=0;
while($j<$ARGV[2]){
	srand(time()+$i);
	$num=int(rand($c));
	$i++;
	next if(defined $data{$num});
	$data{$num}++;
	$j=(keys %data);		
}
$t=0;
my $out1=$ARGV[0];
$out1=~s/\.gz//g;
open OUT," | gzip -f - > $out1.$ARGV[2].gz"||die;
open IN,"gzip -dc $ARGV[0] |"||die;
while(<IN>){
	$t++;	
	$id=$_;
        $seq=<IN>;
        <IN>;
        $qual=<IN>;#print "$seq\n";
	if(defined $data{$t}){
        	print OUT "$id$seq+\n$qual";
	}
}
$t=0;
my $out2=$ARGV[1];
$out2=~s/\.gz//g;
open IN1,"gzip -dc $ARGV[1]|"||die;
open OUT1,"| gzip -f - >$out2.$ARGV[2].gz"||die;
while(<IN1>){
	$t++;
	$id=$_;
	$seq=<IN1>;
	<IN1>;
	$qual=<IN1>;
	if(defined $data{$t}){
		print OUT1 "$id$seq+\n$qual";
        }
}
