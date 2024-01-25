#!/usr/bin/perl
use Data::Dumper;
my $file=shift;
my $outfile=shift;

open IN, "<$file";
open OU, ">$outfile";

my $first=<IN>;
my @array=split/\t/,$first;
my $chr=$array[1];
#my $start=$array[2];
my $end=$array[3];
my $out=join "\t",@array;
print OU "$out";
while(<IN>){
	chomp;
	my @unit=split;
	#print STDERR "$_\n";
	if($unit[1]==$chr && $unit[2]>$end){
		print OU "$_\n";
		$chr=$unit[1];
		#$start=$unit[2];
		$end=$unit[3];
	}elsif($unit[1]==$chr && $unit[2]<=$end && $unit[3]>$end){
		$start=$end+1;
		print OU "$unit[0]\t$chr\t$start\t$unit[3]\t$unit[4]\t$unit[5]\n";
		$chr=$unit[1];
		$end=$unit[3];
	}elsif($unit[1]==$chr && $unit[2]<=$end && $unit[3]<=$end){
		next;
		#print STDERR "$_\n";
	}elsif($unit[1]!=$chr){
		print OU "$_\n";
		$chr=$unit[1];
		#$start=$unit[2];
		$end=$unit[3];
	}else{
		print STEDRR "$_\n";
	}
	
}
close IN;
close OU;


