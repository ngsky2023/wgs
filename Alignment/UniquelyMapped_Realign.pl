#!/usr/bin/env perl
use strict;
use warnings;

die "usage: perl $0 <Input.bam> > <Output.sam>\n" unless ($#ARGV>=0);
my $fr=shift;
open(IN,"samtools view -h -q 10 -F 4 $fr|")||die "Cannot open file $fr\n";
while(my $a=<IN>){
  if ($a =~ /^\@/){
	print $a;
	next;
  }
  chomp $a;
  my @as=();
  my @xs=();
  my $i=0;
  my @b=split(/\t/,$a);
	while($i<@b){
		if($b[$i]=~/^AS:/){
			@as=split(/:/,$b[$i]);
		}	# AS:i:[Score] alignment score
		if($b[$i]=~/^XS:/){
			@xs=split(/:/,$b[$i]);
		}	# XS:i:[Score] suboptimal alignment score
		$i++;
	}
  my $diff=$as[2]-$xs[2];
  if($diff>10){
    print  "$a\n";
  }
}
close(IN);

print STDERR "$fr Done!\n";
