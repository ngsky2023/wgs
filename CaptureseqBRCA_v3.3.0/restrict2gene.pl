#! /usr/bin/perl -w
#perl restrict2bed.pl [.gff] [dbsnp141] > [restricted_dbsnp]
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	@tt=split/\s+/;
	push @{$info{$tt[0]}},[$tt[3],$tt[4]];
}
my %out;
open IN1,$ARGV[1]||die;
while(<IN1>){
	chomp;
	if(/^#/){print "$_\n";next};
	$line=$_;
	if($line!~/^chr/){
        	$line="chr".$line;
        }
        @tt=split(/\s+/,$line);
	$key=$tt[0];
	$kk=$tt[0]."\t".$tt[1];	
	if(defined $info{$key}){
		@ss=@{$info{$key}};
	}else{
		next;
	}
	foreach my $i(0..$#ss){
		if($tt[1]>=$ss[$i][0] && $tt[1]<=$ss[$i][1] && (! defined $out{$kk})){
			#print "$line\n";
			print "$line\n";
			$out{$kk}++;
		}
	}
}
