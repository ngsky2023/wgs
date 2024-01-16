#! /usr/bin/perl -w
#perl locate.pl hg19.fa.fai zz hg19.fa >h19.panel.fa
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	my @tt=split/\s+/;
	$tt[8]=~s/\;.*//g;
        $tt[8]=~s/ID\=//g;
        $s=$tt[3]-1;
	push @{$info{$tt[0]}},[$tt[3],$tt[4],"$tt[0]\_$s\_$tt[8]"];
}
open IN1,$ARGV[1]||die;
while(<IN1>){
	chomp;
	my $line=$_;
	if(/^#/){
		print "$line\n";
	}else{
	my @tt=split/\s+/;
	@ss=@{$info{$tt[0]}};
        $kk=$tt[0]."\t".$tt[1];
	foreach my $i(0..$#ss){
                if($ss[$i][0]<=$tt[1] && $tt[1]<=$ss[$i][1] && (! defined $out{$kk})){
                        $tt[1]=$tt[1]-$ss[$i][0]+1;
			$tt[0]=$ss[$i][2];
			$ou=join("\t",@tt);
			print "$ou\n";
                        $out{$kk}++;
                }
        }
}}
