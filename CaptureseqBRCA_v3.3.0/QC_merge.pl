#! /usr/bin/perl -w
## input is   [paste   -d "\t" */6.QC/QC/*_QC.xls  >zz]
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	my @tt=split(/\t/,$_);
	$index =((scalar @tt)/2)-1;	
	push @tit,$tt[0];
	foreach my $i(0..$#tt){
		if($i % 2 ==1){
			push @{$info{$tt[0]}},$tt[$i];	
			
		}
	}
}
$out1=join("\t",@tit);
print "$out1\n";
foreach my $j (0..$index){
	my @out;
	foreach(@tit){
		push @out,$info{$_}[$j];
	}	
	$transform=join("\t",(@out));
	print "$transform\n";
}
