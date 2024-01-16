#! /usr/bin/perl -w
#perl overlap_var.pl [SNPs from unifiedgenotyper] [SNP from haplotypecaller]  [SNP from atlas2] > SNP output
#perl overlap_var.pl [SNPs from haplotypecaller] [SNP from atlas2]  [SNP from platypus]         > indel output
open IN,$ARGV[0]||die;
while(<IN>){
	chomp;
	if(/^#/){
	#	print "$_\n";
		next;
	}
	@tt=split/\s+/;
	if($tt[0]!~/chr/){
                $tt[0]="chr".$tt[0]
        }	
	my $key=$tt[0]."\t".$tt[1];
	$info1{$key}=$_;;
	$site{$key}++;
}
open IN1,$ARGV[1]||die;
while(<IN1>){
        chomp;
        if(/^#/){
                next;
        }
        @tt=split/\s+/;
	if($tt[0]!~/chr/){
                $tt[0]="chr".$tt[0]
        }
	my $key=$tt[0]."\t".$tt[1];
        $info2{$key}=$_;
	$site{$key}++;
}
open IN2,$ARGV[2]||die;
while(<IN2>){
        chomp;
        if(/^#/){
                next;
        }
        @tt=split/\s+/;
        if($tt[0]!~/chr/){
		$tt[0]="chr".$tt[0]
        }
	my $key=$tt[0]."\t".$tt[1];
	$info3{$key}=$_;
	$site{$key}++;
}
foreach my $i(keys %site){
	if(! defined $info1{$i} && ! defined $info2{$i}){
			print "$info3{$i}\n";
	}elsif(! defined $info1{$i} && ! defined $info3{$i}){
                print "$info2{$i}\n";
        }
	elsif(! defined $info2{$i} && ! defined  $info3{$i}){
		print "$info1{$i}\n";
	}
}#print "$m\t$n\t$s\n";
