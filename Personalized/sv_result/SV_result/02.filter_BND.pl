use strict;
use warnings;
use Cwd;
use Data::Dumper;

my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
my %hash;
foreach my $sample(@sample){
	open F,"../SV_result_raw/$sample.somaticSV.xls";
	my %id;
	while(<F>){
		chomp;
		next unless /^chr/;
		my @line=split/\t/,$_;
		next if $line[10] eq 'INS';
		#if ($_=~/SVLEN=(\d+)/){
		#	next if $1<500;
		#	# print "$1\n";
		#}
		if($line[10]){
			if($line[10] eq 'BND'){

				my @id=split/:/,$line[6];
				my $id=join "-",@id[0..6];
				if(defined $id{$id}){
					next;
				}else{
					$id{$id}=1;
				}
			}
		}else{
			print "$sample\t$.\n";
		}
		my $chr=join "\t",@line[0..5];
		if(defined $hash{$chr}){
			$hash{$chr}+=1;
		}else{
			$hash{$chr}=1;
		}
	}
	close F;
}
my %type;
foreach my $sample(@sample){
	open F,"../SV_result_raw/$sample.somaticSV.xls";
	open OUT,">$sample.somaticSV.xls";
    my %id;
	while(<F>){
		chomp;
		if (/^#/){
			print OUT "$_\n";
			next;
		}
		next unless /^chr/;
		my @line=split/\t/,$_;
		next if $line[10] eq 'INS';
		#if (/SVLEN=(\d+)/){
		#	next if $1<500;
		#}
		if($line[10] eq "BND"){
			my @id=split/\:/,$line[6];
			my $id=join ":",@id[0..6];
			if(defined $id{$id}){
				next;
			}else{
					$id{$id}=1;
			}
		}
		my $chr=join "\t",@line[0..5];
		#if($hash{$chr} >=5){
		#	next;
		#}
		print OUT "$_\n";
		if(defined $type{$sample}{$line[10]}){
			$type{$sample}{$line[10]}+=1;
		}else{
			$type{$sample}{$line[10]}=1;
		}
	}
		
}
close OUT;

my @class1=("BND","DEL","DUP","INV");
open OUT,">sv_filter_BND_number.xls";
my $title1=join "\t",@class1;
print OUT "sample\t$title1\n";
foreach my $sample(@sample){
        print OUT "$sample\t";
        foreach my $type(@class1){
			$type{$sample}{$type}||=0;
			print OUT "$type{$sample}{$type}\t"
		}
		print OUT "\n";
}




