use strict;
use warnings;
use Cwd;
use Data::Dumper;
my %hash;
my %count;
my %length;
my $path="./result";
my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);

chomp(my $head =readpipe("cat head.txt"));
open S,">ac_kata_number.xls";
print S "$head\n";
my %base;
my @base=("C>G","C>T","C>A","T>C","T>G","T>A");
#my @base=("C>G\tC>T\tC>A\tT>C\tT>G\tT>A)
my %count2;
foreach my $sample(@sample){
	my $kata_result="$path/${sample}_Kataegis.tsv";
	if(-e $kata_result){
		# chomp(my $lane=`wc -l $kata_result|cut -d" " -f1`);
		# my $col_num=$lane-1;
		# print S "$sample\t$col_num\n";
		my @type;
		open F,$kata_result;
		while(<F>){
			chomp;
			my @line=split/\t/,$_;
			if (/^Chromosome/){
				@type=@line[8..$#line-1];
				next;
			}
			if (defined $count2{$sample}){
				$count2{$sample}+=1;
			}else{
				$count2{$sample}=1;
			}
			my $out=join "\t",@line[0..6];
			print S "$out\n";
			if(defined $hash{$sample}{$line[0]}){
				$hash{$sample}{$line[0]}+=1;
				$count{$sample}{$line[0]}+=$line[3];
				$length{$sample}{$line[0]}+=$line[4];
				
			}else{
				$hash{$sample}{$line[0]}=1;
				$count{$sample}{$line[0]}=$line[3];
				$length{$sample}{$line[0]}=$line[4];
			}
			my $num=8;
			foreach my $type(@type){
				$num+=1;
				my $base_num=$line[$num];
				next if $base_num eq 'NA';
				# print "$type\t$base_num\n";
				if(defined $base{$sample}{$type}){
					$base{$sample}{$type}+=$base_num;
				}else{
					$base{$sample}{$type}=$base_num;
				}
				
			}
		}
	}else{
		print S "0\t0\t0\t0\t0\t0\t$sample\n";
	}
	
}
close S;
#
open CN,">ac_kata_sample_seq.xls";
print CN "sample\tseq\n";
open OUT,">ac_kata_chr_count.xls";
open CH,">ac_kata_chr.xls";
print CH "sample\tchr\tmutation\n";
print OUT "Tumor_Sample_Barcode\tChromosome\tChrom_count\tnMuts\tAvg_intermutation_dist\n";
foreach my $sample(@sample){
	$count2{$sample}||=0;
	print CN "$sample\t$count2{$sample}\n";
	if (exists $hash{$sample}){
		my $hash2=$hash{$sample};
		# print Dumper \%$hash2;
		foreach my $chr(keys %$hash2){
			print CH "$sample\t$chr\tKATA\n";
			print OUT "$sample\t$chr\t$hash{$sample}{$chr}\t$count{$sample}{$chr}\t$length{$sample}{$chr}\n";
		}
	}else{
		print OUT "$sample\t0\t0\t0\t0\n";
	}
}
close OUT;
#碱基
open BA,">ac_kata_base_num.xls";
my $title =join "\t",@base;
print BA "sample\t$title\n";
foreach my $sample(@sample){
	print BA "$sample\t";
	foreach my $type(@base){
		$base{$sample}{$type}||=0;
		print BA "$base{$sample}{$type}\t";
	}
	print BA "\n";
}







