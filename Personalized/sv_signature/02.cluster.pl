use strict;
use warnings;
use Cwd;
use Data::Dumper;
my $samplelist="../clinical/sample.list";
my $out_path="./bedtools";
open IN,$samplelist;
chomp(my @sample=<IN>);
close IN;
my @clust=("clustered","non-clustered");
my @type=("del","dup","inv","tra");#DEL DUP INV TRA
my @class2=("1-10kb","10-100kb","100kb-1Mb","1Mb-10Mb",">10Mb");
open OUT,">class.xls";
print OUT "type\t";
my %type;
foreach my $sample(@sample){
	system("sort -k1,1V -k2,2n  $out_path/$sample.xls> $out_path/$sample.sort.xls");
	system("/share/public/software/bedtools2/2.27/bin/bedtools cluster -i $out_path/$sample.sort.xls -d 1000000 >$out_path/$sample.bedtools.xls");
	my %hash;
	open M,"$out_path/$sample.bedtools.xls";
	while(<M>){
		chomp;
		next if /^\s+$/;
		my $flag=(split/\t/,$_)[-1];
		#print "$flag\t$_\n";
		if(defined $hash{$flag}){
			$hash{$flag}+=1;
		}else{
			$hash{$flag}=1;
		}
	}
	close M;
	open G,"$out_path/$sample.bedtools.xls";
	while(<G>){
		chomp;
		next if /^\s+$/;
		my @line=split/\t/,$_;
		my $class1;
		if(exists $hash{$line[-1]} && $hash{$line[-1]} >=10){
			$class1="clustered";
		}else{
			$class1="non-clustered";
		}
		#chr1    960303  1960303 TRA     TRA     1
		#chr1    1204260 1246577 INV     42316   1
		my ($class2);
		$line[3]=lc $line[3];
		if($line[3] eq "tra"){
			$class2="T";
		}else{
			my $length=$line[4]/1000;
			if($length<=10){
			$class2="1-10kb";
			}elsif($length<=100){
				$class2="10-100kb";
			}elsif($length<=1000){
				$class2="100kb-1Mb";
			}elsif($length <=100000){
				$class2="1Mb-10Mb";
			}elsif($length >10000){
				$class2=">10Mb";
			}
		}
		my $class="$class1\[$line[3]\]$class2";#class type
		if(defined $type{$sample}{$class}){
				$type{$sample}{$class}+=1;
		}else{
				$type{$sample}{$class}=1;
		}

        }
	print OUT "$sample\t";
}
#system("paste clonames.xls *.class.xls >class.xls");
print OUT "\n";
my @uniq_class1=("clustered","non-clustered");
my @uniq_class2=("del","dup","inv","tra");
my @uniq_class3;
foreach my $class1(@uniq_class1){
        foreach my $class2(@uniq_class2){
			if($class2 eq 'tra'){
				@uniq_class3="T";
			}else{
				@uniq_class3=("1-10kb","10-100kb","100kb-1Mb","1Mb-10Mb",">10Mb");
			}
			foreach my $class3(@uniq_class3){
				my $key1="$class1\[$class2\]$class3";#class type
				print "$key1\n";
				print OUT "$key1\t";
				foreach my $sample(@sample){
						$type{$sample}{$key1} ||=0;
						print OUT "$type{$sample}{$key1}\t";
				}
				print OUT "\n";
			}

        }
}
