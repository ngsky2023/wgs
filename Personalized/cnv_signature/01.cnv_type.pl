use strict;
use warnings;
use Data::Dumper;
open S ,"../clinical/sample.list";
chomp(my @sample=<S>);
close S;
chomp(my $head=readpipe("cat head"));
open OUT,">cnv_class.xls";
print OUT "type\t"; ##biaoti
my (@class1,@class2,@class3,%hash);
foreach my $sample(@sample){
	my $cnv="../cnv_result/$sample.purple.cnv.somatic.tsv";
	open IN,$cnv;
	while(<IN>){
		my ($class,$class1,$class2,$class3);
		chomp;
		next if /^chromosome/;
		my @line=split/\t/,$_;
		next if $line[3]<0;
		my $length=($line[2]-$line[1])/1000;
		my $copyNumber=$line[3];
		my $minorAllelePloidy=sprintf("%.0f",$line[14]);
		my $majorAllelePloidy=sprintf("%.0f",$line[15]);
		if($minorAllelePloidy ==0 && $majorAllelePloidy !=0){
				$class1='LOH';
				if($copyNumber <=1){
					$class2='del'
				}elsif($copyNumber>1 && $copyNumber <=3){
					$class2='neut';
				}elsif($copyNumber>3 && $copyNumber <=4){
					$class2='dup';
				}elsif($copyNumber >=4){
					$class2='amp';
				}else{
					$class2="nrrow\-$copyNumber";
				}
		}elsif($minorAllelePloidy ==0 && $majorAllelePloidy ==0){
				$class1='homdel';
				$class2='del';
				
		}else{
				$class1='het';
				if($copyNumber <=3){
					$class2='neut';
				}elsif($copyNumber>3 && $copyNumber <=4){
					$class2='dup';
				}elsif($copyNumber >=4){
					$class2='amp';
				}else{
					$class2="nrrow\-$copyNumber";
				}
		}
		# if($length<=10){
			# $class3="0-0.01Mb";
		# }elsif($length<=100){
				# $class3="0.01-0.1Mb";
		# }elsif($length<=1000){
				# $class3="0.1-1Mb";
		# }elsif($length <=100000){
				# $class3="1-10Mb";
		# }elsif($length >10000){
				# $class3=">10Mb";
		# }
		if($length<=10){
			$class3="\(-0.01,0.01\]";
		}elsif($length<=100){
				$class3="\(0.01,0.1\]";
		}elsif($length<=1000){
				$class3="\(0.1,1\]";
		}elsif($length <=100000){
				$class3="\(1,10\]";
		}elsif($length >10000){
				$class3="\(10,Inf\]";
		}
		$class="$class2\:$class1\:$class3";#class type
		if(defined $hash{$sample}{$class}){
			$hash{$sample}{$class}+=1;
		}else{
			$hash{$sample}{$class}=1;
		}
		
	}
	print OUT "$sample\t";
}
print OUT "\n";
my @uniq_class1=("het","LOH","homdel");
my @uniq_class2;
my  @uniq_class3=('(-0.01,0.01]','(0.01,0.1]','(0.1,1]','(1,10]','(10,Inf]');
foreach my $class1(@uniq_class1){
	if($class1 eq 'het'){
		@uniq_class2=("amp","dup","neut");
	}elsif($class1 eq 'LOH'){
		@uniq_class2=("amp","del","dup","neut");
	}else{
		@uniq_class2=("del");
	}
	foreach my $class2(@uniq_class2){
		foreach my $class3(@uniq_class3){
			my $key1="$class2\:$class1\:$class3";#class type
			print OUT "$key1\t";
			foreach my $sample(@sample){
				$hash{$sample}{$key1} ||=0;
				print OUT "$hash{$sample}{$key1}\t"; 
			}
			print OUT "\n";
		}
		
	}
}
