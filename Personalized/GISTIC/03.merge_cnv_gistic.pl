use strict;
use warnings;
use Data::Dumper;
#
my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
close IN;

my $cosmic="/share/work1/hanwj4457/softwore/GISTIC2/COSMIC_driver/COSMIC.driver.tier1";
open CO,$cosmic;
chomp(my @cosmic=<CO>);
close CO;

my $CGC="/share/work1/hanwj4457/database/CGC/gene_MoA.txt";
my %cgc;
open C,$CGC;
while(<C>){
	chomp;
	next if /^#/;
	my @line=split/\t/,$_;
	$cgc{$line[0]}=$line[1];

}
close C;
my $OG="/share/work1/hanwj4457/database/TUSON/OGs.txt";
my $TSG="/share/work1/hanwj4457/database/TUSON/TSGs.txt";
open O,$OG;
open T,$TSG;
chomp (my @og=<O>);
chomp (my @tsg=<T>);
close O;
close T;
my (%hgnc, %previous);
open IN, "/share/work1/hanwj4457/database/hg19/HGNC_gnames.txt" or die;
while (<IN>) {
        chomp;
        next if /^$/;
        my @f = split/\t/;
        $hgnc{$f[1]} = 1;
        my @old_name = split/, /, $f[4];
        if ($#old_name >= 0) {
                foreach (@old_name) {
                        $previous{$_} = $f[1];
                }
        }
}
close IN;



sub to_hgnc{
        my $gene = $_[0];
        my @gene = split/,;/, $gene;
        my $hgnc_name = '';
        foreach (@gene) {
                if (exists $hgnc{$_} and !exists $previous{$_}) {
                        $hgnc_name = $_;
                }elsif (exists $previous{$_} and !exists $hgnc{$_}) {
                        $hgnc_name = $previous{$_};
                }else{
                        $hgnc_name = '';
                }
        }
        $hgnc_name = $gene[0] if $hgnc_name eq '';
        return $hgnc_name;
}
###sample CNV
my %cnv;
my $path="/share/work1/hanwj4457/project/kongyan_jiayan_zhiduan+nianmo/analysis/mutation";
my $path1="/share/Onc_KYproject/hanwj4457/analysis/project/beizhong_mucosalMeanoma_20sample/20210414_17sample";
foreach my $sample(@sample){
	#my @file_list=glob("$path/*/cnv/$sample/purple*/$sample.purple.cnv.gene.tsv $path/*/*/cnv/$sample/*/$sample.purple.cnv.gene.tsv $path1/*/cnv/$sample/purple*/$sample.purple.cnv.gene.tsv");
        my @file_list=glob("/share/Onc_KYproject/hanwj4457/analysis/project/kongyan_jiayan_zhiduan+niamo20201209/acral/cnv_result/$sample.purple.cnv.gene.tsv");
	if(@file_list){
			open F,"@file_list" or die "@file_list is error\n";
			while(<F>){
				chomp;
				next if /^chromosome/;
				my @line=split/\t/;
				my $g1=$line[3];
				$g1=to_hgnc($g1);
				$line[4]=$line[4] <0 ? 0 : $line[4];
                                $line[5]=$line[5] <0 ? 0 : $line[5];
                                my $cnv=($line[4]+$line[5])/2;				
				$cnv{$sample}{$g1}=$cnv;
			}
		}
		close F;
}
####

open OUT,">result/my_interst_gene_list.xls";
my $amp="result/amp_genes.conf_95.t.txt";
my $head2;
open N,$amp;
while(<N>){
	chomp;
	my @line=split/\t/;
	if($.==1){
		$head2=join "\t",@line[0..3];
		my $head3=join "\t",@sample;
		print OUT "Type\tGENE\tTUSON\tCGC\tCOSMIC\t$head2\t$head3\n";
	}else{
		# my $gene=join "\t",@line[4..-1];
		foreach my $i (4..$#line){
		# foreach my $gene(@line[4..-1]){
			my $tuson=0;
			my $gene_MoA=0;
			my $cosmic_gene=0;
			my $gene=$line[$i];
			my $g1=$gene;
			$g1=to_hgnc($g1);
			if(grep /^$g1$/,@og){
				$tuson="OG";
			}elsif(grep /^$g1$/,@tsg){
				$tuson="TSG";
			}
			if(defined $cgc{$g1}){
				$gene_MoA=$cgc{$g1};
			}
			if(grep /^$g1$/,@cosmic){
				$cosmic_gene="cosmic";
			}
			# next unless grep{$_ eq $1} @gene;
			my $out1=join "\t",@line[0..3];
			print OUT "Amp\t$g1\t$tuson\t$gene_MoA\t$cosmic_gene\t$out1\t";
			foreach my $sample(@sample){
				$cnv{$sample}{$g1}||=0;
				print OUT "$cnv{$sample}{$g1}\t";
			}
			print OUT "\n";
		}
	}
	
	
}
close N;
my $del="result/del_genes.conf_95.t.txt";
open N,$del;
while(<N>){
	chomp;
	my @line=split/\t/;
	if($.==1){
		next;
	}else{
		# my $gene=join "\t",@line[4..-1];
		foreach my $i (4..$#line){
		# foreach my $gene(@line[4..-1]){
			my $tuson=0;
			my $gene_MoA=0;
			my $gene=$line[$i];
			my $cosmic_gene=0;
			my $g1=$gene;
			$g1=to_hgnc($g1);
			if(grep /^$g1$/,@og){
				$tuson="OG";
			}elsif(grep /^$g1$/,@tsg){
				$tuson="TSG";
			}
			if(defined $cgc{$g1}){
				$gene_MoA=$cgc{$g1};
			}
			if(grep /^$g1$/,@cosmic){
				$cosmic_gene="cosmic";
			}
			# next unless grep{$_ eq $gene} @gene;
			my $out1=join "\t",@line[0..3];
			print OUT "Del\t$g1\t$tuson\t$gene_MoA\t$cosmic_gene\t$out1\t";
			foreach my $sample(@sample){
				$cnv{$sample}{$g1}||=0;
				print OUT "$cnv{$sample}{$g1}\t";
			}
			print OUT "\n";
			
		}
	}
	
	
}





#
=cut
my %cnv;
my $path="/share/work1/hanwj4457/project/kongyan_jiayan_zhiduan+nianmo/analysis/mutation";
my $path1="/share/Onc_KYproject/hanwj4457/analysis/project/beizhong_mucosalMeanoma_20sample/20210414_17sample";
foreach my $sample(@sample){
	#my @file_list=glob("$path/*/cnv/$sample/purple*/$sample.purple.cnv.gene.tsv $path/*/*/cnv/$sample/*/$sample.purple.cnv.gene.tsv $path1/*/cnv/$sample/purple*/$sample.purple.cnv.gene.tsv");
        my @file_list=glob("/share/Onc_KYproject/hanwj4457/analysis/project/kongyan_jiayan_zhiduan+niamo20201209/acral/cnv_result/$sample.purple.cnv.gene.tsv");
	if(@file_list){
			open F,"@file_list";
			while(<F>){
				chomp;
				my @line=split/\t/;
				if(grep { $_ eq $line[3] } @gene){
					my $cn=($line[4]+$line[5])/2;
					$cnv{$sample}{$line[3]}=$cn;
				}
			}
		}
		close F;
}

my $gistic="./result/all_lesions.conf_95.txt";
open OUT,">combin_cnv_gistic.xls";
open F,$gistic;
my @sample2;
while(<F>){
	chomp;
	my @line=split/\t/;
	
	if($.==1){
		my $head2=join "\t",@line[0..8];
		@sample2=@line[9..-1];
		print OUT "$head1\t$_\n";
		next;
	}else{
		next unless $_=~/CN values/;
		$region{$line[1]}||="NA";
		$hash{$line[1]}||="NA";
		# print "$region{$line[1]}\n";
		# print OUT "$hash{$line[1]}\t@line[0..8]\t";
		# my $relation=$region{$line[1]};
		# foreach my $sample2(@sample2){
			# print OUT "$cnv{$sample2}{$relation}\t";
		# }
		# print OUT "\n";
	}
	
	
}
=cut
