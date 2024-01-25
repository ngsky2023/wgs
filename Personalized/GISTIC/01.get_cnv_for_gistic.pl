use strict;
use warnings;
use Cwd;
use Data::Dumper;
use Math::Complex;
my %hash;
my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
#chdir "/share/work1/hanwj4457/project/kongyan_jiayan_zhiduan+nianmo/mucosal/kataegis/result";
open OUT,">01.cn.txt";
print OUT "SampleName\tChr\tStart\tEnd\tmak.num\tmedian(logCN-1)\n";
#print OUT "SampleName\tChr\tStart\tEnd\tmedian(logCN-1)\n";
foreach my $sample(@sample){
        #`ln -s $snv_path/$sample.final.xls $sample.final.xls`;
	open S,"../../acral/cnv_result/$sample.purple.cnv.somatic.tsv";
	while(<S>){
		chomp;
		next if /^chromosome/;
		my @line=split/\t/,$_;
		$line[0]=~s/chr//g;
		$line[3]=$line[3] < 0 ? 0.01 :$line[3];
		my $lrr=logn($line[3],2)-1;
		#my $out=join "\t",@line[0,1,2,10];
		my $out=join "\t",@line[0,1,2,4];
		print OUT "$sample\t$out\t$lrr\n";
		#print $line[3]."\t$lrr\n";
	}
}













