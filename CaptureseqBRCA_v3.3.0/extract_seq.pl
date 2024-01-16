#! /usr/bin/perl -w
# to extract wrong mapped sequence while most of the hg19 is masked by "N",except of genes in the panel
system("samtools view $ARGV[0] > ./sample1.sam");
system("samtools view $ARGV[1] > ./sample2.sam");
open IN,"./sample1.sam"||die; ### before 
while(<IN>){
        chomp;
        @tt=split/\s+/;
        $key=$tt[0]."\t".$tt[2]."\t".$tt[3];
        $tit{$key}++
}
open IN1,"./sample2.sam"||die; ## after
while(<IN1>){
        chomp;
        @tt=split/\s+/;
        $key=$tt[0]."\t".$tt[2]."\t".$tt[3];
        if( ! defined $tit{$key}){
                $extra{$tt[0]}++;
        }
}
open OUT1,">$ARGV[4]"||die;
open OUT2,">$ARGV[5]"||die;
open IN2,"gzip -dc $ARGV[2] |"||die; ##fq file
while(<IN2>){
        chomp;
        @tt=split/\s+/;
	$tt[0]=~s/^\@//;
        $seq=<IN2>;
        <IN2>;
        $qual=<IN2>;
        if( defined $extra{$tt[0]}){
                if(length $seq >100){
			$seq=substr($seq,0,100);
			$qual=substr($qual,0,100);
			print OUT1 "\@$tt[0] $tt[1]\n$seq\n\+\n$qual\n";
		}else{
			print OUT1 "\@$tt[0] $tt[1]\n$seq\+\n$qual";
        	}
	}
}
open IN3,"gzip -dc $ARGV[3] |"||die; ##fq file
while(<IN3>){
        chomp;
        @tt=split/\s+/;
        $tt[0]=~s/^\@//;
        $seq=<IN3>;
        <IN3>;
        $qual=<IN3>;
        if( defined $extra{$tt[0]}){
                if((length $seq) >100){
                        $seq=substr($seq,0,100);
                        $qual=substr($qual,0,100);
			print OUT2 "\@$tt[0] $tt[1]\n$seq\n\+\n$qual\n";
		}else{
                	print OUT2 "\@$tt[0] $tt[1]\n$seq\+\n$qual";
        	}
	}
}

