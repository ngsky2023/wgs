use strict;
use warnings;
use Cwd;
use Data::Dumper;
my %hash;
my $snv_path="/path/snv-indel/result";
my $samplelist="../clinical/sample.list";
open IN,$samplelist;
chomp(my @sample=<IN>);
my $result= "./result";
`mkdir $result` unless -d $result;
chdir $result;
foreach my $sample(@sample){
	next if (-s "${sample}_Kataegis.tsv");
	open S,">$sample.kata.sh";
	open OUT,">IN/$sample.maf";
	print OUT "Hugo_Symbol\tChromosome\tStart_position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n";
        my (@title,@snv);
	my $num=0;
        open F,"$snv_path/$sample/$sample.final.xls" or die "$sample file is error \n";
        while(<F>){
                chomp;
                next if /^\s+$/;
		next if /^SampleName/;
		my @line=split/\t/,$_;
		$line[14]=~s/.*synonymous SNV/Missense_Mutation/g;
		$line[14]=~s/nonframeshift deletion/In_Frame_Del/g;
		$line[14]=~s/frameshift insertion/Frame_Shift_Ins/g;
		$line[14]=~s/frameshift deletion/In_Frame_Del/g;
		$line[14]=~s/exonic;splicing/Splice_Site/g;
		$line[14]=~s/intergenic/IGR/g;
		$line[14]=$line[11] if $line[14] eq '.';
		my $type;
		if($line[14] =~/Ins/){
			$type='INS';
		}elsif($line[14]=~/Del/){
			$type='DEL';
		}else{
			$type='SNP';
		}
		# next unless (length($line[4]) ==1 && length($line[5])==1);
		# next unless $line[4]=~/[ATCG]/;
		# next unless $line[5]=~/[ATCG]/;
		print OUT "$line[12]\t$line[1]\t$line[2]\t$line[3]\t\+\t$line[14]\t$type\t$line[4]\t$line[4]\t$line[5]\t$line[0]\n";
		
        }
	print S "./kataegis.R ./IN/$sample.maf $sample\n";
	# print "$sample\n";
	print my $jobs=`qsub -cwd -l vf=2G $sample.kata.sh`;
}













