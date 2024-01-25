#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
my ($input , $out, $help);

GetOptions(
	"i|input=s"     =>      \$input,
	"o|out=s" => \$out,
	"help"  =>      \$help,
);

if ($help or ! $input or ! $out){
	&help;
	exit;
}


open FILE,$input || die "Can't open the file: $!\n";
open OUT,">$out" || die "Can't write: $!\n";

print OUT "gene\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tdbSNP_Val_Status\tpatient\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\tValidation_Status\tMutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\n";
my %hash;
my $title1 = <FILE>;
chomp($title1);
my @title1 = split /\t/, $title1;
while(<FILE>){
	next if (/^\s$/);
	chomp;
	my @line = split /\t/,$_;
	my $a = (split /_/, $line[0])[-1];
	$line[0] = $a;
	$hash{$title1[$_]} = $line[$_] foreach (0..$#line);
	my $gene_symbol = $hash{"Gene.refGene"};
	$gene_symbol = (split ';', $hash{"Gene.refGene"})[0] if($hash{"Gene.refGene"} =~ /;/);
	$gene_symbol = (split ',', $hash{"Gene.refGene"})[0] if($hash{"Gene.refGene"} =~ /,/);
	my $chr = $hash{"Chr"};
	my $start = $hash{"Start"};
	my $end = $hash{"End"};
	my $variant_type;
	if($hash{"Ref"}=~/[ATGC]/ && $hash{"Alt"}=~/[ATGC]/){
		$variant_type = "SNP"
	}elsif($hash{"Ref"}=~/[ATGC]+/ && $hash{"Alt"}=~/-/){
		$variant_type = "DEL";
	}elsif($hash{"Ref"}=~/-/ && $hash{"Alt"}=~/[ATGC]+/){
		$variant_type = "INS";
	}
	my $variant_classification;
	if($hash{"Func.refGene"}=~/splicing/){
		$variant_classification="Splice_Site";
	}elsif($hash{"ExonicFunc.refGene"}eq"frameshift deletion"){
		$variant_classification="Frame_Shift_Del";
	}elsif($hash{"ExonicFunc.refGene"}eq"frameshift insertion"){
		$variant_classification="Frame_Shift_Ins";
	}elsif($hash{"ExonicFunc.refGene"}eq"nonframeshift deletion"){
		$variant_classification="In_Frame_Del";
	}elsif($hash{"ExonicFunc.refGene"}eq"nonframeshift insertion"){
		$variant_classification="In_Frame_Ins";
	}elsif($hash{"ExonicFunc.refGene"}eq"nonsynonymous SNV"){
		$variant_classification="Missense_Mutation";
	}elsif($hash{"ExonicFunc.refGene"}eq"stopgain" or $hash{"ExonicFunc.refGene"}eq"stoploss"){
		$variant_classification="Nonsense_Mutation";
	}elsif($hash{"ExonicFunc.refGene"}eq"synonymous SNV"){
		$variant_classification="Silent";
	}elsif($hash{"ExonicFunc.refGene"} eq "unknown"){
		next;
	}elsif($hash{"ExonicFunc.refGene"}eq"."){
		if($hash{"Func.refGene"}eq"intronic"){
			$variant_classification="Intron";
		}elsif($hash{"Func.refGene"}=~/ncRNA/){
			$variant_classification="RNA";
		}elsif($hash{"Func.refGene"}=~/UTR5/){
			$variant_classification="5'UTR";
		}elsif($hash{"Func.refGene"}eq"UTR3"){
			$variant_classification="3'UTR";
		}elsif($hash{"Func.refGene"}eq"upstream"){
			$variant_classification="upstream";
		}elsif($hash{"Func.refGene"}eq"downstream"){
			$variant_classification="downstream";
		}elsif($hash{"Func.refGene"}eq"upstream;downstream"){
			$variant_classification="upstream;downstream";
		}elsif($hash{"Func.refGene"}eq"intergenic"){
			$variant_classification="IGR";
		}
	}
	my $ref = $hash{"Ref"};
	my $alt = $hash{"Alt"};
	my $sample = $hash{"SampleName"};
	$sample =~ s/-/_/g;
	my $blank ="-\t-";
	my $blank2 ="\t-"x 18;
	print OUT "$gene_symbol\t0\t$blank\t$chr\t$start\t$end\t+\t$variant_classification\t$variant_type\t$ref\t$ref\t$alt\t$blank\t$sample$blank2\n";
}
close OUT;
close FILE;
