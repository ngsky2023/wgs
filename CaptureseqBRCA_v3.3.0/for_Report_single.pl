use warnings;
use strict;

die "perl $0 <in.level> <out.core> <out.candidate>" unless @ARGV==3;

open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
open OUT2,">$ARGV[2]";
print OUT "SampleName\tGene\tType\tPosition\tTranscript\tExon\tHGVS(.c)\tHGVS(.p)\tGenotype\tAlleleFrequency\tDisease\tClinicalLevel\tEvidence\tAAchange\n";
print OUT2 "SampleName\tGene\tType\tPosition\tTranscript\tExon\tHGVS(.c)\tHGVS(.p)\tGenotype\tAlleleFrequency\tDisease\tClinicalLevel\tEvidence\tAAchange\n";
while(<IN>){
	chomp;
	my @atm=split /\t/;
	my ($gene,$chr,$start,$end,$ExonicFunc,$AAchange,$genotype,$frequence,$clinvar,$level,$evidence)=($atm[14],$atm[3],$atm[4],$atm[5],$atm[16],$atm[17],$atm[8],$atm[28],$atm[26],$atm[0],$atm[1]);
	my ($sample,$transcript,$funcref,$exon,$hc,$hg)=($atm[2],$atm[130],$atm[13],$atm[131],$atm[132],$atm[133]);
	if($level =~ /pathogenic/i){
		print OUT "$sample\t$gene\t$ExonicFunc\t$chr:$start-$end\t$transcript\t$funcref\t$exon\t$hc\t$hg\t$genotype\t$frequence\t$clinvar\t$level\t$evidence\t$AAchange\n";
	}elsif($level =~ /Undefined/ || $level =~ /Uncertain/){
		#print OUT2 "$_\n";
		print OUT2 "$sample\t$gene\t$ExonicFunc\t$chr:$start-$end\t$transcript\t$funcref\t$exon\t$hc\t$hg\t$genotype\t$frequence\t$clinvar\t$level\t$evidence\t$AAchange\n";
	}
}
close IN;
close OUT;
close OUT2;
