use warnings;
use strict;

die "perl $0 <in.level> <out.level1> <out.level2><out.level3>" unless @ARGV==4;
`sort -k6 -rn -o $ARGV[0] $ARGV[0]`;
open IN,"<$ARGV[0]";
open OUT,">$ARGV[1]";
open OUT2,">$ARGV[2]";
open OUT3,">$ARGV[3]";
print OUT "SampleName\tGene\tTranscript\tExon\tHGVS(.c)\tHGVS(.p)\tGenotype\tlevel\tAlleleFrequency\tEvidence\ttype\tPos\tAAchange\tClinicalLevel\n";
print OUT2 "SampleName\tGene\tTranscript\tExon\tHGVS(.c)\tHGVS(.p)\tGenotype\tlevel\tAlleleFrequency\tEvidence\ttype\tPos\tAAchange\tClinicalLevel\n";
print OUT3 "SampleName\tGene\tTranscript\tExon\tHGVS(.c)\tHGVS(.p)\tGenotype\tlevel\tAlleleFrequency\tEvidence\ttype\tPos\tAAchange\tClinicalLevel\n";

while(<IN>){
	chomp;
	my @atm=split /\t/;
	my ($function,$gene,$chr,$start,$end,$genedetail,$ExonicFunc,$AAchange,$genotype,$frequence,$clinvar,$level,$evidence)=($atm[13],$atm[14],$atm[3],$atm[4],$atm[5],$atm[15],$atm[16],$atm[17],$atm[8],$atm[28],$atm[26],$atm[0],$atm[1]);
	my ($sample,$transcript,$exon,$hc,$hp)=($atm[2],$atm[96],$atm[97],$atm[98],$atm[99]);
        if($gene=~/BRCA1/ && $function=~/exonic|splicing/)
        {
	if($level =~ /pathogenic/i){
           print OUT "$sample\t$gene\t$transcript\t$exon\t$hc\t$hp\t$genotype\t$level\t$frequence\t$evidence\t$ExonicFunc\t$chr:$start-$end\t$AAchange\t$clinvar\n";
	}
        elsif($level =~ /Undefined/ || $level =~ /Uncertain/){
	   print OUT2 "$sample\t$gene\t$transcript\t$exon\t$hc\t$hp\t$genotype\t$level\t$frequence\t$evidence\t$ExonicFunc\t$chr:$start-$end\t$AAchange\t$clinvar\n";
	}
        elsif($level =~ /benign/i){
           print OUT3 "$sample\t$gene\t$transcript\t$exon\t$hc\t$hp\t$genotype\t$level\t$frequence\t$evidence\t$ExonicFunc\t$chr:$start-$end\t$AAchange\t$clinvar\n";
       }
       }
}
close IN;
`sort -k6 -n -o $ARGV[0] $ARGV[0]`;
open IN,"<$ARGV[0]";
while(<IN>){
        chomp;
        my @atm=split /\t/;
        my ($function,$gene,$chr,$start,$end,$genedetail,$ExonicFunc,$AAchange,$genotype,$frequence,$clinvar,$level,$evidence)=($atm[13],$atm[14],$atm[3],$atm[4],$atm[5],$atm[15],$atm[16],$atm[17],$atm[8],$atm[28],$atm[26],$atm[0],$atm[1]);
        my ($sample,$transcript,$exon,$hc,$hp)=($atm[2],$atm[96],$atm[97],$atm[98],$atm[99]);
        if($gene=~/BRCA2/ && $function=~/exonic|splicing/)
        {
        if($level =~ /pathogenic/i){
           print OUT "$sample\t$gene\t$transcript\t$exon\t$hc\t$hp\t$genotype\t$level\t$frequence\t$evidence\t$ExonicFunc\t$chr:$start-$end\t$AAchange\t$clinvar\n";
        }
        elsif($level =~ /Undefined/ || $level =~ /Uncertain/){
           print OUT2 "$sample\t$gene\t$transcript\t$exon\t$hc\t$hp\t$genotype\t$level\t$frequence\t$evidence\t$ExonicFunc\t$chr:$start-$end\t$AAchange\t$clinvar\n";
        }
        elsif($level =~ /benign/i){
           print OUT3 "$sample\t$gene\t$transcript\t$exon\t$hc\t$hp\t$genotype\t$level\t$frequence\t$evidence\t$ExonicFunc\t$chr:$start-$end\t$AAchange\t$clinvar\n";
       }
       }
}
close IN;

close OUT;
close OUT2;
close OUT3;
