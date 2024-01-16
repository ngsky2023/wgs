#!/usr/bin/perl -w
use strict;

my ($input,$gff,$output) = @ARGV;
die "ERRO!\nUsage:\n\tperl $0 Sample.predSV.txt genes.gff output.xls\n\n" if (@ARGV<3);
my (%data,%genes) = ();
my ($gene,$gene_name) = ("","");

open GFF,"$gff" || die "$!";
while(my $line = <GFF>){
	chomp($line);
	next if ($line =~ /^\s*$/);
	next if ($line =~ /^#/);
	my @arr = split /\t/,$line;
	if ($arr[2] =~ /gene/i){
		$gene = "$arr[3]-$arr[4]";
		if ($arr[8] =~ /ID=(.+?);/){
			$gene_name = $1;
		}else{
			$gene_name = "-";
		}
		$genes{$arr[0]}{$arr[3]}{$arr[4]} = $gene_name;
		next;
	}else{
		my $fuc = "$arr[2]-$arr[3]-$arr[4]";
		$data{$arr[0]}{$gene}{$fuc} = $gene_name;
	}
}
close GFF;

open IN,"$input" || die "$!";
open OUT,">$output" || die "$!";
#print OUT "SV_type\tleft_chr\tleft_pos\tleft_strand\tleft_soft-clipped_reads\tright_chr\tright_pos\tright_strand\tright_soft-clipped_reads\tregion\tgene\n";
while(my $line=<IN>){
	chomp($line);
	my @arr = split /\t/,$line;
	###filter
	if ($arr[0] eq $arr[4]){
		next if ($arr[3]<4 || $arr[7]<4);
		my $length = abs($arr[5] - $arr[1]);
		next if ($arr[8] ne 'INS' and ($length < 50 || $length > 1000000));
		next if ($arr[8] eq 'INS' and $length > 1000);
	}
	###annotate
	if ($arr[0] ne $arr[4]){
		print OUT "$arr[8]\t$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[5]\t$arr[6]\t$arr[7]\t-\t-\n";
		next;
	}else{
		my $flag = 0;
		my ($tmp_region,$tmp_gene) = ('','');
		foreach my $key1(keys %{$data{$arr[0]}}){
			my ($start,$end) = split /-/,$key1;
			if (($arr[1] >= $start && $arr[1] < $end) || ($arr[5] > $start && $arr[5] < $end) || ($arr[8] eq 'DEL' and $arr[1] <= $start and $end <= $arr[5])){
				$flag = 1;
				$tmp_gene .= "$genes{$arr[0]}{$start}{$end};";
				my $region = "";
				foreach my $key2(keys %{$data{$arr[0]}{$key1}}){
					my ($re,$s,$e) = split /-/,$key2;
					if (($arr[1] >= $s && $arr[1] < $e) || ($arr[5] > $s && $arr[5] < $e) || ($arr[8] eq 'DEL' and $arr[1] <= $s and $e <= $arr[5])){
						$region .= "$re,";
					}
				}
				if ($region eq ""){
					$tmp_region .= "Intron;";
				}else{
					if ($region =~ /CDS/i){
						$tmp_region .= "CDS;";
					}elsif($region =~ /UTR/i){
						$tmp_region .= "UTR;";
					}elsif($region =~ /exon/i){
						$tmp_region .= "exon;";
					}elsif($region =~ /mRNA/i){
						$tmp_region .= "mRNA;";
					}else{
						$tmp_region .= "$region;";
					}
				}
			}
		}
		if ($flag == 0){
			print OUT "$arr[8]\t$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[5]\t$arr[6]\t$arr[7]\tIntergenic\t-\n";
		}else{
			$tmp_region =~ s/[;,]$//;
			$tmp_gene =~ s/[;,]$//;
			print OUT "$arr[8]\t$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[5]\t$arr[6]\t$arr[7]\t$tmp_region\t$tmp_gene\n";
		}
	}
}
close IN;
close OUT;

