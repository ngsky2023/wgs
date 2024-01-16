#!/usr/bin/perl -w
use strict;

if (@ARGV < 2){
        die "perl $0 input_variation.vcf prefix_for_out\n";
}
my ($input,$prefix) = @ARGV;

open IN,"$input" || die "$!";
open SNP,">$prefix.snp.vcf";
open INDEL,">$prefix.indel.vcf";
while(my $line = <IN>){
        chomp($line);
        if ($line =~ /^#/){
                print SNP "$line\n";
                print INDEL "$line\n";
                next;
        }
        my @arr = split /\t/,$line;
        my $len1 = length $arr[3];
        my $len2 = length $arr[4];
        if ($len1 != $len2 ){
                print INDEL "$line\n";
        }else{
                print SNP "$line\n";
        }
}
close IN;
close SNP;
close INDEL;
