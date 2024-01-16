#! /usr/bin/perl -w
#perl annovar.3mode.pl  $sample.annovar.3mode.tmp.xls $sample.annovar.3mode.xls\n";
die "perl $0 <in> <out>" unless @ARGV==2;

my $in=shift;
my $out=shift;
open IN, "<$in" or die $!;
open OUT,">$out" or die $!;
while (<IN>){
    chomp;
    if ($_ =~ /SampleName/){
        print OUT "$_\n";
        next;
    }
    my @dat =split /\t/,$_;
    if($dat[15] eq ",,,,,,"){
         $dat[15] = ".";
         print  OUT join ("\t",@dat),"\n";
    }elsif ($dat[15] =~/(.*)\,$/){
        #print $1;
        $dat[15] = $1;
        print  OUT join ("\t",@dat),"\n";
    }else{
        print OUT "$_\n";
    }
}
close IN;
close OUT;