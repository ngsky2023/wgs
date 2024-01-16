#!/usr/bin/env perl
use warnings;
use strict;
use FindBin '$Bin';
use Getopt::Long;
use File::Basename;

my $usage =
"
Usage:

  Options:
  -i|input      <FILE>   Input bam-file.
  -o|out        <FILE>   Output file.
  -s|samtools   <STR>    The all path of samtools,default is '/home/leiyoubing/software/samtools-0.1.19/samtools'.
  -r|read       <FILE>   R1 clean data for stat clean-reads, mut-lane is split ',',such as 'JMW-3_1_R1.fq.gz,JMW-3_2_R1.fq.gz'.
  -t|type       <STR>    'PE' or 'SE',default is 'PE'.
  -l|leng       <INT>    The length of reads,default is 100.
  -h|help                Help

For example:
        perl $0 -i sample.bam -o sample_QC.xls -s /home/leiyoubing/software/samtools-0.1.19/samtools
";

my ($input,$output,$samtools,$clean_data,$type,$leng,$help);
GetOptions(
  "i|input=s"=>\$input,
  "o|out=s"=>\$output,
  "s|samtools=s"=>\$samtools,
  "r|read=s"=>\$clean_data,
  "t|type=s"=>\$type,
  "l|leng=s"=>\$leng,
  "h|help"=>\$help
);
$samtools ||= "/share/public/software/samtools-1.4/samtools";
$type ||= 'PE';
$leng ||= 150;
die "$usage\n" if (!$input || !$output || !$clean_data || $help);

my $clean_reads=0;
my $clean_bases=0;
my $mapped_reads=0;
my $mapped_bases=0;
my $mapping_rate=0;
my $dup_reads=0;
my $dup_rate=0;
my $mismatch_bases=0;
my $mismatch_rate=0;

my @r1 = split /,/,$clean_data;
foreach my $r1(@r1){
	if($r1 =~ /.gz$/){
		my $tmp = `gzip -dc $r1|wc -l `;
		$clean_reads += $tmp;
		open IN , "gzip -dc $r1|";
	}else{
		my $tmp = `wc -l $r1`;
		$tmp = $1 if ($tmp =~ /^\s*(\d+)\s*/);
		$clean_reads += $tmp;
		open IN , "$r1";
	}
	<IN>;
	my $seq = <IN>;
	close IN;
	chomp $seq;
	$leng = length($seq);
}
if ($type =~ /PE/i){
	$clean_reads=$clean_reads/2;
}elsif($type =~ /SE/i){
	$clean_reads=$clean_reads/4;
}
$clean_bases = $clean_reads*$leng;

my ($insl , $insn) = (0 , 0);
open BAM,"$samtools view $input | " or die;
while(<BAM>)
{
	chomp;
	my $info=$_;
	my @info=split /\s+/,$info;
	unless($info[1] & 0x800){
		unless($info[1] & 0x100){
			unless($info[1] & 0x4){
				$mapped_reads++;
				if($info[1] & 0x400){$dup_reads++;}
				if($info=~/XC:i:(\d+)/){
					$mapped_bases+=$1;
				}else{
					$mapped_bases+=length($info[9]);
				}
				my $temp_m = 0;
				if($_=~/NM:i:(\d+)/){
					$temp_m = $1;
					$temp_m -= $1 if($info[5] =~ /(\d+)[ID]/);
				}
				$mismatch_bases+=$temp_m;
				if ($info[8] > 0 and $info[8] < 1000){
					$insl += $info[8];
					$insn++;
				}
			}
		}
	}
}
close BAM;

$mismatch_rate=$mismatch_bases/$mapped_bases;
$mapping_rate=$mapped_reads/$clean_reads;
$dup_rate=$dup_reads/$mapped_reads;
my $insert = int($insl/$insn);

my $name=basename($input);
my $sample=(split /\./,$name)[0];

open OUT,">$output" or die;
print OUT "Sample\t$sample\n";
print OUT "Clean reads\t$clean_reads\n";
print OUT "Clean bases(bp)\t$clean_bases\n";
print OUT "Mapped reads\t$mapped_reads\n";
print OUT "Mapped bases(bp)\t$mapped_bases\n";
print OUT "Mapping rate\t",sprintf("%.2f%%",100*$mapping_rate),"\n";
print OUT "Duplicate reads\t$dup_reads\n";
print OUT "Duplicate rate\t",sprintf("%.2f%%",100*$dup_rate),"\n";
print OUT "Mismatch bases(bp)\t$mismatch_bases\n";
print OUT "Mismatch rate\t",sprintf("%.2f%%",100*$mismatch_rate),"\n";
print OUT "Average insert size(bp)\t$insert\n";
close OUT;

