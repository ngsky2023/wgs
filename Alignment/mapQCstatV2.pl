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
my %soft;
open IN , "$Bin/soft.list";
while (<IN>){
	chomp;
	my ($k , $d) = split /\t/ , $_;
	$soft{$k} = $d;
}
close IN;

$samtools ||= $soft{samtools};
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

open IN , "$clean_data";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	if ($F[0] eq 'total_bases'){
		$clean_bases = $F[1];
	}elsif ($F[0] eq 'total_reads'){
		$clean_reads = $F[1];
	}
}
close IN;


my %ins;
my ($insl , $insn) = (0 , 0);
my %insize;
srand();
my @col = qw/mapped_reads dup_reads mapped_bases mismatch_bases/;
my %stat;
open BAM,"$samtools view $input | " or die;
while(<BAM>)
{
	chomp;
	my $info=$_;
	my @info=split /\s+/,$info;
	unless($info[1] & 0x800){
		unless($info[1] & 0x100){
			unless($info[1] & 0x4){
				my $group = 'nogroup';
				if (/RG:Z:(\S+)/){
					$group = $1;
				}
				$mapped_reads++;
				$stat{$group}{mapped_reads}++;
				if($info[1] & 0x400){
					$dup_reads++;
					$stat{$group}{dup_reads}++;
				}elsif ($info[8] > 30 and $info[8] < 1000){
					$insl += $info[8];
					$ins{$info[8]}++;
					$insn++;
					if ($info[2] =~ /chr[\dXY]+/){
						my @stat = stat($input);
						my $size = $stat[7];
						my $randb = int(($size/1000000000) * (1000/100));
						$randb = 1 if $randb < 1;
						my $randnum = int(rand($randb));
						if ($randnum == 0){
							push @{$insize{Nuclear}} , $info[8];
						}
					}elsif ($info[2] =~ /chrM/){
						$info[2] = "MT";
						push @{$insize{$info[2]}} , $info[8];
					}
				}
				my $temp_m = 0;
				if($_=~/NM:i:(\d+)/){
					$temp_m = $1;
				}

				while ($info[5] =~ /(\d+)([MID])/g){
					my ($mn , $nty) = ($1 , $2);
					if ($nty eq 'M'){
						$mapped_bases += $mn;
						$stat{$group}{mapped_bases} += $mn;
					}elsif ($nty eq 'I'){
						$mapped_bases += $mn;
						$stat{$group}{mapped_bases} += $mn;
						$temp_m -= $mn;
					}elsif ($nty eq 'D'){
						$temp_m -= $mn;
					}
				}
				$mismatch_bases+=$temp_m if $temp_m > 0;
				$stat{$group}{mismatch_bases} +=$temp_m if $temp_m > 0;
			}
		}
	}
}
close BAM;

open INSTXT , ">$output.insert.txt";
print INSTXT "Class\tInsert\n";
for my $chr (sort keys %insize){
	for my $size (@{$insize{$chr}}){
		print INSTXT "$chr\t$size\n";
	}
}
close INSTXT;

$mismatch_rate=$mismatch_bases/$mapped_bases;
$mapping_rate=$mapped_reads/$clean_reads;
$dup_rate=$dup_reads/$mapped_reads;
my $insert = int($insl/$insn);
my $peakIns = (sort {$ins{$b}<=>$ins{$a}} keys %ins)[0];

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
print OUT "Peak of insert size(bp)\t$peakIns\n";
close OUT;

open OUT , ">$output.table";
print OUT "group\t" , join("\t" , @col) , "\n";
for my $group (sort keys %stat){
	print OUT $group;
	for my $tt (@col){
		if (exists $stat{$group}->{$tt}){
			print OUT "\t" , $stat{$group}->{$tt};
		}else{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close OUT;

