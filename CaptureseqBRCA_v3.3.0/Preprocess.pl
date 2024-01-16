#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use Cwd;

my $pwd = getcwd;
my $usage =
"
Usage:

  Options:
  -i|input      <FILE>   Input file record sample name and flowcell name.
  -o|outdir     <DIR>    Output file direcory, if this is not exists created, default(./).
  -f|flexbar	<software> direction of software flexbar
  -r|rd	        <RD dir> direction of RD
  -a|adapter	<STR>    adapter sequence
  -d|random     <STR>    example:1000000,2000000,5000000
  -h|help                Help

For example:
        perl $0 -i flowcell.info -o ./analysis -f flexbar -r rd dir -a adapter -d random
";

my ($input,$outdir,$config,$step,$cut,$flexbar,$rd,$adapter,$random,$help);
GetOptions(
  "i|input=s"=>\$input,
  "o|outdir=s"=>\$outdir,
  "f|flexbar=s"=>\$flexbar,
  "r|rd=s"=>\$rd,
  "d|random=s"=>\$random,
  "a|adapter=s"=>\$adapter,
  "help"=>\$help
);
$outdir ||= $pwd;
$outdir = abs_path($outdir);
$flexbar ||="/share/software/software/flexbar_v2.4_linux64/flexbar";
$flexbar ||="/share/work1/staff/xuxiong/software/flexbar_v2.5_linux64/flexbar";
$rd ||="/share/seq_dir1/Item1/RD";
$adapter ||="/home/sunfl/workdir/adaptor/nextera";
$random ||="NO";
if (!$input  || $help){
        die "$usage\n";
}
system("mkdir -m 755 $outdir/Input_fq") unless (-d "$outdir/Input_fq");
system("mkdir -m 755 -p $outdir/Input_fq/1.Trim") unless (-e "$outdir/Input_fq/1.Trim");
system("mkdir -m 755 -p $outdir/Input_fq/2.Merge") unless (-e "$outdir/Input_fq/2.Merge");
system("mkdir -m 755 -p $outdir/Input_fq/3.Random") unless (-e "$outdir/Input_fq/3.Random");

my $name=$input;
$name=~s/.*\///g;
open INFO,">$outdir/Input_fq/sample.$name.txt";
if($random eq "NO"){
	my @samp=();
	open FL,"$input" || die "$!";
	while(my $line = <FL>){
		chomp($line);
		my @items = split /\s+/,$line;
		my $i=0;
		while($i<=$#items-1){
			push @samp,[$items[$i],$items[$i+1]];
			$i+=2;
		}
		if($#items>1){
			&Preprocess($rd,$outdir,$flexbar,$adapter,$random,@samp);
			print INFO "$items[1]\tL1\t$outdir/Input_fq/2.Merge/$items[1].L1_1.fastq.gz\t$outdir/Input_fq/2.Merge/$items[1].L1_2.fastq.gz\n";
		}else{
			&Preprocess($rd,$outdir,$flexbar,$adapter,$random,@samp);
			print INFO "$items[1]\tL1\t$outdir/Input_fq/1.Trim/$items[1].L1_1.fastq.gz\t$outdir/Input_fq/1.Trim/$items[1].L1_2.fastq.gz\n";
		}@samp=();
	}
	close FL;
}else{
	my @ran=split(/\,/,$random);
	foreach my $j(0..$#ran){
		my @samp=();
		open FL,"$input" || die "$!";
		while(my $line = <FL>){
			chomp($line);
			my @items = split /\s+/,$line;
			my $i=0;
			while($i<=$#items-1){
				push @samp,[$items[$i],$items[$i+1]];
				$i+=2;
			}
			if($#items>1){
				&Preprocess($rd,$outdir,$flexbar,$adapter,$ran[$j],@samp);
				print INFO "$items[1].$ran[$j]\tL1\t$outdir/Input_fq/3.Random/$items[1].L1_1.fastq.$ran[$j].gz\t$outdir/Input_fq/3.Random/$items[1].L1_2.fastq.$ran[$j].gz\n";
			}else{
				&Preprocess($rd,$outdir,$flexbar,$adapter,$ran[$j],@samp);
				print INFO "$items[1].$ran[$j]\tL1\t$outdir/Input_fq/3.Random/$items[1].L1_1.fastq.$ran[$j].gz\t$outdir/Input_fq/3.Random/$items[1].L1_2.fastq.$ran[$j].gz\n";
			}@samp=();
		}
		close FL;
	}
}
close INFO;



sub Preprocess{
	my ($rd,$outdir,$flexbar,$adapter,$random,@samp)=@_;
	my ($index,$fqlist1,$fqlist2)=("","","");
	if($random ne "NO"){
		$index=$samp[0][1]."\.$random";
	}else{
		$index =$samp[0][1];
	}
	open TA,">$outdir/Input_fq/1.Trim/Preprocess_$index.sh";
        print TA "#!/bin/bash\n";
        print TA "echo ===start at : `date` ===\n";
	foreach my $i (0..$#samp){
		my $fq1="$rd/*$samp[$i][0]/$samp[$i][1]"."_*.R1.clean.fastq.gz";
       		my $fq2="$rd/*$samp[$i][0]/$samp[$i][1]"."_*.R2.clean.fastq.gz";
		if($random ne "NO" && $#samp==0){
			print TA "$flexbar -at 2 -u 1 -n 5 -m 36 -ao 5 -f i1.8 -r $fq1  -p $fq2   $adapter  -t $outdir/Input_fq/1.Trim/$samp[$i][1].L1 -z GZ\n";
			print TA "perl $Bin/random_read.pl $outdir/Input_fq/1.Trim/$samp[$i][1].L1_1.fastq.gz   $outdir/Input_fq/1.Trim/$samp[$i][1].L1_2.fastq.gz $random\n";
			print TA "mv $outdir/Input_fq/1.Trim/$samp[$i][1].L1_1.fastq.$random.gz $outdir/Input_fq/3.Random/\n";
			print TA "mv $outdir/Input_fq/1.Trim/$samp[$i][1].L1_2.fastq.$random.gz $outdir/Input_fq/3.Random/\n";
			#print TA "gzip $outdir/Input_fq/3.Random/$samp[$i][1].L1_1.fastq.$random\n";
			#print TA "gzip $outdir/Input_fq/3.Random/$samp[$i][1].L1_2.fastq.$random\n";
			#print TA "gzip $outdir/Input_fq/1.Trim/$samp[$i][1].L1_1.fastq.$random\n";
			#print TA "gzip $outdir/Input_fq/1.Trim/$samp[$i][1].L1_2.fastq.$random\n";
		}elsif($random eq "NO" && $#samp==0){
			print TA "$flexbar -at 2 -u 1 -n 5 -m 36 -ao 5 -f i1.8 -r $fq1  -p $fq2   $adapter  -t $outdir/Input_fq/1.Trim/$samp[$i][1].L1 -z GZ\n";
		}elsif($#samp >0){
			print TA "$flexbar -at 2 -u 1 -n 5 -m 36 -ao 5 -f i1.8 -r $fq1  -p $fq2   $adapter  -t $outdir/Input_fq/1.Trim/$samp[$i][1].L1 \n";
			#print TA "perl $Bin/cut_seq.pl $fq1  $flexbar $outdir/Input_fq/1.Trim/$samp[$i][1].L1_1.fastq 1\n";
			#print TA "perl $Bin/cut_seq.pl $fq2  $flexbar $outdir/Input_fq/1.Trim/$samp[$i][1].L1_2.fastq 1\n";
			$fqlist1.="$outdir/Input_fq/1.Trim/$samp[$i][1].L1_1.fastq ";
			$fqlist2.="$outdir/Input_fq/1.Trim/$samp[$i][1].L1_2.fastq ";
		}
	}
	if($#samp >0){
		print TA "cat $fqlist1|gzip - >$outdir/Input_fq/2.Merge/$samp[0][1].L1_1.fastq.gz\n";
		print TA "cat $fqlist2|gzip - >$outdir/Input_fq/2.Merge/$samp[0][1].L1_2.fastq.gz\n";
		#print TA "gzip   $outdir/Input_fq/2.Merge/$samp[0][1].L1_1.fastq\n";
                #print TA "gzip   $outdir/Input_fq/2.Merge/$samp[0][1].L1_2.fastq\n";
		if($random ne "NO"){
			print TA "perl $Bin/random_read.pl $outdir/Input_fq/2.Merge/$samp[0][1].L1_1.fastq.gz   $outdir/Input_fq/2.Merge/$samp[0][1].L1_2.fastq.gz $random\n";
			print TA "mv $outdir/Input_fq/2.Merge/$samp[0][1].L1_1.fastq.$random.gz $outdir/Input_fq/3.Random/\n";
			print TA "mv $outdir/Input_fq/2.Merge/$samp[0][1].L1_2.fastq.$random.gz $outdir/Input_fq/3.Random/\n";
			#print TA "gzip $outdir/Input_fq/3.Random/$samp[0][1].L1_1.fastq.$random\n";
			#print TA "gzip $outdir/Input_fq/3.Random/$samp[0][1].L1_2.fastq.$random\n";
		}
	}
	print TA "echo \"trim adapter of $samp[0][1] finish\" >Preprocess_$index.mark\n";
        print TA "echo ===end at : `date` ===\n";
        close TA;
}
