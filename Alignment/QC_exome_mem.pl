#!/usr/bin/env perl
use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my ($bamfile,$outdir,$regionfile,$Plot,$clean_data,$samtools,$bedtools,$type,$help);
GetOptions(
		"i:s"=>\$bamfile,
		"o:s"=>\$outdir,
		"r:s"=>\$regionfile,
		"c=s"=>\$clean_data,
		"s:s"=>\$samtools,
		"b:s"=>\$bedtools,
		"t:s"=>\$type,
		"plot"=>\$Plot,
		"h"=>\$help,
);

my $usage=<<USAGE;
usage:perl $0
		-i <bamfile>
		-r <region file>
		-c <clean_data>
		-o <outdir>
		-s <samtools path>
		-b <bedtools path>
		-t <'PE' or 'SE'>
		-plot
		-h help
USAGE

my %soft;
open IN , "$Bin/soft.list";
while (<IN>){
	chomp;
	my ($k , $d) = split /\t/ , $_;
	$soft{$k} = $d;
}
close IN;

$samtools ||= $soft{samtools};
$bedtools ||= $soft{bedtools};
$type ||= 'PE';
die $usage if (!$bamfile || $help || !$outdir || !$regionfile || not -f $clean_data );

mkdir ("$outdir", 0755) unless (-d $outdir);

my $sample_name = basename $outdir;
my $Initial_bases_on_target=0;
my $Initial_bases_near_target=0;
my $Initial_bases_on_or_near_target=0;
my $Total_effective_yield=0;
my $Effective_sequences_on_target=0;
my $Effective_sequences_near_target=0;
my $Uniquely_mapped_to_target=0;
my $Effective_sequences_on_or_near_target=0;
my $Uniquely_mapped_to_genome=0;
my $Fraction_of_effective_bases_on_target=0;
my $Fraction_of_uniq_mapping_on_target=0;
my $Average_sequencing_depth_on_target=0;
my $Average_sequencing_depth_near_target=0;
my $Fraction_of_effective_bases_on_or_near_target=0;
my $Mismatch_rate_in_target_region=0;
my $Mismatch_rate_in_all_effective_sequence=0;
my $Base_covered_on_target=0;
my $Coverage_of_target_region=0;
my $Base_covered_near_target=0;
my $Coverage_of_flanking_region=0;
my $Fraction_of_target_covered_with_at_least_200x=0;
my $Fraction_of_target_covered_with_at_least_150x=0;
my $Fraction_of_target_covered_with_at_least_100x=0;
my $Fraction_of_target_covered_with_at_least_50x=0;
my $Fraction_of_target_covered_with_at_least_20x=0;
my $Fraction_of_target_covered_with_at_least_10x=0;
my $Fraction_of_target_covered_with_at_least_4x=0;
my $Mismatch_base_in_target_region=0;
my $Mismatch_base_in_all_effective_sequence=0;
my $clean_reads=0;
my $mapped_reads=0;
my $mapping_rate=0;
my $dup_reads=0;
my $unique_reads=0;
my $unique_rate=0;
my $dup_rate=0;
my %hash;

open IN , "$clean_data";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	if ($F[0] eq 'total_reads'){
		$clean_reads = $F[1];
	}
}
close IN;

#get unique reads number of uniquely mapped to genome
my $key = (split /\./, basename $bamfile)[0];
open BAM,"$samtools view -q 10 -F 0x0004 $bamfile | " or die $!;
while(my $info = <BAM>){
	chomp($info);
	my @info=split /\s+/,$info;
	unless($info[1] & 0x800){ #0x800   SUPPLEMENTARY
		unless($info[1] & 0x100){ #0x100   SECONDARY
			unless($info[1] & 0x4){ #0x4     UNMAP
				$mapped_reads++;
				if($info[1] & 0x400){
					$dup_reads++;
					next;
				}
				if($info[4] >= 5)
				{
					$Uniquely_mapped_to_genome++;
				}
				if($info=~/NM:i:(\d+)/)
				{
					my $temp_m = $1;
					$temp_m -= $1 if($info[5] =~ /(\d+)[ID]/);
					$Mismatch_base_in_all_effective_sequence+=$temp_m;
				}
			}
		}
	}
}
close BAM;

$unique_reads = $mapped_reads - $dup_reads;
$unique_rate = sprintf("%.2f%%", 100*$unique_reads/$mapped_reads);

$mapping_rate=$mapped_reads/$clean_reads;
$mapping_rate=sprintf("%.2f%%",100*$mapping_rate);
$dup_rate=$dup_reads/$mapped_reads;
$dup_rate=sprintf("%.2f%%",100*$dup_rate);

open REG,"$regionfile" or die $!;
open FREG,">$outdir/freg_tmp.txt" or die $!;
while(<REG>)
{
	chomp;
	next if(/^\s*$/ || /^#/);
	my @info=split (/\t/,$_);
	next if (@info < 3);
	$hash{$info[0]}+=$info[2]-$info[1];
	$Initial_bases_on_target+=$info[2]-$info[1];
	my $pri_beg_pos=$info[1]-200;
	my $pri_end_pos=$info[1];
	my $next_beg_pos=$info[2];
	my $next_end_pos=$info[2]+200;

	$pri_beg_pos=0 if($pri_beg_pos < 0);
	if ($info[1] <= 1)
	{
		print FREG "$info[0]\t$next_beg_pos\t$next_end_pos\n";
	}
	else
	{
		print FREG "$info[0]\t$pri_beg_pos\t$pri_end_pos\n";
		print FREG "$info[0]\t$next_beg_pos\t$next_end_pos\n";
	}
}
close(FREG);
close(REG);

`for i in \$(seq 1 22) X Y;do grep -w chr\$i $outdir/freg_tmp.txt|sort -k2n >>$outdir/tmp.txt;done`;
`mv $outdir/tmp.txt $outdir/freg_tmp.txt`;
`$bedtools/mergeBed -i $outdir/freg_tmp.txt >$outdir/tmp.txt`;
`$bedtools/subtractBed -a $outdir/tmp.txt -b $regionfile >$outdir/freg.txt`;
`rm $outdir/freg_tmp.txt $outdir/tmp.txt`;

$Initial_bases_near_target=`awk '{total+=\$3-\$2};END{print total}' $outdir/freg.txt`;
chomp($Initial_bases_near_target);

`$samtools depth -d 50000 -Q 10 -b $regionfile $bamfile >$outdir/target.depth`;
`$samtools depth -d 50000 -Q 10 -b $outdir/freg.txt $bamfile >$outdir/flanking.depth`;

$Effective_sequences_on_target=`awk '{total+=\$3};\$3 != 0 {cov+=1};END{print total"\t"cov}' $outdir/target.depth`;
chomp($Effective_sequences_on_target);
($Effective_sequences_on_target,$Base_covered_on_target)=split(/\t/,$Effective_sequences_on_target);

$Effective_sequences_near_target=`awk '{total+=\$3};\$3 != 0 {cov+=1};END{print total"\t"cov}' $outdir/flanking.depth`;
chomp($Effective_sequences_near_target);
if($Effective_sequences_near_target!~/\d/){$Effective_sequences_near_target="0\t0";}
($Effective_sequences_near_target,$Base_covered_near_target)=split (/\t/,$Effective_sequences_near_target);

$Initial_bases_on_or_near_target=$Initial_bases_on_target+$Initial_bases_near_target;

`$samtools depth -Q 10 -d 50000 $bamfile >$outdir/whole_genome.depth`;
$Total_effective_yield=`awk '{total+=\$3};END{print total}' $outdir/whole_genome.depth`;
chomp $Total_effective_yield;


$Fraction_of_effective_bases_on_target=$Effective_sequences_on_target/$Total_effective_yield;
$Effective_sequences_on_or_near_target=$Effective_sequences_on_target + $Effective_sequences_near_target;
$Fraction_of_effective_bases_on_or_near_target=$Effective_sequences_on_or_near_target/$Total_effective_yield;

open TMP,"$samtools view -q 10 -F 0x4 -L $regionfile $bamfile | " or die $!;
while(<TMP>)
{
	chomp;
	my @arry = split /\t/, $_;
	my $fflag=$arry[1];
	unless($fflag & 0x800){ #0x800   SUPPLEMENTARY
		unless($fflag & 0x100){ #0x100   SECONDARY
			unless($fflag & 0x400){ #0x400     DUPLICATE
				if($arry[4] >= 5)
				{
					$Uniquely_mapped_to_target++;
				}
				my $temp_m = 0;
				if($_=~/NM:i:(\d+)/){
					$temp_m = $1;
					$temp_m -= $1 if($arry[5] =~ /(\d+)[ID]/);
				}
				$Mismatch_base_in_target_region+=$temp_m;
			}
		}
	}
}
close(TMP);

$Mismatch_rate_in_target_region=$Mismatch_base_in_target_region/$Effective_sequences_on_target;
$Mismatch_rate_in_all_effective_sequence=$Mismatch_base_in_all_effective_sequence/$Total_effective_yield;
$Fraction_of_uniq_mapping_on_target = $Uniquely_mapped_to_target/$Uniquely_mapped_to_genome;

$Average_sequencing_depth_on_target=$Effective_sequences_on_target/$Initial_bases_on_target;
$Average_sequencing_depth_near_target=$Effective_sequences_near_target/$Initial_bases_near_target;
	
$Coverage_of_target_region=$Base_covered_on_target/$Initial_bases_on_target;
$Coverage_of_flanking_region=$Base_covered_near_target/$Initial_bases_near_target;

my $tmp1=`awk '\$3 >=200 {total7++};\$3 >=150 {total6++};\$3 >=100 {total5++};\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};\$3 >=50 {total4++};END{print total1"\t"total2"\t"total3"\t"total4"\t"total5"\t"total6"\t"total7}' $outdir/target.depth`;
chomp($tmp1);
my @info1;
@info1=split /\t/,$tmp1;
$info1[0]=0 unless($info1[0] =~ /\d+/);
$info1[1]=0 unless($info1[1] =~ /\d+/);
$info1[2]=0 unless($info1[2] =~ /\d+/);
$info1[3]=0 unless($info1[3] =~ /\d+/);
$info1[4]=0 unless($info1[4] =~ /\d+/);
$info1[5]=0 unless($info1[5] =~ /\d+/);
$info1[6]=0 unless($info1[6] =~ /\d+/);
$Fraction_of_target_covered_with_at_least_200x=$info1[6]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_150x=$info1[5]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_100x=$info1[4]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_50x=$info1[3]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_20x=$info1[0]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_10x=$info1[1]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_4x=$info1[2]/$Initial_bases_on_target;

my $name=basename($bamfile);
$name=~s/\.bam//g;
my $sample=$name;
$sample =~ s/.H1.rmdup//;

`$samtools view -q 10 -f66 $bamfile |cut -f 9|sed 's/^-//' >$outdir/$sample.insert.txt\n"`;
Plot_insert_size($outdir,"$outdir/$sample.insert.txt",$sample);
my $mean_insert_size=`$samtools view -q 10 -f66 $bamfile |awk '{if (\$9 > 0 && \$9 < 500) {S+=\$9; T+=1}}END{print  S/T}' `;


open STAT,">$outdir/$key\_QC.xls" or die $!;
print STAT "Sample\t$sample\n";
print STAT "Clean reads\t$clean_reads\n";
print STAT "Mapped reads\t$mapped_reads($mapping_rate)\n";
print STAT "Duplicate reads\t$dup_reads($dup_rate)\n";
print STAT "Unique reads\t$unique_reads($unique_rate)\n";
print STAT "Reads uniquely mapped to target\t$Uniquely_mapped_to_target\n";
print STAT "Reads uniquely mapped to genome\t$Uniquely_mapped_to_genome\n";
print STAT "Total bases on target\t$Initial_bases_on_target\n";
print STAT "Total bases near target\t$Initial_bases_near_target\n";
printf STAT "Total sequences on target(Mb)\t%.2f\n",$Effective_sequences_on_target/1000000;
printf STAT "Total sequences near target(Mb)\t%.2f\n",$Effective_sequences_near_target/1000000;
printf STAT "Total effective yield(Mb)\t%.2f\n",$Total_effective_yield/1000000;
printf STAT "Fraction of effective bases on target\t%.2f%%\n",100*$Fraction_of_effective_bases_on_target;
printf STAT "Fraction of effective bases on or near target\t%.2f%%\n",100*$Fraction_of_effective_bases_on_or_near_target;
printf STAT "Fraction of uniquely mapped on target\t%.2f%%\n",100*$Fraction_of_uniq_mapping_on_target;
#printf STAT "Average sequencing depth on target\t%.2f\n",$Average_sequencing_depth_on_target;
printf STAT "Average sequencing depth on target\t%d\n",int($Average_sequencing_depth_on_target);
printf STAT "Average sequencing depth near target\t%.2f\n",$Average_sequencing_depth_near_target;
printf STAT "Average insert size of the library\t%.2f\n",$mean_insert_size;

print STAT "Base covered on target\t$Base_covered_on_target\n";
printf STAT "Coverage of target region\t%.2f%%\n",100*$Coverage_of_target_region;
print STAT "Base covered near target\t$Base_covered_near_target\n";
printf STAT "Coverage of near target\t%.1f%%\n",100*$Coverage_of_flanking_region;
printf STAT "Mismatch rate in target\t%.2f%%\n",100*$Mismatch_rate_in_target_region;
printf STAT "Mismatch rate in total\t%.2f%%\n",100*$Mismatch_rate_in_all_effective_sequence;
printf STAT "Fraction of target covered at least 4x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_4x;
printf STAT "Fraction of target covered at least 10x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_10x;
printf STAT "Fraction of target covered at least 20x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_20x;
printf STAT "Fraction of target covered at least 50x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_50x;
printf STAT "Fraction of target covered at least 100x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_100x;
printf STAT "Fraction of target covered at least 150x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_150x;
printf STAT "Fraction of target covered at least 200x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_200x;
close STAT;

open DB,"<$outdir/target.depth" or die $!;
open DF,">$outdir/depth_frequency.xls" or die $!;
my %depth=();
while(<DB>)
{
	chomp;
	my @tmp=split;
	$depth{$tmp[2]}++;
}
close(DB);
my $maxCov=0;
#$depth{"0"}+=$Initial_bases_on_target-$Base_covered_on_target;
	
foreach my $depth (sort {$a<=>$b} keys %depth)
{
	next if($depth==0);
	my $per=$depth{$depth}/$Initial_bases_on_target;
	$maxCov = $per if($per > $maxCov);
	print DF "$depth\t$per\t$depth{$depth}\n";
}
close(DF);
	
open CU,">$outdir/cumu.xls" or die $!;
print CU "Depth\tTRPercent\n";
my @depth= sort {$a<=>$b} keys %depth;
foreach my $depth1 (sort {$a<=>$b} keys %depth)
{
	my $tmp=0;
	next if($depth1==0);
	foreach my $depth2 (@depth)
	{
		if($depth2 >= $depth1)
		{
			$tmp+=$depth{$depth2};
		}
	}
	$tmp = $tmp/$Initial_bases_on_target;
	print CU "$depth1\t$tmp\n";
}
close(CU);

if($Plot)
{
	my $ylim = 100*$maxCov;
	my ($xbin,$ybin);
	$ylim= int($ylim) + 1;
	if($ylim <= 3)
	{
		$ybin = 0.5;
	}else{
		$ybin=1;
	}
	my $xlim=0;
	if($Average_sequencing_depth_on_target<30)
	{
		$xlim=100;
		$xbin=20;
	}elsif($Average_sequencing_depth_on_target < 50)
	{
		$xlim=160;
		$xbin=20;
	}elsif($Average_sequencing_depth_on_target < 120)
	{
		$xlim=250;
		$xbin=50;
	}else{
		$xlim=600;
		$xbin=100;
	}
	histPlot($outdir,"$outdir/depth_frequency.xls",$ylim,$ybin,$xlim,$xbin, $key);
	cumuPlot($outdir,"$outdir/cumu.xls",$xlim,$xbin, $key);
}

my %dep=();
my %cov=();
open IN,"$outdir/target.depth" or die "Can't open the $outdir/target.depth:$!";
while (<IN>)
{
	chomp;
	my ($chr,$dd)=(split /\t/,$_)[0,2];
	$cov{$chr}++;
	$dep{$chr}+=$dd;
}
close IN;
open OUT,">$outdir/chrall.stat" or die $!;
print OUT "Chr\tCoverage\tDepth\n";
foreach my $ch (sort keys %hash)
{
	my $cp = 100*$cov{$ch}/$hash{$ch};
	my $dp = $dep{$ch}/$hash{$ch};
	printf OUT "$ch\t%.2f%%\t%.2f\n",$cp,$dp;
}
close OUT;

`rm $outdir/freg.txt $outdir/whole_genome.depth $outdir/target.depth $outdir/flanking.depth`;

sub cumuPlot {
        my ($outdir, $dataFile, $xlim, $xbin, $sample) = @_;
        my $figFile = "$outdir/$sample\_cumuPlot.pdf";
        my $Rline=<<Rline;
        pdf(file="$figFile",w=8,h=6)
        rt <- read.table("$dataFile",header=T)
        opar <- par()
        x <- rt\$Depth[1:($xlim+1)]
        y <- 100*rt\$TRPercent[1:($xlim+1)]
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
        par(opar)
        dev.off()
        png(filename="$outdir/$sample\_cumuPlot.png",width = 480, height = 360)
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
        par(opar)
        dev.off()

Rline
        open (ROUT,">$figFile.R");
        print ROUT $Rline;
        close(ROUT);

        system("Rscript $figFile.R");
}


sub histPlot {
	my ($outdir, $dataFile, $ylim, $ybin, $xlim, $xbin, $sample) = @_;
	my $figFile = "$outdir/$sample\_histPlot.pdf";
	my $Rline=<<Rline; 
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$dataFile")
	opar <- par()
	t=sum(rt\$V2[($xlim+1):length(rt\$V2)])
	y=c(rt\$V2[1:$xlim],t)
	y <- y*100
	x <- rt\$V1[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/$sample\_histPlot.png",width = 480, height = 360)
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("/share/public/software/R-3.0.1/lib64/R/bin/Rscript $figFile.R");
	#system("rm  $figFile.R  $figFile.Rout");
}

sub Plot_insert_size {
	my ($outdir, $dataFile,$sample) = @_;
	my $figFile = "$outdir/$sample\_insert_size.pdf";
	my $Rline=<<Rline; 
	data <- scan("$dataFile")
	pdf("$figFile",width=8)
	if(range(data)[2] > 1000){
        	plot(table(data), xlab="Insert Size", ylab="Reads Number",xlim=c(0,1000))
        	data <- data[data<1000]
    		hist(data, xlim=c(0,1000))
	}else{
        	plot(table(data), xlab="Insert Size", ylab="Reads Number")
        	data <- data[data<1000]
		hist(data)
	}
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);
	system("/share/public/software/R-3.0.1/lib64/R/bin/Rscript $figFile.R");
}

