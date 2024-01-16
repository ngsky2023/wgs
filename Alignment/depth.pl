#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use File::Basename;

#hg19:2897310462 hg18:2885270602 mm9:2725765481
#my $bed18 = "/home/zhanghk/software/CancerPipe/script/hg18/nonN_region";
#my $bed19 = "/home/zhanghk/software/CancerPipe/script/hg19/nonN_region";
		
my $usage =
"
Usage:

  Options:
  -i|input      <FILE>   Input bam-file.
  -o|out        <DIR>    Output file.
  -s|samtools   <STR>    The all path of samtools,default is '/home/leiyoubing/software/samtools-0.1.19/samtools'.
  -q            <INT>    base quality threshold, default: 0.
  -Q            <INT>    mapping quality threshold, default: 0.
  -r|ref        <FILE>   Reference,such as '/home/zhanghk/software/CancerPipe/hg19/hg19.fa'.
  -b            <BED>    list of positions or regions for reference,such as hg19:'/home/zhanghk/software/CancerPipe/script/hg19/nonN_region', hg18:'/home/zhanghk/software/CancerPipe/script/hg18/nonN_region'.
  -h|help                Help

For example:
        perl $0 -i sample.bam -o ./JMW_QC.xls -s /home/leiyoubing/software/samtools-0.1.19/samtools -q 0 -Q 0 -r /home/~/hg19/hg19.fa
";

my ($bam,$outfile,$samtools,$bed,$basethres,$mapQthres,$ref,$help);
GetOptions(
  "i|input=s"=>\$bam,
  "o|out=s"=>\$outfile,
  "s|samtools=s"=>\$samtools,
  "b|bed=s"=>\$bed,
  "q=i"=>\$basethres,
  "Q=i"=>\$mapQthres,
  "r|ref=s"=>\$ref,
  "h|help"=>\$help
);
$basethres ||= 0;
$mapQthres ||= 0;
$samtools ||= "/share/public/software/samtools-1.4/samtools";
my $R = "/share/public/software/R-3.2.5/lib64/R/bin/R";

die "$usage\n" if (!$bam || !$outfile || !$ref || $help);
my $key = (split /\./, basename $bam)[0];
my $outdir = dirname $outfile;

##read reference for get total-number of base.
my $total_chr = 0;
open IN,"$ref" || die "$!";
while(my $line = <IN>){
	chomp($line);
	next if ($line =~ /^>/);
	$line =~ s/[Nn]//g;
	$line =~ s/\s//g;
	my $len = length $line;
	$total_chr += $len;
}
close IN;

##read bam-file for get depth.
my %depth=();
if (defined $bed){
	open DEPTH,"$samtools depth -q $basethres -Q $mapQthres -b $bed $bam | " or die;
}else{
	open DEPTH,"$samtools depth -q $basethres -Q $mapQthres $bam | " or die;
}
while(<DEPTH>)
{
	chomp;
	my @arry = split /\s+/,$_;
	$depth{$arry[2]}+=1;
}
close DEPTH;

my @depth=sort {$a<=>$b} keys %depth;
my $maxCov=0;
my $Average_sequencing_depth=0;
my $Average_sequencing_depth4=0;
my $Average_sequencing_depth10=0;
my $Average_sequencing_depth20=0;
my $Average_sequencing_depth50=0;
my $Average_sequencing_depth100=0;
my $Average_sequencing_depth200=0;
my $Coverage=0;
my $Coverage4=0;
my $Coverage10=0;
my $Coverage20=0;
my $Coverage50=0;
my $Coverage100=0;
my $Coverage200=0;

my $Coverage_bases=0;
my $Coverage_bases_4=0;
my $Coverage_bases_10=0;
my $Coverage_bases_20=0;
my $Coverage_bases_50=0;
my $Coverage_bases_100=0;
my $Coverage_bases_200=0;

my $total_Coverage_bases=0;
my $total_Coverage_bases_4=0;
my $total_Coverage_bases_10=0;
my $total_Coverage_bases_20=0;
my $total_Coverage_bases_50=0;
my $total_Coverage_bases_100=0;
my $total_Coverage_bases_200=0;

open HIS,">$outdir/depth_frequency.txt" or die;
open CUM,">$outdir/cumu.txt" or die;
open OUT,">>$outfile" or die;
foreach my $depth1 (sort {$a<=>$b} keys %depth)
{
	next if($depth1==0);
	my $per=$depth{$depth1}/$total_chr;
	$total_Coverage_bases += $depth1*$depth{$depth1};
	$Coverage_bases += $depth{$depth1};

	if($depth1>=4)	
	{
		$total_Coverage_bases_4 += $depth1 * $depth{$depth1};
		$Coverage_bases_4 += $depth{$depth1};
	}
	if($depth1>=10)
	{
		$total_Coverage_bases_10 += $depth1 * $depth{$depth1};
		$Coverage_bases_10 += $depth{$depth1};
	}
	if($depth1>=20)
	{
		$total_Coverage_bases_20 += $depth1 * $depth{$depth1};
		$Coverage_bases_20 += $depth{$depth1};
	}
	if($depth1>=50)
	{
		$total_Coverage_bases_50 += $depth1 * $depth{$depth1};
		$Coverage_bases_50 += $depth{$depth1};
	}
	if($depth1>=100)
	{
		$total_Coverage_bases_100 += $depth1 * $depth{$depth1};
		$Coverage_bases_100 += $depth{$depth1};
	}
	if($depth1>=200)
	{
		$total_Coverage_bases_200 += $depth1 * $depth{$depth1};
		$Coverage_bases_200 += $depth{$depth1};
	}

	$maxCov=$per if($maxCov<$per);
	my $tmp=0;
	print HIS "$depth1\t$per\n";
	foreach my $depth2(@depth)
	{
		$tmp+=$depth{$depth2} if($depth2 >= $depth1); 
	}
	$tmp=$tmp/$total_chr;
	print CUM "$depth1\t$tmp\n";
}

$Average_sequencing_depth=$total_Coverage_bases/$total_chr;
$Coverage=$Coverage_bases/$total_chr;
$Average_sequencing_depth4=$total_Coverage_bases_4/$total_chr;
$Coverage4=$Coverage_bases_4/$total_chr;
$Average_sequencing_depth10=$total_Coverage_bases_10/$total_chr;
$Coverage10=$Coverage_bases_10/$total_chr;
$Average_sequencing_depth20=$total_Coverage_bases_20/$total_chr;
$Coverage20=$Coverage_bases_20/$total_chr;
$Average_sequencing_depth50=$total_Coverage_bases_50/$total_chr;
$Coverage50=$Coverage_bases_50/$total_chr;
$Average_sequencing_depth100=$total_Coverage_bases_100/$total_chr;
$Coverage100=$Coverage_bases_100/$total_chr;
$Average_sequencing_depth200=$total_Coverage_bases_200/$total_chr;
$Coverage200=$Coverage_bases_200/$total_chr;



print OUT "Average sequencing depth\t",sprintf("%.2f",$Average_sequencing_depth),"\n";
print OUT "Coverage\t",sprintf("%.2f%%",100*$Coverage),"\n";
print OUT "Coverage at least 4X\t",sprintf("%.2f%%",100*$Coverage4),"\n";
print OUT "Coverage at least 10X\t",sprintf("%.2f%%",100*$Coverage10),"\n";
print OUT "Coverage at least 20X\t",sprintf("%.2f%%",100*$Coverage20),"\n";
print OUT "Coverage at least 50X\t",sprintf("%.2f%%",100*$Coverage50),"\n";
print OUT "Coverage at least 100X\t",sprintf("%.2f%%",100*$Coverage100),"\n";
print OUT "Coverage at least 200X\t",sprintf("%.2f%%",100*$Coverage200),"\n";

close HIS;
close CUM;
close OUT;
if(1)
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
	if($Average_sequencing_depth<30)
	{
		$xlim=100;
		$xbin=20;
	}elsif($Average_sequencing_depth < 50)
	{
		$xlim=160;
		$xbin=20;
	}elsif($Average_sequencing_depth  < 120)
	{
		$xlim=250;
		$xbin=50;
	}else{
		$xlim=600;
		$xbin=100;
	}
	histPlot($outdir,"$outdir/depth_frequency.txt",$ylim,$ybin,$xlim,$xbin, $key);
	cumuPlot($outdir,"$outdir/cumu.txt",$xlim,$xbin, $key);
}

sub cumuPlot {
	my ($outdir, $dataFile, $xlim, $xbin, $sample) = @_;
	my $figFile = "$outdir/$sample\_cumuPlot.pdf";
	my $Rline=<<Rline;
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$dataFile")
	opar <- par()
	x <- rt\$V1[1:($xlim+1)]
	y <- 100*rt\$V2[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col=rgb(148,0,211,max=255),type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,100,by=20)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/$sample\_cumuPlot.png",width = 480, height = 360, type = "cairo")
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col=rgb(148,0,211,max=255),type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,100,by=20)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
	
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("$R CMD BATCH  $figFile.R");
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
	plot(x,y,col=rgb(255,51,153,max=255),type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/$sample\_histPlot.png",width = 480, height = 360, type = "cairo")
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col=rgb(255,51,153,max=255),type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("$R CMD BATCH $figFile.R");
}
