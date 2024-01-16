#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use List::Util qw[min max];
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib $Bin;

my ($id, $outfile, $symmetric, $outIGV, $infile, $calls, $help , $minMarkN , $minLen);
$id=0;
$calls = 1;
$symmetric = 1;
$minMarkN = 5;
$minLen = 1000;
GetOptions (
		'id=s' => \$id,
	    'infile|i=s' => \$infile,
	    'symmetric|s=i' => \$symmetric,
        'outfile|o=s' => \$outfile,
        'outIGV|igv=s' => \$outIGV,
		'minMarkN=s'	=>	\$minMarkN,
		'minLen=s'	=>	\$minLen,
        'help|?' => \$help
);
if ($help or ! $infile){
	&help;
	exit;
}

print "Parameters:\nid=$id\ninfile=$infile\noutfile=$outfile\noutIGV=$outIGV\n";

my ($name,$path,$suffix) = fileparse($infile);
my $nameOut = $id;

open(SEGFILE, $infile) || die("Can't open $infile!\n");
my $header = <SEGFILE>;
my $line = <SEGFILE>; chomp($line);
my($chr, $start, $ref, $nRef, $N, $ratio, $logR, $cn, $state, $call, $clust, $clustFreq, @rest) = split(/\t/,$line);
my $end = $start;
my @totalRatio=();
my @totalLogR=();
push(@totalRatio,max($ratio,1-$ratio));
push(@totalLogR,$logR);
my @block;
while($line=<SEGFILE>){
	chomp($line);
	my ($chrS, $startS, $refS, $nRefS, $NS, $ratioS, $logRS, $cnS, $stateS, $callS, $clustS, $clustFreqS, @rest) = split(/\t/,$line);	
	if ($chrS ne $chr || $callS ne $call || $clustS ne $clust or $startS - $end >= 1000000){
		if (scalar(@totalRatio) >= $minMarkN and $end-$start >= $minLen){
			my $medianRatio = max($ratio,1-$ratio); #1 position only
			my $medianLogR = $logR;
			$medianRatio = median(@totalRatio) if (($end-$start+1)>1); #more than 1 position
			$medianLogR = median(@totalLogR) if (($end-$start+1)>1);
			my @majmin;
			if ($symmetric == 0){
				@majmin = &getMinorMajorCN($state);
			}elsif ($symmetric == 1){
				@majmin = &getMinorMajorCNSymmetric($state);
			}
			my $output = $nameOut . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . ($end-$start+1). "\t" . sprintf("%.2f",$medianRatio) . "\t" . sprintf("%.2f",$medianLogR);	
			$output = $output . "\t" .  $state . "\t" . $call if ($calls);
			$output = $output . "\t" . $cn . "\t" . $majmin[0] . "\t" . $majmin[1] . "\t" . $clust . "\t" . $clustFreq;
			my $output2 = $nameOut . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . scalar(@totalRatio). "\t" . sprintf("%.2f",$medianLogR);
		
			push @block , {
				'segs'	=>	$output,
				'igv'	=>	$output2,
				'ratio'	=>	[@totalRatio],
				'logr'	=>	[@totalLogR],
				'bb'	=>	[$chr, $start , $end, $ref, $nRef, $N, $ratio, $logR, $cn, $state, $call, $clust, $clustFreq]
			};
			#reset
		}else{
			if ($#block >= 0){
				my $pre = pop @block;
				($chr, $start , $end, $ref, $nRef, $N, $ratio, $logR, $cn, $state, $call, $clust, $clustFreq) = @{$pre->{bb}};
				@totalRatio = @{$pre->{ratio}};
				@totalLogR = @{$pre->{logr}};
				redo;
			}
		}
		@totalRatio = ();
		@totalLogR = ();
		$start = $startS; 
	}else{
		
		#push(@totalRatio,max($ratioS,1-$ratioS));
		#assign current state to previous variables
		#new end, but start still the same
	}
	($chr, $end, $ref, $nRef, $N, $ratio, $logR, $cn, $state, $call, $clust, $clustFreq) = ($chrS, $startS, $refS, $nRefS, $NS, $ratioS, $logRS, $cnS, $stateS, $callS, $clustS, $clustFreqS);
	push(@totalRatio,max($ratioS,1-$ratioS));
	push(@totalLogR,$logRS);
}

#final output
if (scalar(@totalRatio) >= $minMarkN and $end-$start >= $minLen){
	my $medianRatio = max($ratio,1-$ratio); #1 position only
	$medianRatio = median(@totalRatio) if (($end-$start+1)>1); #more than 1 position
	my $medianLogR = $logR;
	$medianLogR = median(@totalLogR) if (($end-$start+1)>1);
	my $output = $nameOut . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . ($end-$start+1). "\t" . sprintf("%.4f",$medianRatio) . "\t" . sprintf("%.2f",$medianLogR);	
	my @majmin;
	if ($symmetric == 0){
			@majmin = &getMinorMajorCN($state);
	}elsif ($symmetric == 1){
			@majmin = &getMinorMajorCNSymmetric($state);
	}
	$output = $output . "\t" .  $state . "\t" . $call if ($calls);
	$output = $output . "\t" . $cn . "\t" . $majmin[0] . "\t" . $majmin[1] . "\t" . $clust . "\t" . $clustFreq;
	my $output2 = $nameOut . "\t" . $chr . "\t" . $start . "\t" . $end . "\t" . scalar(@totalRatio). "\t" . sprintf("%.2f",$medianLogR);
	push @block , {
		'segs'	=>	$output,
		'igv'	=>	$output2,
		'ratio'	=>	[@totalRatio],
		'logr'	=>	[@totalLogR],
		'bb'	=>	[$chr, $end, $ref, $nRef, $N, $ratio, $logR, $cn, $state, $call, $clust, $clustFreq]
	};
}
open OUTFILE, ">$outfile" || die("Cannot write to $outfile\n");
open IGV, ">$outIGV" || die("Cannot write to $outIGV\n");
print OUTFILE "Sample\tChromosome\tStart_Position(bp)\tEnd_Position(bp)\tLength(bp)\tMedian_Ratio\tMedian_logR\tTITAN_state\tTITAN_call\tCopy_Number\tMinorCN\tMajorCN\tClonal_Cluster\tClonal_Frequency\n";
print IGV "sample\tchr\tstart\tend\tnum.mark\tmedian\n";
for my $bb (@block){
	print OUTFILE $bb->{segs} , "\n";
	print IGV $bb->{igv} , "\n";
}
close IGV;
close OUTFILE;


sub median { 
	my (@array_ref) = @_; 
	my $count = scalar @array_ref; 
	# Sort a COPY of the array, leaving the original untouched 
	my @array = sort { $a <=> $b } @array_ref; 
	if ($count % 2) { 
		return $array[int($count/2)]; 
	} else { 
		return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
	} 
} 


sub getMinorMajorCN {
	my ($state) = @_;
	my @majmin; $majmin[0]="NA"; $majmin[1]="NA";
	if ($state==0){
		$majmin[0] = 0; $majmin[1] = 0;
	}elsif ($state==1 || $state==2){
		$majmin[0] = 0; $majmin[1] = 1;
	}elsif($state==4){
		$majmin[0] = 1; $majmin[1] = 1;
	}elsif ($state==3 || $state==5){
		$majmin[0] = 0; $majmin[1] = 2;
	}elsif ($state==6 || $state==9){
		$majmin[0] = 0; $majmin[1] = 3;
	}elsif ($state==10 || $state==14){
		$majmin[0] = 0; $majmin[1] = 4;
	}elsif ($state==15 || $state==20){
		$majmin[0] = 0; $majmin[1] = 5;
	}elsif ($state==7 || $state==8){
		$majmin[0] = 1; $majmin[1] = 2;
	}elsif ($state==11 || $state==13){
		$majmin[0] = 1; $majmin[1] = 3;
	}elsif ($state==16 || $state==19){
		$majmin[0] = 1; $majmin[1] = 4;
	}elsif ($state==17 || $state==18){
		$majmin[0] = 2; $majmin[1] = 3;
	}elsif ($state==12){
		$majmin[0] = 2; $majmin[1] = 2;
	}
	#my $retVal = ($majmin[0],$majmin[1]);	
	return @majmin;
}

sub getMinorMajorCNSymmetric {
	my ($state) = @_;
	my @majmin; $majmin[0]="NA"; $majmin[1]="NA";
	if ($state==0){
		$majmin[0] = 0; $majmin[1] = 0;
	}elsif ($state==1){
		$majmin[0] = 0; $majmin[1] = 1;
	}elsif($state==2){
		$majmin[0] = 0; $majmin[1] = 2;
	}elsif ($state==3){
		$majmin[0] = 1; $majmin[1] = 1;
	}elsif ($state==4){
		$majmin[0] = 0; $majmin[1] = 3;
	}elsif ($state==5){
		$majmin[0] = 1; $majmin[1] = 2;
	}elsif ($state==6){
		$majmin[0] = 0; $majmin[1] = 4;
	}elsif ($state==7){
		$majmin[0] = 1; $majmin[1] = 3;
	}elsif ($state==8){
		$majmin[0] = 2; $majmin[1] = 2;
	}elsif ($state==9){
		$majmin[0] = 0; $majmin[1] = 5;
	}elsif ($state==10){
		$majmin[0] = 1; $majmin[1] = 4;
	}elsif ($state==11){
		$majmin[0] = 2; $majmin[1] = 3;
	}elsif ($state==12){
		$majmin[0] = 0; $majmin[1] = 6;
	}elsif ($state==13){
		$majmin[0] = 1; $majmin[1] = 5;
	}elsif ($state==14){
		$majmin[0] = 2; $majmin[1] = 4;
	}elsif ($state==15){
		$majmin[0] = 3; $majmin[1] = 3;
	}elsif ($state==16){
		$majmin[0] = 0; $majmin[1] = 7;
	}elsif ($state==17){
		$majmin[0] = 1; $majmin[1] = 6;
	}elsif ($state==18){
		$majmin[0] = 2; $majmin[1] = 5;
	}elsif ($state==19){
		$majmin[0] = 3; $majmin[1] = 4;
	}elsif ($state==20){
		$majmin[0] = 0; $majmin[1] = 8;
	}elsif ($state==21){
		$majmin[0] = 1; $majmin[1] = 7;
	}elsif ($state==22){
		$majmin[0] = 2; $majmin[1] = 6;
	}elsif ($state==23){
		$majmin[0] = 3; $majmin[1] = 5;
	}elsif ($state==24){
		$majmin[0] = 4; $majmin[1] = 4;
	}
	#my $retVal = ($majmin[0],$majmin[1]);	
	return @majmin;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: titanSegmentFiles.pl
#
#        USAGE: ./titanSegmentFiles.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wang Ruiru (wangrr), wangruiru\@berrygenomics.com
# ORGANIZATION: Berry Genomics
#      VERSION: 1.0
#      CREATED: 02/27/18 14:45:02
#     REVISION: ---
#===============================================================================
EOF!
}



