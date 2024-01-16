#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib $Bin;

my ($input , $help);
my $T_depth =50;
my $T_depth_Var =4;
my $N_depth =20;
my $N_depth_Var =1;
my $T_Var_Strand =2;
my $T_Var_Frac =0.05;
my $somaticp =0.05;


GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}

my %pv = ('A'=>4,'C'=>5,'G'=>6,'T'=>7);
if ($input =~ /\.gz/){
	open IN , "zcat $input|";
}else{
	open IN , "$input";
}
while (<IN>){
	if (/^##/){
		print $_;
		next;
	}elsif (/^#/){
		print '#' , join("\t" , qw/CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT	TUMOR	NORMAL/) , "\n";
		next;
	}
	chomp;
	my @F = split /\t/ , $_;
	next unless $F[6] eq 'PASS';
	my @param=split(/;/,$F[7]);
	my $p=0;
	my $nt = '';
	my $sgt = '';
	while($p<=$#param){
		my @info=split(/[=,]/,$param[$p]);
		if ($info[0] eq 'NT'){
			$nt = $info[1];
		}elsif ($info[0] eq 'SGT'){
			$sgt = $info[1];
		}
		$p++;
	}
	next unless $nt eq 'ref';
	my $gt = '';
	if ($sgt =~ /ref->het/){
		$gt = '0/1';
	}else{
		$gt = '1/1';
	}
	#DP:FDP:SDP:SUBDP:AU:CU:GU:TU    60:2:0:0:0,0:0,0:0,0:58,60      294:13:0:0:0,0:2,2:11,16:268,279
	##DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50	22:22:22,23:0,1:0,16:19.25:2.20:0.00:0.11	51:51:36,41:9,15:12,88:51.07:9.64:0.00:0.18
	@F[9,10] = @F[10,9];
	my @n=split(/:/,$F[10]);
	my @t=split(/:/,$F[9]);
	my $nvar = $n[3];
	my $nref = $n[2];
	$nvar =~ s/,\d+//;
	$nref =~ s/,\d+//;
	my $tvar = $t[3];
	my $tref = $t[2];
	$tvar =~ s/,\d+//;
	$tref =~ s/,\d+//;
	my $ndp = $nvar + $nref;
	my $tdp = $tvar + $tref;
	next if $ndp <= 1 or $tdp <= 1;
	if ($tvar < $nvar){
		@F[9,10] = @F[10,9];
		($tvar, $nvar, $tdp, $ndp, $tref, $nref) = ($nvar, $tvar, $ndp, $tdp, $nref, $tref);
	}
	my $tvf = sprintf("%.4f" , $tvar/$tdp);
	my $nvf = sprintf("%.4f" , $nvar/$ndp);
	$F[8] = "GT:AD:AF:$F[8]";
	$F[9] = "$gt:$tref,$tvar:$tvf:$F[9]";
	$F[10] = "0/0:$nref,$nvar:$nvf:$F[10]";

	if($tdp>=$T_depth and $tvar>=$T_depth_Var and $ndp>=$N_depth and $nvar<=$N_depth_Var and $tvf>=$T_Var_Frac and $nvf<$tvf/2){
		print join("\t" , @F) , "\n";
	}
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: FilterStrelka.pl
#
#        USAGE: ./FilterStrelka.pl  
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
#      CREATED: 12/07/17 10:56:55
#     REVISION: ---
#===============================================================================
EOF!
}



