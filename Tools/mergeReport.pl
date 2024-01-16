#!/usr/bin/env perl
use warnings;
use strict;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use lib $Bin;


my %loc;
my $pcgr = pop @ARGV;
for my $file (@ARGV){
	open IN , "$file";
	my $head = <IN>;
	chomp $head;
	print $head , "\t" , join("\t" , @h) , "\n";
	while (<IN>){
		chomp;
		my $out = "$_\t";
		my @F = split /\t/ , $_;
		$loc{"$F[1]:$F[2]:$F[3]"} = $_;
	}
	close IN;
}

my @col = (9,10,11,12,13,16,19,25,26,28,31,32);
my @h = qw/ONCOSCORE	ONCOGENE	TUMOR_SUPPRESSOR	CANCER_CENSUS_SOMATIC	CANCER_CENSUS_GERMLINE	PROTEIN_DOMAIN	CANCER_MUTATION_HOTSPOT	COSMIC_SITE_HISTOLOGY	COSMIC_DRUG_RESISTANCE	CLINVAR_SIG	TIER	TIER_DESCRIPTION/;
my $hn = @h;

if (-f $pcgr){
	open IN , "$pcgr";
	<IN>;
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		#g.chr17:7577121:G>A;
		my ($chr , $str , $ref , $var) = ($F[0] =~ /..(\w+):(\d+):(\w+)>(\w+)/);
		my ($lr , $lv) = (length($ref) , length($var));
		my $end;
		if ($lr == $lv){
			$end = $str + $lr - 1;
		}elsif ($lr > $lv){
			$end = $str++ + $lr - $lv;
		}else{
			$end = $str;
		}
		if (exists $loc{"$chr:$str:$end"}){
			print $loc{"$chr:$str:$end"} , "\t" , join("\t" , @F[@col]) , "\n";
		}
	}
	close IN;
}




