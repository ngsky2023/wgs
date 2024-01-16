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
GetOptions(
	"i|input=s"	=>	\$input,
	"help"	=>	\$help,
);

if ($help or ! $input){
	&help;
	exit;
}
my $T_depth =50;
my $T_depth_Var =15;
my $N_depth =20;
my $N_depth_Var =1;
my $T_Var_Frac =0.01;


if ($input =~ /\.gz/){
	open IN , "zcat $input|";
}else{
	open IN , "$input";
}
while (<IN>){
	if (/^##/){
		print $_;
		next;
	}
	chomp;
	my @F = split /\t/ , $_;
	if ($F[4] eq '<DUP:TANDEM>'){
		$F[4] = '<DUP>';
	}
	if (/^#/){
		$F[9] = "TUMOR";
	}else{
		my @n = (split /:/ , $F[9]);
		my ($nr , $nv , $tr , $tv) = (0 , 0);
		for my $c (@n){
			my ($r , $v) = split /,/ , $c;
			$nr += $r;
			$nv += $v;
		}
		next if $nv <= 2;
		my ($nd) = ($nr+$nv);
		next if $nd == 0;
		my ($naf) = (sprintf("%.4f",$nv/$nd));
		$F[8] = "GT:AF:$F[8]";
		$F[9] = "0/1:$naf:$F[9]";
	}
	print join("\t" , @F) , "\n";
}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#       AUTHOR: Wang Ruiru (wangrr), wangruiru\@berrygenomics.com
# ORGANIZATION: Berry Genomics
#      VERSION: 1.0
#      CREATED: 12/04/17 16:23:01
#     REVISION: ---
#===============================================================================
EOF!
}



