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
my $T_depth = 8;
my $T_depth_Var = 2;
my $N_depth = 8;
my $N_depth_Var = 0;
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
	if (/^#/){
		@F[9,10] = ("NORMAL" , "TUMOR");
		print join("\t" , @F) , "\n";
		next;
	}
	next unless $F[6] eq 'PASS';
	next if $F[7] =~ /IMPRECISE/;
	if ($F[4] eq '<DUP:TANDEM>'){
		$F[4] = '<DUP>';
	}

		my @n = (split /:/ , $F[9]);
		my @t = (split /:/ , $F[10]);
		my ($nr , $nv , $tr , $tv) = (0 , 0 , 0 , 0);
		for my $c (@n){
			my ($r , $v) = split /,/ , $c;
			$nr += $r;
			$nv += $v;
		}
		for my $c (@t){
			my ($r , $v) = split /,/ , $c;
			$tr += $r;
			$tv += $v;
		}

		if ($nv > $tv){
			@F[9,10] = @F[10,9];
			($nr , $nv , $tr , $tv) = ($tr , $tv , $nr , $nv);
		}
		my ($nd , $td) = ($nr+$nv , $tr+$tv);
		next if $nd == 0 or $td == 0;
		my ($naf , $taf) = (sprintf("%.4f",$nv/$nd) , sprintf("%.4f",$tv/$td));
		next if $naf>$taf/5 or $tv < $T_depth_Var or $td < $T_depth or $nd < $N_depth or $nv > $N_depth_Var or $taf < $T_Var_Frac;
		$F[8] = "GT:AF:$F[8]";
		$F[9] = "0/0:$naf:$F[9]";
		$F[10] = "0/1:$taf:$F[10]";
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



