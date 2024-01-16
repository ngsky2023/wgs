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

my %mateid;
open IN , "$input";
while (<IN>){
	if (/^#/ or $_ !~ /SVTYPE=BND;/){
		print $_;
		next;
	}

	chomp;
	my @F = split /\t/ , $_;
	my ($md) = (/MATEID=([^;]+);/);
	if (exists $mateid{$F[2]}){
		if ($mateid{$F[2]} == 1){
			next;
		}else{
			print "$_\n";
			next;
		}
	}

#	chr7    133785004       9:0     C       C[chr7:133798330[       0       PASS    SVTYPE=BND;MATEID=9:1;ASSEMBLED;EVENT=8 G
	if (/\w+\[(\w+):(\d+)\[/){
		my ($chr , $pos) = ($1 , $2);
		if ($chr eq $F[0] and $pos > $F[1]){
			$F[4] = '<DEL>';
			$F[7] =~ s/SVTYPE=BND;/SVTYPE=DEL;/;
			$F[7] = "END=$pos;SVLEN=".($pos-$F[1]).";$F[7]";
			$F[7] =~ s/MATEID=([^;]+);//;
			$mateid{$md} = 1;
#		}elsif ($chr eq $F[0] and $pos < $F[1]){
#			$F[4] = '<DUP>';
#			$F[7] =~ s/SVTYPE=BND;/SVTYPE=DUP;/;
#			$F[7] = "END=$pos;SVLEN=".($F[1]-$pos).";$F[7]";
#			$F[7] =~ s/MATEID=([^;]+);//;
#			$mateid{$md} = 1;
		}else{
			$mateid{$md} = 2;
		}
	}elsif (/\](\w+):(\d+)\]\w+/){
		my ($chr , $pos) = ($1 , $2);
		if ($chr eq $F[0] and $pos < $F[1]){
			$F[4] = '<DEL>';
			$F[7] =~ s/SVTYPE=BND;/SVTYPE=DEL;/;
			$F[7] = "END=$pos;SVLEN=".($F[1]-$pos).";$F[7]";
			$F[7] =~ s/MATEID=([^;]+);//;
			$mateid{$md} = 1;
#		}elsif ($chr eq $F[0] and $pos > $F[1]){
#			$F[4] = '<DUP>';
#			$F[7] =~ s/SVTYPE=BND;/SVTYPE=DUP;/;
#			$F[7] = "END=$pos;SVLEN=".($pos-$F[1]).";$F[7]";
#			$F[7] =~ s/MATEID=([^;]+);//;
#			$mateid{$md} = 1;
		}else{
			$mateid{$md} = 2;
		}
	}elsif (/\w+\](\w+):(\d+)\]/){
		my ($chr , $pos) = ($1 , $2);
		if ($chr eq $F[0] and $pos > $F[1]){
			$F[4] = '<INV>';
			$F[7] =~ s/SVTYPE=BND;/SVTYPE=INV;/;
			$F[7] = "END=$pos;SVLEN=".($pos-$F[1]).";$F[7]";
			$F[7] =~ s/MATEID=([^;]+);//;
			$mateid{$md} = 1;
		}else{
			$mateid{$md} = 2;
		}
	}elsif (/\[(\w+):(\d+)\[\w+/){
		my ($chr , $pos) = ($1 , $2);
		if ($chr eq $F[0] and $pos < $F[1]){
			$F[4] = '<INV>';
			$F[7] =~ s/SVTYPE=BND;/SVTYPE=INV;/;
			$F[7] = "END=$pos;SVLEN=".($F[1]-$pos).";$F[7]";
			$F[7] =~ s/MATEID=([^;]+);//;
			$mateid{$md} = 1;
		}else{
			$mateid{$md} = 2;
		}
	}else{
		$mateid{$md} = 2;
	}
	if ($mateid{$md} == 1){
		print join("\t" , @F) , "\n";
	}else{
		print "$_\n";
	}


}
close IN;

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: shiftBND.pl
#
#        USAGE: ./shiftBND.pl  
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
#      CREATED: 01/23/18 16:14:43
#     REVISION: ---
#===============================================================================
EOF!
}



