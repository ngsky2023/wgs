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

my %p;
open IN , "$ARGV[0]";
while (<IN>){
	chomp;
	next if /^#/;
	my @F = split /\t/ , $_;
	#3,12,13,15
	my $tran = (split / / , $F[11])[0];
	next unless $tran;
	my @loc = split /\// , $F[14];
        if($loc[2] =~ /(.*?)\*fs\*\d+$/){
            $loc[2] = $1. "Xfs";
	}elsif ($loc[2] =~ /(.*)\*\d+$/){
            $loc[2] = $1;
        }elsif($loc[2] =~ /(.*)\*$/){
            $loc[2] = $1. "X";
        }

	$p{$F[2]}->{$tran} = [@loc[1,2]];
}
close IN;

open IN , "$input";
my $h = <IN>;
print $h;
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	$F[9] = $F[7] if $F[5] eq 'splicing';
	if (exists $p{$F[-8]}){
		my @ann = split /[,;]/ , $F[9];
		my $fix = '';
		if ($#ann == -1){
			print "$_\n";
			next;
		}
		for my $ann (@ann){
			my ($nm , $cds , $aa) = ('' , '' , '');
			for my $an (split /:/ , $ann){
				if ($an =~ /^NM/){
					$nm = $an;
					$fix .= "$an:";
				}elsif ($an =~ /^c\./){
					if (exists $p{$F[-8]}->{$nm}){
						$cds = $p{$F[-8]}->{$nm}->[0];
					}else{
						$cds = $an;
					}
					$fix .= "$cds:";
				}elsif ($an =~ /^p\./){
					if (exists $p{$F[-8]}->{$nm}){
						$aa = $p{$F[-8]}->{$nm}->[1];
					}else{
						$aa = $an;
					}
					$fix .= "$aa:";
				}else{
					$fix .= "$an:";
				}
			}
			if ($cds and not $aa and exists $p{$F[-8]}->{$nm}){
				$aa = $p{$F[-8]}->{$nm}->[1];
				$fix .= "$aa:";
			}	
			$fix =~ s/:$//;
			$fix .= ',';
		}
		$fix =~ s/,$//;
		$F[9] = $fix;
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
#         FILE: fixAnn.pl
#
#        USAGE: ./fixAnn.pl  
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
#      CREATED: 11/21/2018 11:31:08 AM
#     REVISION: ---
#===============================================================================
EOF!
}



