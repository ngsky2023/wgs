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

my %ss;
for my $file (@ARGV){
	open IN , "$file";
	while (<IN>){
		chomp;
		next if /^#/;
		my @F = split /\t/ , $_;
		my $strand;
		if ($F[7] =~ /DEL/){
			$strand = ["+" , "-"];
		}elsif ($F[7] =~ /DUP/){
			$strand = ["-" , "+"];
		}elsif ($F[7] =~ /INV/){
			if ($F[7] =~ /INV3;/){
				$strand = ["+" , "+"];
			}elsif ($F[7] =~ /INV5;/){
				$strand = ["-" , "-"];
			}
		}elsif ($F[7] =~ /BND/){
			if ($F[4] =~ /^\[/){
				$strand = ["-" , "-"];
			}elsif ($F[4] =~ /^\]/){
				$strand = ["-" , "+"];
			}elsif ($F[4] =~ /\]$/){
				$strand = ["+" , "+"];
			}elsif ($F[4] =~ /\[$/){
				$strand = ["+" , "-"];
			}
		}
		unless ($strand){
			print STDERR "$file\t$_\n";
		}else{
			$ss{$F[2]} = $strand;
		}
	}
	close IN;
}

open IN , "$input";
while (<IN>){
	chomp;
	my @F = split /\t/ , $_;
	if (exists $ss{$F[6]}){
		@F[8,9] = @{$ss{$F[6]}};
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
#         FILE: addSTRAND.pl
#
#        USAGE: ./addSTRAND.pl  
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
#      CREATED: 08/17/2020 01:35:38 PM
#     REVISION: ---
#===============================================================================
EOF!
}



