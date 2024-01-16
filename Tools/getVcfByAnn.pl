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

my (@input , $help , $prx , $vcf);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"v|vcf=s"	=>	\$vcf,
	"help"	=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

my %s;
for my $file (@input){
	open IN , "$file";
	while (<IN>){
			chomp;
			my @F = split /\t/ , $_;
			if ($F[5] eq '-'){
				$F[2]--;
			}
			$s{$F[1]}->{$F[2]} = 1;
	}
	close IN;
}

open IN , "$vcf";
while (<IN>){
	if (/^#/){
		print $_;
		next;
	}
	my @F = split /\t/ , $_;
	print $_ if exists $s{$F[0]}->{$F[1]};
}
close IN;



sub help{
print << "EOF!";
#===============================================================================
#===============================================================================
EOF!
}



