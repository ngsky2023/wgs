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

my (@input , $help);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"help"	=>	\$help,
);

if ($help or $#input < 0){
	&help;
	exit;
}

print "sample\tgene\ttranscript\texon\tc.hgvs\tp.hgvs\tZygosity\tsignificance\n";
for my $file (@input){
	open IN , "$file";
	<IN>;
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		next if $F[12] =~ /intron/ or $F[11] eq '' or $F[12] eq '.';
		next if $F[20] =~ /enign/;
		print "$F[0]\t$F[7]\t$F[11]\t$F[12]\t$F[13]\t$F[15]\t$F[16]\t$F[20]\t$F[2]\t$F[18]\n";
	}
	close IN;
}

sub help{
print << "EOF!";
#===============================================================================
#===============================================================================
EOF!
}



