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

if ($help or $#input<0){
	&help;
	exit;
}

my %m;
my %chr;
for my $file (@input){
	my $base = basename $file;
	my ($sample) = ($base =~ /^([^\.]+)/);
	open IN , "$file";
	<IN>;
	while (<IN>){
#		next unless /^HBV_gt[BC]/ or /^chr/;
		chomp;
		my @F = split /\t/ , $_;
		$m{$sample}->{$F[0]} = $F[1];
		$chr{$F[0]} = 1;
	}
	close IN;
}
my @chr = sort keys %chr;
my (%dchr , @zchr);
for my $c (@chr){
	if ($c =~ /^(\d+)/){
		$dchr{$c} = $1;
	}else{
		push @zchr , $c;
	}
}
@chr = sort {$dchr{$a}<=>$dchr{$b} or $a cmp $b} keys %dchr;
push @chr , @zchr;
print "sample\t" , join("\t" , @chr) , "\n";
for my $sample (sort keys %m){
	print $sample;
	for my $chr (@chr){
		print "\t" , $m{$sample}->{$chr};
	}
	print "\n";
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: mergeCount.pl
#
#        USAGE: ./mergeCount.pl  
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
#      CREATED: 04/10/18 09:13:20
#     REVISION: ---
#===============================================================================
EOF!
}



