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
use Spreadsheet::WriteExcel;

my (@input , $help , $out);
GetOptions(
	"i|input=s{1,}"	=>	\@input,
	"o|output=s"	=>	\$out,
	"help"	=>	\$help,
);

if ($help or $#input<0){
	&help;
	exit;
}

my $workbook = Spreadsheet::WriteExcel->new("$out.xls");
my $format = $workbook->add_format();
$format->set_font('Times New Roman');

for my $file (@input){
	open IN , "$file";
	my $base = basename $file;
	$base =~ s/\..*$//;
	my $worksheet = $workbook->add_worksheet($base);
	my ($row) = (0);
	while (<IN>){
		chomp;
		my @F = split /\t/ , $_;
		for my $col (0..$#F){
			my $bb = Encode::decode("utf8" , $F[$col]);
			$worksheet->write($row , $col , $bb , $format); 
		}
		$row++;
	}
	close IN;
}

sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: tab2Excel.pl
#
#        USAGE: ./tab2Excel.pl  
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
#      CREATED: 11/16/17 10:48:24
#     REVISION: ---
#===============================================================================
EOF!
}



