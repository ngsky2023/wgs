#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin '$Bin';
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Path qw( mkpath );
use Math::CDF;
use lib $Bin;

print Math::CDF::pbinom($ARGV[0],$ARGV[1],$ARGV[2]) , "\n";
print Math::CDF::pnbinom($ARGV[0],$ARGV[1],$ARGV[2]) , "\n";


sub help{
print << "EOF!";
#===============================================================================
#
#         FILE: 1.pl
#
#        USAGE: ./1.pl  
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
#      CREATED: 11/22/2018 05:03:41 PM
#     REVISION: ---
#===============================================================================
EOF!
}



