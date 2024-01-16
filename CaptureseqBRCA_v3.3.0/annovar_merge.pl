#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
#---------------------------------------------------------
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));

my $ver="1.0";
my %opts;
GetOptions(\%opts,"f1=s","f2=s","o=s","h");

#&help()if(defined $opts{h});
if(!defined($opts{f1}) || !defined($opts{f2}) || !defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver


	Usage:

		-f1           Samplename_snp.annovar.3mode.xls             <file>         must be given
		-f2           Samplename_indel.annovar.3mode.xls           <file>         must be given
 		-o            Samplename_annovar.xls                       <file>         must be given	
		-h            Help document

	Usage End.

	exit;
}
#------------------------------
print "Program Start Time:$Time_Start\n";
my $f1 = $opts{f1};
my $f2 = $opts{f2};
my $out = $opts{o};

open F1,"<$f1"  or die "couldn't creat $f1,$!\n";
open F2,"<$f2"  or die "couldn't creat $f2,$!\n";
open OUT,">$out"  or die "couldn't creat $out,$!\n";

while (<F1>){
    chomp;
    next if ($_ =~ /^#/);
    next if ($_ =~ /^\s*$/);
	print OUT "$_\n";
}
close F1;
while (<F2>){
    chomp;
    next if ($_ =~ /^#/);
    next if ($_ =~ /^\s*$/);
	next if ($_ =~ /^SampleName/);
	print OUT "$_\n";
}
close F2;

close OUT;


my $DONE_TIME=time();
my $WORKDONE_TIME = &sub_format_datetime(localtime($DONE_TIME));
print "Program Done Time:$WORKDONE_TIME\n";

#-------------------------------
sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}



sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub show_log()
{
	my ($txt) = @_ ;
	my $time = time();
	my $Time = &sub_format_datetime(localtime($time));
	print "$Time:\t$txt\n" ;
	return ($time) ;
}
