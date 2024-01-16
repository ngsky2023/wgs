#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
#---------------------------------------------------------
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));

my $ver="1.0";
my %opts;
GetOptions(\%opts,"f=s","v=s","o=s","h");

#&help()if(defined $opts{h});
if(!defined($opts{f}) || !defined($opts{v}) ||!defined($opts{o}) || defined($opts{h}))
{
	print <<"	Usage End.";
	Description:
		
		Version: $ver
		 writer: yangrutao


	Usage:

		-f          level.xls           <file>         must be given
		-v          avinput file        <file>         must be given
 		-o          level_1base.xls     <file>         must be given	
		-h          Help document

	Usage End.

	exit;
}
#------------------------------
print "Program Start Time:$Time_Start\n";
my $level = $opts{f};
my $avin = $opts{v};
my $level_1b = $opts{o};
my ($genesymb,$alleid,$name);
my ($hgvsc,$hgvsp);
my (%hs_pos,%hs_ref,%hs_alt);

open AV,"<$avin" or die $!;
while (<AV>){
	chomp;
	next if ($_ =~ /^\s*$/);
	my @dat1=split /\t/,$_;
	$hs_pos{$dat1[0]}{$dat1[1]}{$dat1[3]}{$dat1[4]}=$dat1[6];
	$hs_ref{$dat1[0]}{$dat1[1]}{$dat1[3]}{$dat1[4]}=$dat1[8];
	$hs_alt{$dat1[0]}{$dat1[1]}{$dat1[3]}{$dat1[4]}=$dat1[9];
	
}
close AV;
#print Dumper(%hs_ref);
#print Dumper(%hs_alt);
open IN,"<$level"  or die "couldn't creat $level,$!\n";
open OUT,">$level_1b"  or die "couldn't creat $level_1b,$!\n";

while (<IN>){
    chomp;
    if ($_ =~ /^Level/){
		print OUT "$_\n";
		next;
	};
    next if ($_ =~ /^\s*$/);
    my @dat = split /\t/,$_;
    my $chr = $dat[3];
	my $pos1 = $dat[4];
    my $pos2 = $dat[5];
	my $ref = $dat[6];
	my $alt = $dat[7];
	#my $delsize = length($ref) - length($alt);
    if (($ref eq "-" || $alt eq "-") && defined $hs_pos{$chr}{$pos1}{$ref}{$alt}){
		$dat[4] = $hs_pos{$chr}{$pos1}{$ref}{$alt};
		
        if (length($ref) >= 1 && $alt eq "-") {
            $dat[5] = $pos1 + length($ref) -1;
        }
        if (length($alt) >= 1 && $ref eq "-") {
            $dat[5] = $pos1 ;
        }
		$dat[6]=$hs_ref{$chr}{$pos1}{$ref}{$alt};
		$dat[7]=$hs_alt{$chr}{$pos1}{$ref}{$alt};
		#print "$dat[6]\t$dat[7]/n";
	}
	print OUT join ("\t",@dat),"\n";
}
close IN;

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

sub run_or_die()
{
	my ($cmd) = @_ ;
	&show_log($cmd);
	my $flag = system($cmd) ;
	if ($flag != 0){
		&show_log("Error: command fail: $cmd");
		exit(1);
	}
	&show_log("done.");
	return ;
}

