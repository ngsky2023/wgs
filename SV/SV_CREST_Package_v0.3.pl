#!/usr/bin/perl -w
use strict;
use File::Basename;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use Getopt::Long;

my $crest_pl        = "perl $Bin/crest-4f6952e0af48/CREST.pl";
my $extractSClip_pl = "perl $Bin/crest-4f6952e0af48/extractSClip.pl";
my $env_sh          = "$Bin/CREST.Activation.sh";
my $faToTwoBit      = "$Bin/crest-4f6952e0af48/DependenceSoft/Blat_34/faToTwoBit";
my $gfServer        = "$Bin/crest-4f6952e0af48/DependenceSoft/Blat_34/gfServer";

my ($outdir, $shdir, $ref, $gff, $range, $port, $bamlist, $host, $Copt, $Eopt, $help, $Ehelp, $Chelp );
GetOptions(
	"h!"        => \$help,
	"outdir=s"  => \$outdir,
	"shdir=s"   => \$shdir,
	"ref=s"     => \$ref,
	"gff=s"     => \$gff,
	"gfSport=i" => \$port,
	"gfShost=s" => \$host,
	"bamlist=s" => \$bamlist,
	"range=s"   => \$range,
	"Copt=s"    => \$Copt,
	"Eopt=s"    => \$Eopt,
	"Ehelp!"    => \$Ehelp,
	"Chelp!"    => \$Chelp
);
die `source $env_sh ; $crest_pl -h`        if $Chelp;
die `source $env_sh ; $extractSClip_pl -h` if $Ehelp;
die &usage() if !$ref || !$gff || (!$bamlist );

$outdir ||= "./";
`mkdir -p $outdir` if ( !-d $outdir );
$outdir = abs_path($outdir);
$shdir ||= $outdir;
`mkdir -p $shdir` if ( !-d $shdir );
$Copt    ||= "";
$Eopt    ||= "";
$host    ||= "10.10.203.201";
$port    ||= int( rand(6000) ) + 9000;
$bamlist ||= "";

$ref = abs_path($ref);
my $ref_name = basename($ref);

my $tmp_range = $range ? "-r $range" : "";
my @bamfiles = "";
&getbamfiles($bamlist, \@bamfiles );

open SH, ">$shdir/crest_gfServer.sh" or die "Can not open file $shdir/crest_gfServer.sh $! .";
print SH "#!/bin/bash\n";
print SH "echo ===start at : `date` ===\n";
print SH "source $env_sh \n";
print SH "cd $outdir \n";
print SH "$faToTwoBit $ref $outdir/$ref_name.2bit \n";
print SH "$gfServer start $host $port -canStop $outdir/$ref_name.2bit >$outdir/gfServer.log 2>&1  \n\n";
print SH "echo \"crest_gfServer finish\" >crest_gfServer.mark\n";
print SH "echo ===end at : `date` ===\n";
close SH;

foreach (<@bamfiles>) {
	my $tmp_ele_bam = $_;
	my $bam_name = basename($tmp_ele_bam);
	my $prefix = $bam_name;
	$prefix =~ s/\.sort\.rmdup\.bam//;
	my $tmp_sample_outdir = "$outdir/$prefix";
	system("mkdir $tmp_sample_outdir") unless(-e $tmp_sample_outdir);

	open IN,">$tmp_sample_outdir/split_bam.sh";
	print IN "$Bin/split_bam -i $tmp_ele_bam -p $prefix -m 10000000\n";
	print IN "echo \"split_bam finish\" >split_bam.mark\n";
	close IN;
	
	open CREST,">$tmp_sample_outdir/sv_$prefix.sh";
	print CREST "perl $Bin/run_CREST.pl $tmp_sample_outdir $ref $outdir/$ref_name.2bit $host $port $gff\n";
	close CREST;
}
open SH, ">$shdir/stop_gfServer.sh" or die "$!";
print SH "$gfServer stop $host $port\n"; ##close gfserver
close SH;

######---------------Sub usage(), contains detail Help info-------------------######
sub usage() {
	print "Name: SV_CREST_Package.pl
Description: The script is to run structure variation use CREST software
Version: 0.2  Date: 2014-11-10
Connector: leiyoubing\@berrygenomics.com
Usage: perl $0 -ref <Reference_genome_fa > -bam <bam file>
	-ref	the ref genome fasta file
	-gff    the gff annotate file.
	-outdir	output directory, defult [./]
	-shdir	shell script output directory,defult same as [-outdir]
	-range	the range of positions need to be extracted. Format: chr1 or chr1:500-5000. defult all the .bam && ref data, []
	-gfSport	the part number of blat seaver, defult [int(rand(6000))+9000]
	-gfShost	the host ip address for blat seaver, defult [10.10.203.201]
	-bamlist	the list file which contains multiple bam file, each bam file is a line
	-Eopt	the options for extractSClip step,defult no special parameter, and you can see detail info by use -Ehelp,[]
	-Ehelp	print out the help info of extractSClip step
	-Copt	the options for CREST step,defult no special parameter, and you can see detail info by use -Chelp,[]
	-Chelp	print out the help info of CREST step
for exsample: 
	perl $0 -ref -gff -bamlist \n";
	exit;
}

######--------Sub getbamfiles(),get all the bam files for next process---------######
sub getbamfiles() {
	my ( $bamlist, $bamfiles ) = @_;
	if ( $bamlist eq "" ) {
		&usage();
	}
	if ( -s $bamlist ) {
		open BAMs, $bamlist or die "Can not open file $bamlist $! .";
		while (<BAMs>) {
			chomp;
			push( @$bamfiles, $_ );
		}
		close BAMs;
	}
}

