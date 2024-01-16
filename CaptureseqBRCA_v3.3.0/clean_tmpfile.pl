#!/usr/bin/perl -w
use strict;
use File::Find;

if (@ARGV < 1){
	die "Usage:\n\tperl $0 project_analysis_path\n";
}

find(\&print_name,$ARGV[0]);

sub print_name
{
	my $file = $_;
	if($file =~ /.*txt$/){
		my $file_path = $File::Find::name;
		print "rm -rf $file_path\n";
		system("rm -rf $file_path");
	}elsif($file =~ /.*bam$/){
		unless($file =~ /.sort.rmdup.bam$/){
			my $file_path = $File::Find::name;
			print "rm -rf $file_path\n";
			system("rm -rf $file_path");
		}
	}elsif($file =~ /.*.sai$/){
		my $file_path = $File::Find::name;
		print "rm -rf $file_path\n";
		system("rm -rf $file_path");
	}elsif($file =~ /.*.bam.bai$/){
		unless($file =~ /.sort.rmdup.bam.bai$/){
			my $file_path = $File::Find::name;
			print "rm -rf $file_path\n";
			system("rm -rf $file_path");
		}
	}elsif($file =~ /.*.mark$/){
		my $file_path = $File::Find::name;
		print "rm -rf $file_path\n";
		system("rm -rf $file_path");
	}elsif($file =~ /.*\.sh\.o\d+$/){
		my $file_path = $File::Find::name;
		print "rm -rf $file_path\n";
		system("rm -rf $file_path");
	}elsif($file =~ /.*\.sh\.e\d+$/){
		my $file_path = $File::Find::name;
		print "rm -rf $file_path\n";
		system("rm -rf $file_path");
	}elsif($file =~ /.*.vcf$/){
		my $file_path = $File::Find::name;
		#print "rm -rf $file_path\n";
		#system("rm -rf $file_path");
	}
}

1;
