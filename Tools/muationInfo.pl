#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my ($input , @bam , $prx);
GetOptions(
	"input=s"	=>	\$input,
	"prx=s"	=>	\$prx,
	"bam=s{1,}"	=>	\@bam,
);


print STDERR  "\nStart Time\t:[" , sub_format_datetime(localtime(time())) , "]\n\n";

#C52     1       267108  G       T       1/1     0,21    H_175   M
#C52     1       267659  G       -       -;A;T;G;C;      0;0;0;12;0;     S       M

my %fh;
open $fh{snp} , ">$prx.snp.txt";
open $fh{indel} , ">$prx.indel.txt";
open $fh{xxx} , ">$prx.xxx.txt";


open IN , "$input";
while (<IN>){
	chomp;
	next if /^SampleName/;
	my $ll = $_;
	my ($chr , $str , $end , $ref , $alt) = (split /\t/ , $_)[0..4];
	if ($chr !~ /chr/){
		$chr = "chr$chr";
	}
	my $loc = $str;
	my $type = '';
	$ref =~ s/-//;
	$alt =~ s/-//;
	my $rl = length($ref);
	my $al = length($alt);

	if ($str == $end and $rl == $al and $rl == 1){
		$type = 'snp';
	}elsif ($rl == 0){
		$type = 'ins';
	}elsif ($al == 0){
		$type = 'del';
	}else{
		$type = 'xxx';
	}

	my %dd;
	my $depth;
	for my $file (@bam){
		open BAM , "-|" , "samtools view -F 1024 $file $chr\:$str-$end";
		while (my $sam = <BAM>){
			my @sp = split /\t/ , $sam;
			my ($readn , $refn) = (-1 , $sp[3]-1);
			$depth++;
			while ($sp[5] =~ /(\d+)([IDMS])/g){
				my ($d , $w) = ($1 , $2);
				if ($w =~ /S/){
					$readn += $d;
				}elsif ($w =~ /I/){
					my $bseq = substr($sp[9] , $readn , $d);
					if ($type eq 'ins' and $bseq eq $alt){
						$dd{'ins'}++;
					}
					$readn += $d;
				}elsif ($w =~ /D/){
#					if ($refn+1<=$loc and $loc<=$refn+$d and $d == length($ref)){
					if ($refn+1==$loc and $d == length($ref)){
						if ($type eq 'del'){
							if ($d == $end-$str+1){
								$dd{'del'}++;
							}
						}
						last;
					}
					$refn += $d;
				}elsif ($w =~ /M/){
					if ($refn+1<=$loc and $loc<=$refn+$d){
						if ($type eq 'snp'){
							$dd{substr($sp[9] , $loc-$refn+$readn , 1)}++;
						}elsif ($type eq 'xxx'){
							if (substr($sp[9] , $str-$refn+$readn , $rl) eq $ref){
								$dd{'ref'}++;
							}elsif (substr($sp[9] , $str-$refn+$readn , $al) eq $alt){
								$dd{'alt'}++;
							}
						}
						last;
					}
					$refn += $d;
					$readn += $d;
				}else{
					print STDERR "$d$w\n";
				}
			}
		}
		close BAM;
	}
	my $fh;
	if ($type eq 'snp'){
		$fh = $fh{snp};
	}elsif ($type eq 'ins' or $type eq 'del'){
		$fh = $fh{indel};
	}else{
		$fh = $fh{xxx};
	}
	print {$fh} "$ll";
	if ($type eq 'snp'){
		my $altn = 0;
		my $depth = 0;
		for my $base ('A','T','G','C'){
			my $nn = 0;
			if (exists $dd{$base}){
#				print {$fh} "\t" , $dd{$base};
				$nn = $dd{$base};
			}else{
#				print {$fh} "\t0";
				$nn = 0;
			}
			$depth += $nn;
			if ($base eq $alt){
				$altn = $nn;
			}
		}
		print {$fh} "\t$depth";
        if ($depth==0){
                print {$fh} "\t$altn\t0";
        }else{
                print {$fh} "\t$altn\t" , $altn/$depth;
        }
	}elsif ($type eq 'del'){
		print {$fh}  "\t$depth\t";
		if (exists $dd{del}){
			print {$fh} $dd{del} , "\t" , $dd{del}/$depth;
		}else{
			print {$fh} "0\t0";
		}
	}elsif ($type eq 'ins'){
		print {$fh}  "\t$depth\t";
		if (exists $dd{ins}){
			print {$fh} $dd{ins} , "\t" , $dd{ins}/$depth;
		}else{
			print {$fh} "0\t0";
		}
	}elsif ($type eq 'xxx'){
		for my $base ('ref','alt'){
			if (exists $dd{$base}){
				print {$fh} "\t" , $dd{$base};
			}else{
				print {$fh} "\t0";
			}
		}
		if (exists $dd{alt}){
			print {$fh} "\t" , $dd{alt}/($dd{ref}+$dd{alt});
		}else{
			print {$fh} "\t0"
		}
	}
	print {$fh} "\n";
}
close IN;
close $fh{snp};
close $fh{indel};
close $fh{xxx};







print STDERR  "\nEnd Time\t:[" , sub_format_datetime(localtime(time())) , "]\n\n";

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


