my %loc;
open IN , "$ARGV[0]";
while (<IN>){
	next if /^#/;
	my @F = split /\t/ , $_;
	$loc{"$F[0]\t$F[1]\t$F[3]\t$F[4]"} = 1;
}
close IN;

my $number = 0;
open IN , "$ARGV[1]";
while (<IN>){
	next if /^#/;
	my @F = split /\t/ , $_;
	if (exists $loc{"$F[0]\t$F[1]\t$F[3]\t$F[4]"}){
		$number++;
	}
}
close IN;

print "$ARGV[1]\t$number\n";
