my ($miranda) = @ARGV;
my $unique = ();
print qq{Seq1\tSeq2\tTot Score\tTot Energy\tMax Score\tMax Energy\tStrand\tLen1\tLen2\tPositions\n};
my  $flag = 0;
open EXO, $miranda or die "Can't open $miranda\n";
while (<EXO>) {
	chomp;
	$flag = 1 if /Score for this Scan:/;
	$flag = 0 if /^Complete/;
	next if $flag == 0;
	if (/>>/) {
		s/^>>//;
		my @arr = split /\s+/, $_, 10;
		next if exists $unique{$arr[0]}{$arr[1]};
		$arr[$#arr] =~ s/\s+/\|/g;
		$unique{$arr[0]}{$arr[1]} = "-";
		my @result = @arr;
		$result[0] = $arr[1];
		$result[1] = $arr[0];
                $result[7] = $arr[8];
                $result[8] = $arr[7];
		my $res = join "\t", @result;
		print qq{$res\n};
	}

}
close EXO;
