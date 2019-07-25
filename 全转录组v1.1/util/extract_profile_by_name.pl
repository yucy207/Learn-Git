#!/usr/bin/env perl

my ($count, $sample) = @ARGV;

my @samples = split /,/, $sample;

my $cnt   = 0;
my $id    = '';
my @names = ();
my @genes = ();
my %meta  = ();
open EXO, $count or die "Can't open $count!\n";
while (<EXO>) {
	chomp;
	$cnt++;
	my @arr    = split /\t/;
	if ($cnt == 1) {
		$id    = shift @arr;
		@names = @arr; 
	} else {
		push @genes, $arr[0];
		foreach my $x (1..$#arr) {
			$meta{$arr[0]}{$names[$x - 1]} = $arr[$x];

		}
	}

}
close EXO;


my $title = join "\t", @samples;
#print qq{$id\t$title\n};
foreach  my $x (@genes) {
	my @res = ();
	foreach my $y (@samples) {
		push @res, $meta{$x}{$y};
	}
	my $line = join "\t", @res;
	print qq{$x\t$line\n};
}
