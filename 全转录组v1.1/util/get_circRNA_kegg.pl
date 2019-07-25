#!/usr/bin/env perl

my ($xls) = @ARGV;

open EXO, $xls;
while (<EXO>) {
	chomp;
	next if /^ID/;
	my @arr = split /\t/;
	next if $arr[4]  >0.05;
	foreach my $x (split /\//,$arr[7]) {
		print qq{$arr[0]\t$x\tup\n};
	} 
}
close EXO;

