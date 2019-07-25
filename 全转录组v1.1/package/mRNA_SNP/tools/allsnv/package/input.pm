package package::input;
use strict;
use warnings;

sub run{
	my $input=<STDIN>;
	$input=~ s/[\r\n]//g;
	return $input;
}


1