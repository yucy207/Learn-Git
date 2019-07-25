use strict;
use warnings;

my %hash;
my $IN = shift @ARGV;
my $OUT = shift @ARGV;
open IN , $IN or die "can't open $IN!\n