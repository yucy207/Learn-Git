use strict;
use warnings;

my %hash;
my $IN = shift @ARGV;
my $OUT = shift @ARGV;
open IN , $IN or die "can't open $IN!\n";
open OUT, ">$OUT" or die "can't open $OUT!\n";
my $head= <IN>;
$head=~ s/[\r\n]//g;
my @heads= split "\t" , $head;
while (<IN>){
    $_=~ s/[\r\n]//g;
    # last if 
    # last if $.==3;
    # my %hash=();
    my @line=split "\t";
    map {$hash{$heads[$_]}{$.-1}=exists $line[$_]?$line[$_]:""} 0..$#heads;
    # map {$hash{$heads[$_]}=$line[$_]} 0..$#line;
    map {print OUT "$_\t"} @heads;
    print OUT "\n";
    map {print OUT "$hash{$_}{$.-1}\t"} @heads;
    print OUT "\n";
}
close IN;
close OUT;
print "$hash{'Position'}{'2'}";