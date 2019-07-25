use strict;
use warnings;

open F ,"primer.txt";
open O ,">gc";
my %gene;
while (<F>){
    my $line=$_;
    my $no = <F>;
    my ($g,$tmp)=split "\t" ,$line,2;
    my ($fg)=$g=~ /^([^_]+)/;
    ($g)=$g=~ /^(.+)_[RF]/;
    $gene{$fg}.="$g,";
    # print O "Gene = TBX2 = TBX2_1";
}
close F;
map {print O "Gene = $_ = $gene{$_}\n"} keys %gene;
# map {print $_ ,"\n"}keys %gene;
close O;