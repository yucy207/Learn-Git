use strict;
use warnings;

open F ,"1";
open O ,">primer.txt";
while (<F>){
    my $line1=$_;
    my $line2=<F>;
    my @l1=split "\t" ,$line1;
    my @l2=split "\t" ,$line2;
    $l1[0]=~ s/\_$// ;
    # $l1[2]=lc($l1[2]);
    # $l2[2]=lc($l2[2]);
    # print $#l1;
    # print $l2[2];
    print O "$l1[0]\_$l1[1]\t$l1[2]";
    print O "$l1[0]\_$l2[1]\t$l2[2]";
}
close F;
close O;