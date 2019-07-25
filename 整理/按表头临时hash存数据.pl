use strict;
use warnings;

# my %hash;
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
    my %hash=();
    my @line=split "\t";
    # map {$hash{$heads[$_]}=exists $line[$_]?$line[$_]:""} 0..$#heads;
    map {$hash{$heads[$_]}=exists $line[$_]?$line[$_]:""} 0..$#heads;#按表头处理后更换为新表头
    # map {$hash{$heads[$_]}=$line[$_]} 0..$#line;
    map {print OUT "$_\t"} @heads;
    print OUT "\n";
    map {print OUT "$hash{$_}\t"} @heads;
    print OUT "\n";
}
close IN;
close OUT;

# # foreach my $title(sort keys %hash){
# # }

# my @arr=("a","b","","d");
# my $a=join "&&", @arr;
# my @b=split "&&", $a;
# print $a,"\n";
# print $#b;