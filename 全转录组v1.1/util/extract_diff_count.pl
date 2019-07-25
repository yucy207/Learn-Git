#!usr/bin/perl

use strict;
use warnings;

my($diff,$count,$out) = @ARGV;
my %gene = ();
open(DE,$diff) || die "cannot open diff.xls\n";
my @head = split/\t/,<DE> ;
my @samples = ();
foreach(@head){
	next if($_!~/\w/);
	push(@samples, $_) if($_!~/baseMean/);
	last if($_=~/baseMean/);
}
while(<DE>){
	my @line = split("\t",$_);
	next if $line[$#line] =~/Not DEG/;
	$gene{$line[0]} = "-";
}
close DE;

open(COUNT,$count) || die "cannot open count.rlog.xls\n";
open(OUT,">$out") || die "cannot make diff_count.rlog.xls\n";
my $a = <COUNT>;
$a =~s/[\r\n]//g;
my @title = split/\t/,$a;
my @index = ();
foreach my $sample(@samples){
	my ($index) = grep{$title[$_] eq $sample}0..$#title;
	push @index,$index;
}
my $title = join("\t",@samples);
print OUT "\t$title\n";
while(<COUNT>){
	$_ =~s/[\r\n]//g;
	my @arr = split("\t",$_);
	if(exists $gene{$arr[0]}){
		my @data;
		foreach my $x(@index){
			push(@data,$arr[$x]);
		}
		my $line = join("\t",$arr[0],@data);
		print OUT "$line\n";
	}
}
close COUNT;
close OUT;