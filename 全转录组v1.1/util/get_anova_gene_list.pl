#!usr/bin/perl
use strict;
use warnings;

my ($ano,$outpath) = @ARGV;
open(ANOVA,$ano) || die "cannot open ANOVA_Test_Result.xls\n";
open(OUT,">$outpath/de_genes.list") || die "cannot make de_genes.list\n";
my @title = split/\t/,<ANOVA>;
my ($index) = grep{$title[$_] eq 'Gene'} 0..$#title;
#print "$index";
my $num = 1;
while(<ANOVA>){
	$_=~s/[\r\n]//g;
	my @arr = split/\t/,$_;
	next if $arr[$#arr] =~/Not DEG/;
	if ($num <= 1000){
		print OUT "$arr[$index]\n";
	}
	$num ++;
}
close ANOVA;
close OUT;