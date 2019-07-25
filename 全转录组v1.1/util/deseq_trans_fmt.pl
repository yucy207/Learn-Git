#!/usr/bin/env perl

my ($gtf, $txt) = @ARGV;


my %hash = ();
open GTF, $gtf or die "Can't open $gtf!\n";
while (<GTF>) {
	chomp;
	my @arr = split /\t/;
	my ($trans) = $_ =~ /transcript_id\s+(.+?);/;
	my ($gene)  = $_ =~ /gene_id\s+(.+?);/;
	$trans      =~ s/"//g;
	$gene       =~ s/"//g;

	$hash{$trans} = $gene;

}
close GTF;


my $cnt = 0;
open TXT, $txt or die "Can't open $txt!\n";
while (<TXT>) {
	chomp;
	$cnt++;
	my @arr = split /\t/;
	my $id = shift @arr;
	my $line = join "\t", @arr;
	my $val = exists $hash{$id} ? $hash{$id} : 'gene_id';
	print qq{$id\t$val\t$line\n};
}
close TXT;