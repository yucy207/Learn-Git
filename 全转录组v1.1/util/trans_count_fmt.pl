#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($gtf, $count,  $help);

GetOptions(
	'gtf|g=s'     => \$gtf,
	'count|c=s'    => \$count,
	'help|h!'       => \$help
);

if ($help or not $gtf or not $count) {
	usage();
	exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl $0 -gtf gtf -count count.tsv
	-gtf   -g gtf file [required]
	-count -c count.tsv  count file [required]
	-help  -h  print help message
EOF

print $help;

}


my %hash = ();

pasre_gtf();

my $cnt = 0;

open COUNT,  $count or die "can't open $count!\n";
while (<COUNT>) {
	chomp;
	$cnt++;
	my @arr = split /\t/;
	my $id;
	my $trans = shift @arr;
	my $line = join "\t", @arr;
	if ($cnt == 1) {
		$id = "gene_name";
	} else {
		
		
		$id = $hash{$trans};
		
		
	}
	print qq{$trans\t$id\t$line\n};
}
close COUNT;

sub pasre_gtf
{

	open GTF, $gtf or die "Can't open $gtf!\n";
	while (<GTF>) {
		chomp;
		my @arr = split /\t/;
		my ($gene_id) = $_ =~ /gene_id\s+\"(.+?)\"/;
		#my ($gene_name) = $_ =~ /gene_name\s+\"(.+?)\"/;
		my ($trans_id) = $_ =~ /transcript_id\s+\"(.+?)\"/;
		$hash{$trans_id} = $gene_id;

	}
	close GTF;
}
