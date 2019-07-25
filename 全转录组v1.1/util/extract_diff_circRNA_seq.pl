#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($input, $fasta, $output,  $help);

GetOptions(
	'input|i=s'     => \$input,
	'fasta|f=s'     => \$fasta,
	'output|o=s'    => \$output,
	'help|h!'       => \$help
);

if ($help or not $input or not $output or not $fasta) {
	usage();
	exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl  -input input -fasta fasta -output output
	-input     -i  circRNA different expression file [required]
	-fasta     -f  circRNA spliced transcipt fasta [required]
	-output    -o  output file [required]
	-help      -h  print help message
EOF

print $help;

}


my %ids = prase_txt($input);

my %seqs = parse_fasta($fasta);

open SAVE, qq{>$output} or die "Can't open $output!\n";
foreach my $x (keys %ids) {
	next if not exists $seqs{$x};
	print SAVE qq{>$x\n$seqs{$x}\n};
}
close SAVE;


sub prase_txt
{
	my $txt  = shift;
	my %hash = ();

	open TXT, $txt or die "Can't open $txt!\n";
	while (<TXT>) {
		s/\s+$//g;
		my @arr = split /\t/;
		next if $arr[$#arr] eq 'type';
		next if $arr[$#arr] eq 'Not DEG';
		my @res = split /:/, $arr[0];
		$res[1] =~ s/-/|/;

		my $key = qq{$res[0]:$res[1]:$res[2]};
		#print qq{$key\n};
		$hash{$key} = $arr[1];


	}
	close TXT;

	return %hash;
}

sub parse_fasta
{
	my $fasta = shift;
	my %hash  = ();
	local $/  = ">";
	open FASTA, $fasta or die "Can't open $fasta!\n";
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($symbol, $seq) = split /\n/, $_, 2;
		my @arr    = split /\|/, $symbol;
		my $id     = qq{$arr[0]|$arr[1]};
		$seq =~ s/\s+//g;
		$hash{$id} = $seq;
	}
	close FASTA;
	return %hash;
}
