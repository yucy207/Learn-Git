#!/usr/bin/env perl

my ($lncRNA_target,  $out_dir) = @ARGV;

system qq{mkdir -p $out_dir} if not -d $out_dir;

# extract target mRNA gene
my %mRNA_gene = ();


open EXO, $lncRNA_target or die  "Can't open $lncRNA_target!\n";
while (<EXO>) {
	chomp;
	next if /^Query/;
	my @arr    = split /\t/;
	next if $arr[5] + 0.1 > 0;
	my ($gene) = $arr[2] =~ /\((.+?)\)/;
	$mRNA_gene{$gene} = "-"; 
}
close EXO;




write_to_file(\%mRNA_gene,  qq{$out_dir/de_genes.list});

sub write_to_file
{
	my $hash = shift;
	my $out  = shift;
	open SAVE, qq{>$out} or die "Can't open $out!\n";
	foreach my $x (keys %{$hash}) {
		my $val = $hash->{$x};
		print SAVE qq{$x\n};
	}
	close SAVE;
}
