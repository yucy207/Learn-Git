
my ($gene_diff, $outdir) = @ARGV;

open DE, $gene_diff or die "Can't open $gene_diff\n";
open ALL, ">$outdir/de_genes.list" or die "Can't make $outdir/de_genes.list\n";

my %unique = ();
my $cnt = 0;
while (<DE>) {
	chomp;
	$cnt++;
	my @arr = split /\t/;
	next if $cnt == 1;
	next if $arr[$#arr] =~ /^Not/;
	next if $arr[1] eq 'n/a';
	my @ids = split /,/, $arr[1];
	foreach my $x (@ids) {
		next if exists $unique{$x};
		print ALL "$x\n";
		$unique{$x} = "-";
	}
}
