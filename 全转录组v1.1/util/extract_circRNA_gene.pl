
my ($gene_diff, $outdir) = @ARGV;

open DE, $gene_diff or die "Can't open $gene_diff\n";
open ALL, ">$outdir/de_genes.list" or die "Can't make $outdir/de_genes.list\n";
open UP, ">$outdir/up.list" or die "Can't make $outdir/up.list\n";
open DOWN, ">$outdir/down.list" or die "Can't make $outdir/down.list\n";

my %unique = ();
my $cnt = 0;
while (<DE>) {
	chomp;
	$cnt++;
	my @arr = split /\t/;
	next if $cnt == 1;
	next if $arr[$#arr] =~ /^Not/;
	my $index = $#arr - 8;
	next if $arr[$index] eq 'n/a';
	my @ids = split /,/, $arr[$index];
	foreach my $x (@ids) {
		next if exists $unique{$x};
		print ALL "$x\n";
		if($arr[$#arr] =~/Up/){
			print UP "$x\n";
		}else{
			print DOWN "$x\n";
		}
		$unique{$x} = "-";
	}
}
