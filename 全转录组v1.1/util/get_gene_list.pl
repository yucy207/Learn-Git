

my ($gene_diff, $out_dir) = @ARGV;
my @temp = ();
open DE, $gene_diff or die "Can't open $gene_diff\n";

open UP,   qq{>$out_dir/up.list} or die "Can't open $out_dir/up.list\n";
open DOWN, qq{>$out_dir/down.list} or die "Can't open $out_dir/down.list\n";
open ALL,  qq{>$out_dir/de_genes.list} or die "Can't open $out_dir/de_genes.list\n";

my $cnt = 0;
while (<DE>) {
	chomp;
	$cnt++;
	my @arr = split;
	next if $cnt == 1;
	next if $arr[$#arr] =~ /DEG/;
	print UP   qq{$arr[$0]\n} if $arr[$#arr] eq 'Up';
	print DOWN qq{$arr[0]\n} if $arr[$#arr] eq 'Down';
	print ALL  qq{$arr[0]\n};
}
close DE;
close UP;
close DOWN;
close ALL;
