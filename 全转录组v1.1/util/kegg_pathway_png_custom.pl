

my ($db, $list,$out) = @ARGV;

my %hash = ();
open TXT, $db;
while (<TXT>) {
	chomp;
	my @arr = split /\t/;
	$hash{$arr[0]} = $arr[1];
}
close TXT;

my %res = ();
open TXT, $list;
while (<TXT>) {
	chomp;
	my @arr = split /\t/;
	$arr[2] = lc($arr[2]);
	#print qq{$arr[2]\n};
	my $color = $arr[2] eq 'up' ? 'red' : 'blue';
	push @{$res{$arr[0]}}, [$arr[1], $color];
	
	
}
close TXT;

foreach my $x (keys %res) {
	my @urls = ();
	foreach my $y (@{$res{$x}}) {
		my $gene = $y->[0];
		my $id   = $hash{$gene};
		my $col  = $y->[1];
		my $val  = qq{$id%09$col};
		push @urls, $val;
	}
	my $line = join "+", @urls;
	print qq{"http://www.genome.jp/kegg-bin/show_pathway?$x/$line"\n};
}


