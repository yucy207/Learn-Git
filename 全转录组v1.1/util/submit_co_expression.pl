
use Parallel::ForkManager;
#my ($source_count, $target_count,  $out_dir, $source, $target, $rscript, $co_expression) = @ARGV;

my $source_count = shift @ARGV;
my $target_count = shift @ARGV;
my $out_dir      = shift @ARGV;
my $source       = shift @ARGV;
my $target       = shift @ARGV;
my $co_expression = shift @ARGV;
my $rscript = join " ", @ARGV;

system qq{mkdir -p $out_dir/split} if not -d qq{$out_dir/split};
system qq{mkdir -p $out_dir/correlation} if not -d qq{$out_dir/correlation};

split_count();

run();

sub split_count
{
	system qq{split -a 4 -l 50 $source_count $out_dir/split/x};
}

sub run
{
	my $max_threads = 50;
	my $pm = Parallel::ForkManager->new($max_threads);
	foreach my $x (`ls $out_dir/split`) {
		my $pid = $pm->start and next;
		chomp($x);
		system qq{$rscript $co_expression  $out_dir/split/$x $target_count $source $target $out_dir/correlation/$x.xls\n};
		$pm->finish;
	}
	$pm->wait_all_children;
	system qq{(echo -e "$source\\t$target\\tcor\\tpvalue";cat $out_dir/correlation/*.xls | grep -v "cor") > $out_dir/correlation.xls};
}
