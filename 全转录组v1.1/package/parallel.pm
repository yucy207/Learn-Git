package package::parallel;

use Parallel::ForkManager;


sub run_sub
{
	my $max_threads = scalar @_;
	my $pm = Parallel::ForkManager->new($max_threads);
	foreach my $x (@_) {
		my $pid = $pm->start and next;
		&$x;
		$pm->finish;
	}
	$pm->wait_all_children;
}


sub run_command
{
	my $max_threads = scalar @_;
	my $pm = Parallel::ForkManager->new($max_threads);
	foreach my $x (@_) {
		my $pid = $pm->start and next;
		system $x;
		$pm->finish;
	}
	$pm->wait_all_children;	
}

1;
