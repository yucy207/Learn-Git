#!/usr/bin/perl

use Parallel::ForkManager;


my ($chr_length, $bam, $out_dir, $strand) = @ARGV;

#my $chr_length = shift @ARGV;
#my $bam        = shift @ARGV;
#my $out_dir    = shift @ARGV;
#my $strand     = shift @ARGV;
#my $samtools   = join " ", @ARGV;

system qq{mkdir -p $out_dir} if not -d $out_dir;

my %chr = ();

my $MAX_PROCESS = 10;
my $pm = Parallel::ForkManager->new($MAX_PROCESS);

read_length();
submit();


sub read_length
{
	open EXO, $chr_length;
	while (<EXO>) {
		chomp;
		my @arr = split /\s+/;
		$chr{$arr[0]} = $arr[1];
	}
	close EXO;
}


sub cal_depth
{
	my $chr   = shift;
	my $start = shift;
	my $end   = shift;
	#print qq{$start\t$end\n};
	system qq{$samtools samtools depth $bam -r $chr:$start-$end > $out_dir/$chr.$start.depth\n};

	my %meta = ();
	open EXO, qq{$out_dir/$chr.$start.depth\n};
	while (<EXO>) {
		chomp;
		my @arr = split /\s+/;
		my $pos = int($arr[1] / 1000);
		#print qq{$pos\n};
		push @{$meta{$pos}}, $arr[2];
	}
	close EXO;
	
	my $start_index = int($start / 1000);
	my $end_index   = int($end / 1000);
	open SAVE, qq{>$out_dir/$chr.$start.report\n};
	foreach my $x ($start_index..$end_index) {
		my $pos = $x + 1;
		my $val = exists $meta{$pos} ? cal_median(@{$meta{$pos}}) : 0;
		print SAVE qq{$pos\t$val\t$chr\t$strand\n};
	}
	close SAVE;
}

sub cal_median{
	my @val = sort {$a <=> $b} @_;
	if ( $#val % 2 == 0){
		my $index = $#val - $#val / 2;
		return $val[$index];
	} else {
		my $left_index = ($#val + 1) / 2 - 1;
		my $right_index = ($#val + 1) / 2;
		my $res = ($val[$left_index] + $val[$right_index]) / 2;
		return $res; 
	}
}


sub submit
{
	my @task;
	my $step = 999999;
	foreach my $x (keys %chr) {
		my $start = 0;
		while ($start < $chr{$x}) {
			my $end = $start + $step;
			push @task, [$x, $start, $end];
			#cal_depth($x, $start, $start + $step);
			$start += ($step + 1);
		}
	}

	foreach my $x (@task) {
		my $chr = $x->[0];
		my $start = $x->[1];
		my $end = $x->[2];
		my $pid = $pm->start and next;
		cal_depth($chr, $start, $end);
		$pm->finish;
	}
	$pm->wait_all_children;
	print qq{Finish !\n};
}

