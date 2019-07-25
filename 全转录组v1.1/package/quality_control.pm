package package::quality_control;
use Parallel::ForkManager;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my $adapter  = shift;
	   $adapter  = "truseq" if($adapter eq "");
	my @samples  = split /,/, $metadata->{'samples'};

	my $report   = qq{$metadata->{report}/01_Quality_Statistic/After_QC_Quality_Statistic};
	my $result   = qq{$metadata->{project}/trim};
	my $raw      = qq{$metadata->{raw_data}};

	my $fastp     = qq{$base->{fastp_bin}};
	my $adapter_1 = qq{$base->{adapter_1}{$adapter}};
	my $adapter_2 = qq{$base->{adapter_2}{$adapter}};
	die "adapter_1 isn't exist!\n" if not defined $adapter_1;
	die "adapter_2 isn't exist!\n" if not defined $adapter_2;

	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;
	system qq{mkdir -p $result/log}         if not -d qq{$result/log};
	system qq{mkdir -p $result/run}         if not -d qq{$result/run};
	system qq{mkdir -p $result/trim_galore} if not -d qq{$result/trim_galore};
	system qq{mkdir -p $result/fastx}       if not -d qq{$result/fastx};
	system qq{mkdir -p $result/fastqc}      if not -d qq{$result/fastqc};
	
	my @res_samples = res_check(qq{$result/trim_galore/}, \@samples);

	if (exists $base->{'force_sample'}) {
		foreach my $x (split /,/, $base->{'force_sample'}) {
			push @res_samples, $x;
		}
	}
	if (exists $base->{'force_step'}) {
		my @steps = split /,/, $base->{'force_step'};
		@res_samples = @samples if 2 ~~ @steps;
	}

	if (scalar @res_samples == 0) {

		print qq{质量控制已经运行完成!\n};
		return;		
	}

	pre_check($metadata, $base);

	my $max_threads = qq{$base->{thread_qc}};;
	my $pm = Parallel::ForkManager->new($max_threads);
	
	foreach my $x (@res_samples) {

		my $pid = $pm->start and next;

		system qq{mkdir -p $result/trim_galore/$x} if not -d qq{$result/trim_galore/$x};
		system qq{mkdir -p $result/fastqc/$x} if not -d qq{$result/fastqc/$x};

		open SAVE, qq{>$result/run/$x.trim.sh} or die "Can't open >$result/run/$x.trim.sh\n";

		my $trim_cmd  = qq{$fastp fastp -w 10 -l 40 -3 --in1 $raw/$x\_R1.fastq.gz --in2 $raw/$x\_R2.fastq.gz --out1 $result/trim_galore/$x/$x\_R1_val_1.fq.gz --out2 $result/trim_galore/$x/$x\_R2_val_2.fq.gz --adapter_sequence $adapter_1 --adapter_sequence_r2 $adapter_2 --html $result/trim_galore/$x/$x.html --json $result/trim_galore/$x/$x.json &> $result/trim_galore/$x/$x.log\n};
		my $finish    = qq{touch $result/trim_galore/$x/$x.finish};

		print SAVE qq{$trim_cmd\n};
		print SAVE qq{$finish\n};
		close SAVE;
		
		system qq{bash $result/run/$x.trim.sh > $result/log/$x.trim.log\n};

		$pm->finish;

	}
	
	$pm->wait_all_children;

	print qq{质量控制已经运行完成!\n};

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my $raw_data = qq{$metadata->{raw_data}/};
	# raw data is exist
	foreach my $x (@samples) {

		my $R1 = qq{$raw_data/$x\_R1.fastq.gz};
		my $R2 = qq{$raw_data/$x\_R2.fastq.gz};

		die "$x R1 fastq isn't exist!\n" if not -e $R1;
		die "$x R2 fastq isn't exist!\n" if not -e $R2;

	}
}

sub res_check
{
	my $output = shift;
	my $sample = shift;
	my @temp_sample = ();

	foreach my $x (@{$sample}) {
		next if -e qq{$output/$x/$x.finish};
		push @temp_sample, $x;
	}

	return @temp_sample;
}

1;
