package package::data_eval;
use Parallel::ForkManager;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
			
	my $report   = qq{$metadata->{report}/01_Quality_Statistic/Raw_Reads_Quality_Statistic};
	my $result   = qq{$metadata->{project}/data_eval};
	my $raw      = qq{$metadata->{raw_data}};

	my $rscript  = qq{$base->{rscript_bin}};
	my $util     = qq{$base->{util}};
	my $fastqc   = qq{$base->{fastqc_bin}};
	my $fastx    = qq{$base->{fastx_bin}};
	my $seqtk    = qq{$base->{seqtk_bin}};

	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/fastqc} if not -d qq{$result/fastqc};
	system qq{mkdir -p $result/fastx}  if not -d qq{$result/fastx};

	my @res_samples = res_check(qq{$result/fastx}, \@samples);

	if (exists $base->{'force_sample'}) {
		foreach my $x (split /,/, $base->{'force_sample'}) {
			push @res_samples, $x;
		}
	}
	if (exists $base->{'force_step'}) {
		my @steps = split /,/, $base->{'force_step'};
		@res_samples = @samples if 1 ~~ @steps;
	}

	if (scalar @res_samples == 0) {
		print qq{质量评估已经运行完成!\n};
		return;
	}
	
	pre_check($metadata, $base);

	my $max_threads = qq{$base->{thread_qc}};
	my $pm = Parallel::ForkManager->new($max_threads);
	foreach my $x (@res_samples) {

		my $pid = $pm->start and next;

		system qq{mkdir -p $result/fastqc/$x} if not -d qq{$result/fastqc/$x};
		system qq{mkdir -p $report/$x}        if not -d qq{$report/$x};

		my $fastqc_cmd    = qq{$fastqc fastqc -t 50 -o $result/fastqc/$x --extract $raw/$x\_R1.fastq.gz $raw/$x\_R2.fastq.gz &>> $result/log/$x.log};
		my $fastx_r1_cmd  = qq{$seqtk fqchk $raw/$x\_R1.fastq.gz > $result/fastx/$x.R1.quality.stat.txt};
		my $fastx_r2_cmd  = qq{$seqtk fqchk $raw/$x\_R2.fastq.gz > $result/fastx/$x.R2.quality.stat.txt};
		my $cat           = qq{(grep -v "ALL\\|POS\\|min_len" $result/fastx/$x.R1.quality.stat.txt;grep -v "ALL\\|POS\\|min_len" $result/fastx/$x.R2.quality.stat.txt) > $result/fastx/$x.quality.stat};	
		my $base_com      = qq{$rscript Rscript $util/base_composition.R $result/fastx/$x.quality.stat $report/$x/$x.base.composition.pdf $x &>> $result/log/$x.log};
		my $error_rate    = qq{$rscript Rscript $util/error_rate.R $result/fastx/$x.quality.stat $report/$x/$x.error.rate.pdf $x &>> $result/log/$x.log};
		
		my $cp_r1         = qq{cp $result/fastqc/$x/$x\_R1_fastqc/Images/per_base_quality.png $report/$x/$x.R1.qual.png};
		my $cp_r2         = qq{cp $result/fastqc/$x/$x\_R2_fastqc/Images/per_base_quality.png $report/$x/$x.R2.qual.png};
		my $finish        = qq{touch $result/fastx/$x.finish};

		open SAVE, qq{>$result/run/$x.sh} or die "Can't open >$result/run/$x.sh\n";
		print SAVE qq{$fastqc_cmd\n$cp_r1\n$cp_r2\n};
		print SAVE qq{$fastx_r1_cmd\n$fastx_r2_cmd\n};
		print SAVE qq{$cat\n$base_com\n$error_rate\n};
		print SAVE qq{$finish\n};
		close SAVE;

		system qq{bash $result/run/$x.sh &> $result/log/$x.log\n};
				
		$pm->finish;

	}

	$pm->wait_all_children;
	
	print qq{质量评估已经运行完成!\n};

}


sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my $raw_data = qq{$metadata->{raw_data}};
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
		next if -e qq{$output/$x.finish};
		push @temp_sample, $x;
	}

	return @temp_sample;
}

1;
