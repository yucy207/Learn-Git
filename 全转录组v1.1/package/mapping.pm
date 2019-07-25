package package::mapping;
use Parallel::ForkManager;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};

	my $trim       = qq{$metadata->{project}/trim/trim_galore};
	my $map        = qq{$metadata->{project}/mapping};
	my $organ      = qq{$metadata->{organ}};

	my $util          = qq{$base->{util}};
	my $hisat2        = qq{$base->{hisat2_bin}};
	my $hisat2_index  = qq{$base->{$organ}{hisat2_index}};
	my $samtools      = qq{$base->{samtools_bin}};
	my $sambamba      = qq{$base->{sambamba_bin}};

	if (not exists $metadata->{'organ'}) {

		print "[Warning] : You must set <organ> argument in the config.txt!\n";
		exit;
	}

	
	my @res_samples = res_check(qq{$map/result}, \@samples);

	if (exists $base->{'force_sample'}) {
		foreach my $x (split /,/, $base->{'force_sample'}) {
			push @res_samples, $x;
		}
	}
	if (exists $base->{'force_step'}) {
		my @steps = split /,/, $base->{'force_step'};
		@res_samples = @samples if 4 ~~ @steps;
	}

	if (scalar @res_samples == 0) {
		print qq{比对参考基因组已经运行完成!\n};
		return;
	}

	pre_check($metadata, $base);
	system qq{mkdir -p $map/run}    if not -d qq{$map/run};
	system qq{mkdir -p $map/result} if not -d qq{$map/result};
	system qq{mkdir -p $map/log}    if not -d qq{$map/log};

	my $max_threads = qq{$base->{thread_map}};
	my $pm = Parallel::ForkManager->new($max_threads);
	

	foreach my $x (@res_samples) {

		my $pid = $pm->start and next;
		system qq{mkdir -p $map/result/$x} if not -d qq{$map/result/$x};
		my $sam           = qq{$map/result/$x/$x.sam};
        my $mapping_stats = qq{$map/result/$x/align_summary.txt};

		my $mapping       = qq{$hisat2 -p 1 --dta --dta-cufflink -q -x $hisat2_index -1 $trim/$x/$x\_R1_val_1.fq.gz  -2 $trim/$x/$x\_R2_val_2.fq.gz --rg-id $x --rg SM:$x --rg LB:$x --rg PL:illumina --rg PU:$x 2> $mapping_stats | $samtools view -bh | $sambamba sort --tmpdir $map/result/$x  -t 1 -m 4G -o $map/result/$x/accepted_hits.bam /dev/stdin};
		my $finish        = qq{touch $map/result/$x/$x.finish};

		open  SAVE, qq{>$map/run/$x.sh} or die "Can't open $map/run/$x.sh\n";
		print SAVE qq{$mapping\n};
		print SAVE qq{$finish\n};
		close SAVE;

		system qq{bash $map/run/$x.sh &> $map/log/$x.log\n};

		$pm->finish;
	} 

	$pm->wait_all_children;
	print qq{比对参考基因组已经运行完成!\n};

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my $trim = qq{$metadata->{project}/trim/trim_galore};

	# clean reads is exist
	foreach my $x (@samples) {

		my $R1 = qq{$trim/$x/$x\_R1_val_1.fq.gz};
		my $R2 = qq{$trim/$x/$x\_R2_val_2.fq.gz};

		die "$x R1 clean fastq isn't exist!\n" if not -e $R1;
		die "$x R2 clean fastq isn't exist!\n" if not -e $R2;

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

