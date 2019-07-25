package package::assemble_isoform;
use Parallel::ForkManager;

sub run
{
	my $metadata  = shift;
	my $base      = shift;
	my @samples   = split /,/, $metadata->{'samples'};
	my $organ     = qq{$metadata->{organ}};

	my $map       = qq{$metadata->{project}/mapping/result};
	my $result    = qq{$metadata->{project}/assembly};

	my $stringtie = qq{$base->{stringtie_bin}};
	my $ref_gtf   = qq{$base->{$organ}{genome_gtf}};

	if (not -d $map){

		print "[Warning] : Assemble is failed, you must run mapping firstly!\n";
		exit;

	}

	my @res_samples = res_check(qq{$result/result}, \@samples);

	if (exists $base->{'force_sample'}) {
		foreach my $x (split /,/, $base->{'force_sample'}) {
			push @res_samples, $x;
		}
	}
	if (exists $base->{'force_step'}) {
		my @steps = split /,/, $base->{'force_step'};
		@res_samples = @samples if 6 ~~ @steps;
	}

	if (scalar @res_samples == 0) {
		print qq{转录本组装已经运行完成!\n};
		return;
	} 

	pre_check($metadata, $base);
	
	system qq{mkdir -p $result/run}       if not -d qq{$result/run};
	system qq{mkdir -p $result/result}    if not -d qq{$result/result};
	system qq{mkdir -p $result/log}       if not -d qq{$result/log};

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_samples) {

		my $pid    = $pm->start and next;

		system qq{mkdir -p $result/result/$x} if not -d qq{$result/result/$x};
		my $cmd    = qq{$stringtie stringtie $map/$x/accepted_hits.bam -o $result/result/$x/transcripts.gtf -p 20 -G $ref_gtf};
		my $finish = qq{touch $result/result/$x/$x.finish};

		open SAVE, qq{>$result/run/$x.sh} or die "Can't open $run/$x.cuffquant.sh\n";
		print SAVE qq{$cmd\n$finish\n};
		close SAVE;
		system qq{bash $result/run/$x.sh &>> $result/log/$x.log};

		$pm->finish;
	}

	$pm->wait_all_children;

	print qq{转录本组装已经运行完成!\n};

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};

	my $map      = qq{$metadata->{project}/mapping/result};
	# map bam is exist
	foreach my $x (@samples) {
		my $bam  = qq{$map/$x/accepted_hits.bam};
		die "$x bam file  isn't exist!\n" if not -e $bam;
	} 
	# reference genome gtf
	my $organ    = qq{$metadata->{organ}};
	my $ref_gtf  = qq{$base->{$organ}{genome_gtf}};
	die "reference genome gtf file isn't exist!\n" if not -e $ref_gtf;

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