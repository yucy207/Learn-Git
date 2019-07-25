package package::merge_transcript;
use Parallel::ForkManager;

sub run
{
	my $metadata  = shift;
	my $base      = shift;
	my @samples   = split /,/, $metadata->{'samples'};
	my $organ     = qq{$metadata->{organ}};

	my $assembly  = qq{$metadata->{project}/assembly/result};
	my $merge     = qq{$metadata->{project}/merge_transcript};

	my $cufflinks = qq{$base->{cufflinks_bin}};
	my $ref_gtf   = qq{$base->{$organ}{genome_gtf}};

	if (not -d $assembly) {

		print "[Warning] : Merge transcript is failed, you must assemble firstly!\n";
		exit;

	}

	system qq{mkdir -p $merge}        if not -d $circrna;
	system qq{mkdir -p $merge/log}    if not -d qq{$merge/log};
	system qq{mkdir -p $merge/run}    if not -d qq{$merge/run};
	system qq{mkdir -p $merge/result} if not -d qq{$merge/result};

	if (-e qq{$merge/result/merged.gtf}) {

		print qq{合并转录本已经运行完成!\n};
		return 0;

	}  

	# write merge gtf list
	open SAVE, qq{>$merge/result/merge.list} or die "Can't open $merge/result/merge.list!\n";
	my $gtfs = join "\n", map { qq{$assembly/$_/transcripts.gtf} }  @samples; 
	
	print SAVE $gtfs;
	close SAVE;

	open SAVE, qq{>$merge/run/run.sh} or die "Can't open $merge/run/run.sh\n";
	my $merge_cmd = qq{$cufflinks cuffmerge -p 30 -g $ref_gtf -o $merge/result $merge/result/merge.list};
	my $finish    = qq{touch $merge/result/cuffmerge.finish};
	print SAVE qq{$merge_cmd\n$finish\n};
	close SAVE;

	system qq{bash $merge/run/run.sh &> $merge/log/run.log};

	print qq{合并转录本已经运行完成!\n};

	if (-e qq{$merge/result/merged.gtf}) {
		map  { system qq{rm $assembly/$_/transcripts.gtf} }   @samples;
	}

}

1;
