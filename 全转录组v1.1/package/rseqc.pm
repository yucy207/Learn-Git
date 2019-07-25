package package::rseqc;
use Parallel::ForkManager;
use Encode qw/decode/;
use package::format;
use package::parallel;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};

	my $organ    = qq{$metadata->{organ}};
	my $result   = qq{$metadata->{project}/RNA_evaluate};
	my $report   = qq{$metadata->{report}/02_Mapping_Evaluation};
	my $map      = qq{$metadata->{project}/mapping/result};

	my $rscript  = qq{$base->{rscript_bin}};
	my $util     = qq{$base->{'util'}};
	my $rseqc    = qq{$base->{rseqc_bin}};
	my $hisat2   = qq{$base->{hisat2_bin}};

	my $ref_bed12      = qq{$base->{$organ}{genome_bed12}};
	my $ref_chr_length = qq{$base->{$organ}{genome_chr_length}};

	if (not exists $metadata->{'organ'}) {

		print "[Warning] : You must set <organ> argument in the config.txt!\n";
		exit;
	}

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report/Mapping_Evaluation_Figures}  if not -d qq{$report/Mapping_Evaluation_Figures};
	
	my @res_samples = res_check(qq{$result/result}, \@samples);

	if (exists $base->{'force_sample'}) {
		foreach my $x (split /,/, $base->{'force_sample'}) {
			push @res_samples, $x;
		}
	}

	if (exists $base->{'force_step'}) {
		my @steps = split /,/, $base->{'force_step'};
		@res_samples = @samples if 5 ~~ @steps;
	}

	if (scalar @res_samples == 0) {
		
		print qq{rseqc 已经运行完成!\n};
		return;
	}

	pre_check($metadata, $base);
	
	my $max_threads = $base->{'thread_map'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_samples) {

		my $pid = $pm->start and next;
		system  qq{mkdir -p $result/result/$x}                      if not -d qq{$result/result/$x};
		system  qq{mkdir -p $report/Mapping_Evaluation_Figures/$x}  if not -d qq{$report/Mapping_Evaluation_Figures/$x};

		# rseqc
		my $mkdir         = qq{mkdir -p $result/result/$x};
		my $saturation    = qq{$rseqc RPKM_saturation.py  -i $map/$x/accepted_hits.bam -r $ref_bed12 -o $result/result/$x/$x -l 5 &>> $result/log/$x.log};
		my $plot          = qq{$rscript Rscript $util/sample_rpkm.R  $result/result/$x/$x.eRPKM.xls $report/Mapping_Evaluation_Figures/$x/$x.saturation.pdf &>> $result/log/$x.log};

		my $distribution  = qq{$rseqc read_distribution.py -i $map/$x/accepted_hits.bam -r $ref_bed12 &> $report/Mapping_Evaluation_Figures/$x/$x.distribution.xls };
		my $duplicate     = qq{$rseqc read_duplication.py  -i $map/$x/accepted_hits.bam -o $result/result/$x/$x &>> $result/log/$x.log};
		my $mv_plot       = qq{cp $result/result/$x/$x.DupRate_plot.pdf $report/Mapping_Evaluation_Figures/$x/$x.duplicate.rate.pdf &>> $result/log/$x.log};

		my $coverage      = qq{$rseqc geneBody_coverage.py -r $ref_bed12 -i $map/$x/accepted_hits.bam  -o $result/result/$x/$x &>> $result/log/$x.log};
		my $mv_coverage   = qq{cp $result/result/$x/$x.geneBodyCoverage.curves.pdf $report/Mapping_Evaluation_Figures/$x};

		# depth
		my $forward       = qq{$hisat2 samtools view -F 16 -h -b $map/$x/accepted_hits.bam > $result/result/$x/$x.forward.bam};
		my $forward_index = qq{$hisat2 samtools index $result/result/$x/$x.forward.bam};
		my $reverse       = qq{$hisat2 samtools view -f 16 -h -b $map/$x/accepted_hits.bam > $result/result/$x/$x.reverse.bam};
		my $reverse_inedx = qq{$hisat2 samtools index $result/result/$x/$x.reverse.bam};
		my $split_forward = qq{perl $util/split_depth.pl $ref_chr_length $result/result/$x/$x.forward.bam $result/result/$x/$x\_forward plus};
		my $split_reverse = qq{perl $util/split_depth.pl $ref_chr_length $result/result/$x/$x.reverse.bam $result/result/$x/$x\_reverse minus};

		my $cat           = qq{cat $result/result/$x/$x\_forward/*report  $result/result/$x/$x\_reverse/*report  > $result/result/$x/$x.depth};
		my $depth_plot    = qq{$rscript Rscript  $util/chr_depth.R  $result/result/$x/$x.depth  $report/Mapping_Evaluation_Figures/$x/$x.chr.depth &>> $result/log/$x.log \n};
		my $finish        = qq{touch $result/result/$x/$x.finish};

		open  SAVE, qq{>>$result/run/$x.sh} or die "Can't open $result/run/$x.sh\n";
		print SAVE qq{$mkdir\n$saturation\n$plot\n$distribution\n$duplicate\n$mv_plot\n$rseqc_sort\n$rseqc_index\n$coverage\n$mv_coverage\n};
		print SAVE qq{$finish\n};
		close SAVE;  

		my $saturation_proc = sub {
			system qq{$saturation\n$plot\n};
		};

		my $distribution_proc = sub {
			system qq{$distribution\n};
		};

		my $duplicate_proc = sub {
			system qq{$duplicate\n$mv_plot\n};
		};

		my $geneBodyCoverage_proc = sub {
			system qq{$coverage\n$mv_coverage\n};
		};

		my $depth_proc = sub {
			package::parallel::run_command(qq{$forward\n$forward_index\n$split_forward}, qq{$reverse\n$reverse_inedx\n$split_reverse});
			system qq{$cat\n$depth_plot\n$finish\n};

		};

		package::parallel::run_sub( $saturation_proc, $distribution_proc, $duplicate_proc, $geneBodyCoverage_proc);

		system qq{rm -r $result/result/$x/$x\_forward $result/result/$x/$x\_reverse $result/result/$x/$x.forward.bam $result/result/$x/$x.reverse.bam $result/result/$x/$x.depth} if -e qq{$report/Mapping_Evaluation_Figures/$x/$x.chr.depth.pdf};
		system qq{rm $result/result/$x.saturation.r} if -e qq{$result/result/$x.saturation.pdf};

		$pm->finish;
	} 

	$pm->wait_all_children;

	print qq{rseqc 运行完成!\n};

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my $map = qq{$metadata->{project}/mapping/result};
	# clean reads is exist
	foreach my $x (@samples) {

		my $bam = qq{$map/$x/accepted_hits.bam};
		die "$x bam file isn't exist!\n" if not -e qq{$bam};

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

