package package::mRNA::mRNA_rmats;

use Parallel::ForkManager;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my @groups   = @{$metadata->{'group'}};
	my $organ    = qq{$metadata->{organ}};

	my $map      = qq{$metadata->{project}/mapping/result};
	my $result   = qq{$metadata->{project}/mRNA/rMATs};
	my $report   = qq{$metadata->{report}/03_mRNA_Analysis/04_Isoform_Alternative_Splicing};

	my $util     = qq{$base->{util}};
	my $rmats    = qq{$base->{rmats4_bin}};
	my $rmats2sashimiplot = qq{$base->{rmats2sashimiplot_bin}};
	my $ref_mRNA_gtf      = qq{$base->{$organ}{genome_mRNA_gtf}};

	if (not exists $metadata->{'group'}) {

		print "[Warning] : You must set group in the config.txt!\n";
		exit;
	}

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	
	my @res_groups = res_check(qq{$result/run}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{mRNA 可变剪切分析已经完成!\n};
		return 0;
	} 

	pre_check($metadata, $base);

	my $max_threads = $base->{'max_threads'};
	my $pm  = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $pid = $pm->start and next;

		my $control = $x->[0];
		my $case    = $x->[1];
		my $sample  = $x->[2];
		my @control_samples = split /,/, (split /;/, $sample)[0];
		my @case_samples    = split /,/, (split /;/, $sample)[1];
		my $case_sample     = join ",", map { qq{$map/$_/accepted_hits.bam}} @control_samples;
		my $control_sample  = join ",", map { qq{$map/$_/accepted_hits.bam}} @case_samples;

		my $out_dir = qq{$result/result/$case\_vs_$control};
		system qq{mkdir -p $out_dir} if not -d $out_dir;

		open CONTROL, qq{>$out_dir/b1.txt} or die "Can't open $out_dir/b1.txt!\n";
		print CONTROL qq{$case_sample\n};
		close CONTROL;

		open CASE, qq{>$out_dir/b2.txt} or die "Can't open $out_dir/b2.txt!\n";
		print CASE qq{$control_sample\n};
		close CASE;

		my $cmd        = qq{$rmats python /software/rMATS.4.0.1/rMATS-turbo-Linux-UCS2/rmats.py --b1 $out_dir/b1.txt --b2 $out_dir/b2.txt --od $out_dir --nthread 10 --gtf $ref_mRNA_gtf --tstat 10 --readLength 100 -t paired\n};
		my $report_out = scalar @groups > 1 ? qq{$report/$case\_vs_$control} : qq{$report};
		system qq{mkdir -p $report_out} if not -d $report_out;

		my $add   = qq{bash $util/add_exon_number.4.sh $ref_mRNA_gtf $out_dir/ $out_dir\n};
		my $fmt   = qq{perl $util/fmt_MATS.pl $out_dir $report_out\n};
		my $top10 = qq{perl $util/rMATsPlotsFiles.pl $organ $out_dir\n};
		
		# mkdir
		system qq{mkdir -p $out_dir/A3SS} if not -d qq{$out_dir/A3SS};
		system qq{mkdir -p $out_dir/A5SS} if not -d qq{$out_dir/A5SS};
		system qq{mkdir -p $out_dir/MXE}  if not -d qq{$out_dir/MXE};
		system qq{mkdir -p $out_dir/RI}   if not -d qq{$out_dir/RI};
		system qq{mkdir -p $out_dir/SE}   if not -d qq{$out_dir/SE};
		
		system qq{mkdir -p $report_out/Custom_Figures}      if not -d qq{$report_out/Custom_Figures};
		system qq{mkdir -p $report_out/Custom_Figures/A3SS} if not -d qq{$report_out/Custom_Figures/A3SS};
		system qq{mkdir -p $report_out/Custom_Figures/A5SS} if not -d qq{$report_out/Custom_Figures/A5SS};
		system qq{mkdir -p $report_out/Custom_Figures/MXE}  if not -d qq{$report_out/Custom_Figures/MXE};
		system qq{mkdir -p $report_out/Custom_Figures/RI}   if not -d qq{$report_out/Custom_Figures/RI};
		system qq{mkdir -p $report_out/Custom_Figures/SE}   if not -d qq{$report_out/Custom_Figures/SE};
		
		#rmats2sashimiplot
		my $A3SSplot = qq{$rmats2sashimiplot bash -c "rmats2sashimiplot --b1 $case_sample --b2 $control_sample --exon_s 1 --intron_s 5 -t A3SS -e $out_dir/A3SS.overlap.top10.txt -o $out_dir/A3SS --l1 $control --l2 $case --min-counts 0;chmod -R 775 $out_dir/A3SS/Sashimi_plot;cp $out_dir/A3SS/Sashimi_plot/* $report_out/Custom_Figures/A3SS/"\n};
		my $A5SSplot = qq{$rmats2sashimiplot bash -c "rmats2sashimiplot --b1 $case_sample --b2 $control_sample --exon_s 1 --intron_s 5 -t A5SS -e $out_dir/A5SS.overlap.top10.txt -o $out_dir/A5SS --l1 $control --l2 $case --min-counts 0;chmod -R 775 $out_dir/A5SS/Sashimi_plot;cp $out_dir/A5SS/Sashimi_plot/* $report_out/Custom_Figures/A5SS/"\n};
		my $MXEplot  = qq{$rmats2sashimiplot bash -c "rmats2sashimiplot --b1 $case_sample --b2 $control_sample --exon_s 1 --intron_s 5 -t MXE -e $out_dir/MXE.overlap.top10.txt -o $out_dir/MXE --l1 $control --l2 $case --min-counts 0;chmod -R 775 $out_dir/MXE/Sashimi_plot;cp $out_dir/MXE/Sashimi_plot/* $report_out/Custom_Figures/MXE/"\n};
		my $RIplot   = qq{$rmats2sashimiplot bash -c "rmats2sashimiplot --b1 $case_sample --b2 $control_sample --exon_s 1 --intron_s 5 -t RI -e $out_dir/RI.overlap.top10.txt -o $out_dir/RI --l1 $control --l2 $case --min-counts 0;chmod -R 775 $out_dir/RI/Sashimi_plot;cp $out_dir/RI/Sashimi_plot/* $report_out/Custom_Figures/RI/"\n};
		my $SEplot   = qq{$rmats2sashimiplot bash -c "rmats2sashimiplot --b1 $case_sample --b2 $control_sample --exon_s 1 --intron_s 5 -t SE -e $out_dir/SE.overlap.top10.txt -o $out_dir/SE --l1 $control --l2 $case --min-counts 0;chmod -R 775 $out_dir/SE/Sashimi_plot;cp $out_dir/SE/Sashimi_plot/* $report_out/Custom_Figures/SE/"\n};
		
		my $finish   = qq{touch $result/run/$case\_vs_$control.finish\n};

		open LOG, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $run/$case\_$control.sh\n";
		print LOG $cmd;
		print LOG $add;
		print LOG $fmt;
		print LOG $top10;
		print LOG $A3SSplot;
		print LOG $A5SSplot;
		print LOG $MXEplot;
		print LOG $RIplot;
		print LOG $SEplot;
		print LOG $finish;
		close LOG;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log\n};
		$pm->finish;
	}
	$pm->wait_all_children;
	print qq{mRNA 可变剪切分析已经完成!\n};
}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;

	my $organ        = qq{$metadata->{organ}};
	my $ref_mRNA_gtf = qq{$base->{$organ}{genome_mRNA_gtf}};
	# reference file is exists
	die "reference mRNA gtf file is not exists!\n" if not -e $ref_mRNA_gtf;  
	# group config is exists 
	my @groups = @{$metadata->{'group'}};
	die "group config info is not exists!\n" if scalar @groups == 0;

}

sub res_check
{
	my $out    = shift;
	my $groups = shift;
	my @temp   = ();
	foreach my $x (@{$groups}) {
		my $control = $x->[0];
		my $case    = $x->[1];
		next if  -e qq{$out/$case\_vs_$control.finish};
		push @temp, $x;
	}
	return @temp;
}

1;
