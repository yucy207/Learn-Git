package package::reads_evaluate;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use Encode qw/decode/;
use package::format;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};

	my $report   = qq{$metadata->{report}/01_Quality_Statistic/After_QC_Quality_Statistic};
	my $result   = qq{$metadata->{project}/trim};
	my $raw      = qq{$metadata->{raw_data}};

	my $rscript  = qq{$base->{rscript_bin}};
	my $util     = qq{$base->{util}};
	my $fastqc   = qq{$base->{fastqc_bin}};
	my $seqtk    = qq{$base->{seqtk_bin}};

	system qq{mkdir -p $report}       if not -d $report;
	system qq{mkdir -p $result/fastx} if not -d qq{$result/fastx};
	
	my @res_samples = res_check(qq{$result/fastx}, \@samples);

	if (exists $base->{'force_sample'}) {
		foreach my $x (split /,/, $base->{'force_sample'}) {
			push @res_samples, $x;
		}
	}
	if (exists $base->{'force_step'}) {
		my @steps = split /,/, $base->{'force_step'};
		@res_samples = @samples if 3 ~~ @steps;
	}

	if (scalar @res_samples == 0) {

		print qq{质量控制已经运行完成!\n};
		return;		
	}

	my $max_threads = qq{$base->{thread_qc}};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_samples) {

		my $pid = $pm->start and next;
		system qq{mkdir -p $result/fastqc/$x} if not -d qq{$result/fastqc/$x};
		system qq{mkdir -p $report/$x}        if not -d qq{$report/$x};

		my $fastqc_cmd    = qq{$fastqc fastqc -t 50 -o $result/fastqc/$x --extract  $result/trim_galore/$x/$x\_R1_val_1.fq.gz $result/trim_galore/$x/$x\_R2_val_2.fq.gz &>> $result/log/$x.log};
		my $fastx_r1_cmd  = qq{$seqtk fqchk $result/trim_galore/$x/$x\_R1_val_1.fq.gz >  $result/fastx/$x.R1.quality.stat.txt};
		my $fastx_r2_cmd  = qq{$seqtk fqchk $result/trim_galore/$x/$x\_R2_val_2.fq.gz >  $result/fastx/$x.R2.quality.stat.txt};
		my $cat           = qq{(grep -v "ALL\\|POS\\|min_len" $result/fastx/$x.R1.quality.stat.txt;grep -v "ALL\\|POS\\|min_len" $result/fastx/$x.R2.quality.stat.txt) > $result/fastx/$x.quality.stat};

		my $base_com      = qq{$rscript Rscript  $util/base_composition.R $result/fastx/$x.quality.stat $report/$x/$x.base.composition.pdf $x &>>$result/log/$x.log};
		my $error_rate    = qq{$rscript Rscript  $util/error_rate.R $result/fastx/$x.quality.stat $report/$x/$x.error.rate.pdf $x &>>$result/log/$x.log };
		my $cp_r1         = qq{cp $result/fastqc/$x/$x\_R1_val_1_fastqc/Images/per_base_quality.png $report/$x/$x.R1.qual.png};
		my $cp_r2         = qq{cp $result/fastqc/$x/$x\_R2_val_2_fastqc/Images/per_base_quality.png $report/$x/$x.R2.qual.png};
			
		my $finish        = qq{touch $result/fastx/$x.finish};

		open SAVE, qq{>$result/run/$x.sh};
		print SAVE qq{$fastqc_cmd\n};
		print SAVE qq{$fastx_r1_cmd\n$fastx_r2_cmd\n$cat\n$cp_r1\n$cp_r2\n$base_com\n$error_rate\n};
		print SAVE qq{$finish\n};
		close SAVE;

		system qq{bash $result/run/$x.sh > $result/log/$x.log\n};

		$pm->finish;	

	}

	$pm->wait_all_children;

	my $excel     = qq{$metadata->{report}/01_Quality_Statistic/Quality_Statistic_Summary.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);
	my $worksheet = $workbook->add_worksheet("QC_Table");
		
	my $title = qq{Sample\t# of raw reads\t# of raw bases\tQ20\tQ30\traw reads GC(%)\t# of clean reads\t# of clean bases\tClean Reads(%)\tclean reads GC(%)};
	my @head  = split /\t/, $title;
	$worksheet->write_row( 0, 0, \@head, $format{'title'});

	my $row = 1;

	foreach my $x (@samples) {

		my ($raw_reads, $raw_bases, $raw_q20, $raw_q30, $raw_gc, $clean_reads, $clean_bases, $clean_q20, $clean_q30, $clean_gc) = prase_fastp(qq{$result/trim_galore/$x/$x.json});

		my $clean_frac     = sprintf "%.2f%%", $clean_reads  / $raw_reads     * 100;
		my $raw_gc_frac    = sprintf "%.2f%%", $raw_gc *  100;
		my $clean_gc_frac  = sprintf "%.2f%%", $clean_gc * 100;
		my $q20_frac       = sprintf "%.2f%%", $raw_q20 * 100;
		my $q30_frac       = sprintf "%.2f%%", $raw_q30 * 100;

		my $line = qq{$x\t$raw_reads\t$raw_bases\t$q20_frac\t$q30_frac\t$raw_gc_frac\t$clean_reads\t$clean_bases\t$clean_frac\t$clean_gc_frac\n};
		my @result = split /\t/, $line;
		$worksheet->write_row( $row, 0, \@result, $format{'normal'});
		$row++;
	}

	my $worksheet1 = $workbook->add_worksheet("README");
	my $A1 = decode("utf8", "标题"); my $A2 = decode("utf8", "说明");
	my $B1 = decode("utf8", "# of raw reads"); my $B2 = decode("utf8", "原始序列数");
	my $C1 = decode("utf8", "# of raw bases"); my $C2 = decode("utf8", "原始碱基数");
	my $D1 = decode("utf8", "Q20"); my $D2 = decode("utf8", "在Phred（33）标准下，错误率小于0.01的碱基所占百分比");
	my $E1 = decode("utf8", "Q30"); my $E2 = decode("utf8", "在Phred（33）标准下，错误率小于0.001的碱基所占百分比");
	my $F1 = decode("utf8", "raw reads GC(%)"); my $F2 = decode("utf8", "原始序列GC含量百分比");
	my $G1 = decode("utf8", "# of clean reads"); my $G2 = decode("utf8", "质控后的序列数");
	my $H1 = decode("utf8", "# of clean bases"); my $H2 = decode("utf8", "质控后的碱基数");
	my $I1 = decode("utf8", "Clean reads(%)"); my $I2 = decode("utf8", "质控序列所占百分比");
	my $J1 = decode("utf8", "Clean reads GC(%)"); my $J2 = decode("utf8", "质控后序列GC含量百分比");

	$worksheet1->write_row( 0, 0, [$A1, $A2], $format{'title'});
	$worksheet1->write_row( 1, 0, [$B1, $B2], $format{'normal'});
	$worksheet1->write_row( 2, 0, [$C1, $C2], $format{'normal'});
	$worksheet1->write_row( 3, 0, [$D1, $D2], $format{'normal'});
	$worksheet1->write_row( 4, 0, [$E1, $E2], $format{'normal'});
	$worksheet1->write_row( 5, 0, [$F1, $F2], $format{'normal'});
	$worksheet1->write_row( 6, 0, [$G1, $G2], $format{'normal'});
	$worksheet1->write_row( 7, 0, [$H1, $H2], $format{'normal'});
	$worksheet1->write_row( 8, 0, [$I1, $I2], $format{'normal'});
	$worksheet1->write_row( 9, 0, [$J1, $J2], $format{'normal'});	

	$workbook->close();		

	print qq{数据量统计已经运行完成!\n};

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

sub prase_fastp
{
	my $txt  = shift;

	my $cnt = 0;
	my @res = ();
	open TXT, $txt or die "Can't open $txt!\n";
	while (<TXT>) {
		$cnt++;

		if (/total_reads/) {
			my ($raw_reads)  = $_ =~ /"total_reads":(\d+),/;
			push @res, $raw_reads;
			next;
		}

		if (/total_bases/) {
			my ($raw_bases)  = $_ =~  /"total_bases":(\d+),/;
			push @res, $raw_bases;
			next;
		}



		if (/q20_rate/) {
			my ($raw_q20)    = $_ =~  /"q20_rate":(\S+?),/;
			push @res, $raw_q20;
			next;
		}


		if (/q30_rate/) {
			my ($raw_q30)    = $_ =~  /"q30_rate":(\S+?),/;
			push @res, $raw_q30;
			next;
		}

		if (/gc_content/) {
			my ($raw_gc)     = $_ =~  /"gc_content":(\S+)/;
			push @res, $raw_gc;
			next;
		}



		last if $cnt >= 18;

	}
	close TXT;

	return @res;
}

1;
