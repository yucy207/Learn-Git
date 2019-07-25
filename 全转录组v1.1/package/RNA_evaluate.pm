package package::RNA_evaluate;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use Encode qw/decode/;
use package::format;

sub run
{
	my $metadata   = shift;
	my $base       = shift;
	my @samples    = split /,/, $metadata->{'samples'};

	my $organ      = qq{$metadata->{organ}};
	my $result     = qq{$metadata->{project}/RNA_evaluate};
	my $report     = qq{$metadata->{report}/02_Mapping_Evaluation};
	my $map        = qq{$metadata->{project}/mapping/result};

	my $rscript    = qq{$base->{rscript_bin}};
	my $util       = qq{$base->{'util'}};
	my $rnaseqc    = qq{$base->{rnaseqc_bin}};

	my $ref_fasta   = qq{$base->{$organ}{genome_fasta}};
	my $ref_gencode = qq{$base->{$organ}{genome_gencode}}; 

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
		print qq{RNA_seQC 已经运行完成!\n};
		return;
	}

	pre_check($metadata, $base);
	
	my $max_threads = $base->{'thread_map'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_samples) {

		my $pid = $pm->start and next;
		system qq{mkdir -p $result/result/$x}                      if not -d qq{$result/result/$x};
		system qq{mkdir -p $report/Mapping_Evaluation_Figures/$x}  if not -d qq{$report/Mapping_Evaluation_Figures/$x};

		my $RNA_seQc_cmd = qq{$rnaseqc java -Xmx4g -jar /softwares/RNA-SeQC_v1.1.8.jar -noDoC -t $ref_gencode -s '$x|$map/$x/accepted_hits.bam|$x' -r $ref_fasta -o $result/result/$x};
		
		open  SAVE, qq{>$result/run/$x.sh} or die "Can't open $result/run/$x.sh\n";
		print SAVE qq{$RNA_seQc_cmd\n};
		close SAVE; 
		
		#print qq{$RNA_seQc_cmd\n};
		system qq{$RNA_seQc_cmd &>> $result/log/$x.log};
		system qq{rm -r $result/result/$x/$x $result/result/$x/refGene.txt*} if -e qq{$result/result/$x/metrics.tsv};

		$pm->finish;
	} 

	$pm->wait_all_children;
	
	open SAVE, qq{>$result/result/RNA_seQc.xls};
    print SAVE qq{Sample\tIntragenic Rate\tExonic Rate\tIntronic Rate\tIntergenic Rate\trRNA\trRNA rate\n};
    foreach my $x (@samples) {

		my $res     = qq{$result/result/$x/metrics.tsv};
		my @content = parse_report($res);
		my $result  = join "\t", @content;
		print SAVE qq{$result\n};

    }
    close SAVE;

	foreach my $x (@samples) {

		my $cmd = qq{$rscript Rscript $util/map_region.R $result/result/RNA_seQc.xls $report/Mapping_Evaluation_Figures/$x $x &>> $result/log/$x.log};
		system qq{$cmd\n};

	}

	my $barplot_cmd = qq{$rscript Rscript $util/map_barplot.r $result/result/RNA_seQc.xls $report/Mapping_Evaluation_Figures  &>> $result/log/barplot.log};
	system qq{$barplot_cmd\n};

	my $excel      = qq{$result/result/RNAseQc_Summary.xlsx};
	my $workbook   = Excel::Writer::XLSX->new($excel);
	my %format     = package::format::run($workbook);
	my $worksheet  = $workbook->add_worksheet("Summary");

	my $title =<<EOF;
Sample
Note
End 2 Mapping Rate
Chimeric Pairs
Intragenic Rate
Exonic Rate
Mapping Rate
Genes Detected
Unique Rate of Mapped
Read Length
End 1 Mismatch Rate
Fragment Length StdDev
Estimated Library Size
Mapped
Intergenic Rate
Total Purity Filtered Reads Sequenced
rRNA
Failed Vendor QC Check
Transcripts Detected
Mapped Pairs
Unpaired Reads
Intronic Rate
Mapped Unique Rate of Total
Expression Profiling Efficiency
Mapped Unique
End 2 Mismatch Rate
End 2 Antisense
Alternative Aligments
End 2 Sense
Fragment Length Mean
End 1 Antisense
Split Reads
Base Mismatch Rate
End 1 Sense
End 1 % Sense
rRNA rate
End 1 Mapping Rate
Duplication Rate of Mapped
End 2 % Sense
EOF

	my @rows = split /\n/, $title;
	my $row = 0;
	foreach my $x (0..38) {
		my @row = ();
		foreach my $y (@samples) {

			my $res     = qq{$result/result/$y/metrics.tsv};
			my @content = parse_metrics($res);
			push @row, $content[$x];

		}
		unshift @row, $rows[$x];
		$row == 0 ? $worksheet->write_row( $row, 0, \@row, $format{'title'}) : $worksheet->write_row( $row, 0, \@row, $format{'normal'});
		$row++;
		
	}
	$workbook->close();

	my $excel      = qq{$report/Mapping_Evaluation_Summary.xlsx};
	my $workbook   = Excel::Writer::XLSX->new($excel);
	my %format     = package::format::run($workbook);
	my $worksheet  = $workbook->add_worksheet("Mapping Table");

	my $title =  qq{Sample\t# of total reads\t# of mapped Reads\tMapped reads(%)\tIntragenic Rate(%)\tExonic Rate(%)\tIntronic Rate(%)\tIntergenic Rate(%)\trRNA\trRNA rate(%)};
	my @head  = split /\t/, $title;
	$worksheet->write_row( 0, 0, \@head, $format{'title'});
	my $row = 1;
	foreach my $x (@samples) {

		my ($mapped_reads, $total_reads, $mapped_fraction) = cal_mapping(qq{$map/$x/align_summary.txt});
		$mapped_fraction = sprintf "%.2f%%", $mapped_fraction * 100;
		my $line     = qq{$x\t$total_reads\t$mapped_reads\t$mapped_fraction\t};
		my $res      = qq{$result/result/$x/metrics.tsv};
		my @content  = parse_report($res);
		shift @content;
		$content[0]  = sprintf "%.2f%%", $content[0] * 100;
		$content[1]  = sprintf "%.2f%%", $content[1] * 100;
		$content[2]  = sprintf "%.2f%%", $content[2] * 100;
		$content[3]  = sprintf "%.2f%%", $content[3] * 100;
		#$content[4] = sprintf "%.2f%%", $content[4] * 100;
		$content[5]  = sprintf "%.2f%%", $content[5] * 100;
		my $result   = join "\t", @content;
		$line .= $result;
		my @result = split /\t/, $line;
		$worksheet->write_row( $row, 0, \@result, $format{'normal'});
		$row++;
	}
	my $worksheet1 = $workbook->add_worksheet("README");
	my $A1 = decode("utf8", "标题"); my $A2 = decode("utf8", "说明");
	my $B1 = decode("utf8", "Sample"); my $B2 = decode("utf8", "样本");
	my $C1 = decode("utf8", "# of total reads"); my $C2 = decode("utf8", "QC后总reads数量");
	my $D1 = decode("utf8", "# of mapped reads"); my $D2 = decode("utf8", "比对上的reads数量");
	my $E1 = decode("utf8", "Mapped reads(%)"); my $E2 = decode("utf8", "比对reads所占百分比");
    my $F1 = decode("utf8", "Intragenic Rate(%)"); my $F2 = decode("utf8", "比对到基因区域的reads所占百分比");
	my $G1 = decode("utf8", "Exonic Rate(%)"); my $G2 = decode("utf8", "比对到基因外显子区域的reads所占百分比");
	my $H1 = decode("utf8", "Intronic Rate(%)"); my $H2 = decode("utf8", "比对到基因内含子区域的reads所占百分比");
	my $I1 = decode("utf8", "Intergenic Rate(%)"); my $I2 = decode("utf8", "比对到基因间区域的reads所占百分比");
	my $J1 = decode("utf8", "rRNA"); my $J2 = decode("utf8", "比对到核糖体的reads数量");
	my $K1 = decode("utf8", "rRNA Rate(%)"); my $K2 = decode("utf8", "比对到核糖体的reads所占百分比");
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
	$worksheet1->write_row( 10, 0, [$K1, $K2], $format{'normal'});
	$workbook->close();

	print qq{RNA QC 已经运行完成!\n};

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
		next if -e qq{$output/$x/metrics.tsv};
		push @temp_sample, $x;
	}
	return @temp_sample;

}

sub parse_report {
	my $res = shift;
	open EXO, $res;
	while (<EXO>) {
		chomp;
		next if /^Sample/;
		my @arr  = split /\t/;
		return @arr[0,4,5,21,14,16,35];
	}
	close EXO;
}

sub cal_mapping
{

	my $align  = shift;
	open EXO, $align or die "Can't open $align\n";
	local $/ = undef;
	my $content = <EXO>;
	close EXO;
	my ($pair_reads)     = $content =~ /(\d+)\s+reads; of these/; 
	my ($unmapped_reads) = $content =~ /(\d+)\s+\(.+?\)\s+aligned\s+0\s+times/;
	my $total_reads      = $pair_reads * 2;
	my $mapped_reads     = $total_reads - $unmapped_reads;
	my $fraction         = $mapped_reads / $total_reads;
	return ($mapped_reads, $total_reads, $fraction);
	
}

sub parse_metrics
{
	my $res = shift;
	open EXO, $res;
	while (<EXO>) {
		chomp;
		next if /^Sample/;
		my @arr = split /\t/;
		return @arr;
	}
	close EXO;
}

1;

