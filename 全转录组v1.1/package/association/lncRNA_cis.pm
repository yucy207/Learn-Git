package package::association::lncRNA_cis;
use Excel::Writer::XLSX;
use Parallel::ForkManager;
use package::format;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};

	my $util     = qq{$base->{util}};
	my $rscript  = qq{$base->{rscript_bin}};

	my $mRNA_result   = qq{$metadata->{project}/mRNA/deseq/result};
	my $mRNA_count    = qq{$metadata->{project}/mRNA/gene_profile/result};
	my $lncRNA_result = qq{$metadata->{project}/lncRNA/deseq/result};
	my $lncRNA_count  = qq{$metadata->{project}/lncRNA/gene_profile/result};
	my $lncRNA_target = qq{$metadata->{project}/lncRNA/target_predict/result};

	my $report        = qq{$metadata->{report}/06_Association_Analysis/lncRNA_mRNA/Cis_Analysis};
	my $result        = qq{$metadata->{project}/association/lncRNA_mRNA/Cis_Analysis};
	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;

	my @res_groups = res_check( qq{$report}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{lncRNA Cis 作用已经分析完成!\n};
		return;
	}

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $pid = $pm->start and next;

		my $control  = $x->[0];
		my $case     = $x->[1];
		my $base_out = qq{$case\_vs_$control};

		system qq{mkdir -p $report/$base_out} if not -d qq{$report/$base_out};
		system qq{mkdir -p $result/$base_out} if not -d qq{$result/$base_out};

		my ($control_samples, $case_samples) = split /;/, $x->[2];
		my @control_sample = split /,/, $control_samples;
		my @case_sample    = split /,/, $case_samples;
		my $all_samples    = qq{$control_samples,$case_samples};

		my $mRNA_xls    = qq{$mRNA_count/gene_count_matrix.csv};
		my $mRNA_diff   = qq{$mRNA_result/$base_out/gene/diff.xls};
		my $lncRNA_xls  = qq{$lncRNA_count/gene_count_matrix.csv};	
		my $lncRNA_diff = qq{$lncRNA_result/$base_out/gene/diff.xls};
		my $db          = qq{$lncRNA_target/$base_out/lncRNA.target.fmt.xls};
		
		# extract diff count
		system qq{$rscript Rscript $util/extract_diff_count.R $control_samples $case_samples $mRNA_xls $mRNA_diff $result/$base_out/mRNA.count.xls};
		system qq{$rscript Rscript $util/extract_diff_count.R $control_samples $case_samples $lncRNA_xls $lncRNA_diff $result/$base_out/lncRNA.count.xls};
		#calculate co_expression
		system qq{perl $util/submit_co_expression.pl $result/$base_out/lncRNA.count.xls $result/$base_out/mRNA.count.xls $result/$base_out/ lncRNA mRNA  $util/co_expression.R $rscript Rscript};

		my %diff_lncRNA  = ();
		my %diff_mRNA = ();

		open EXO, $lncRNA_diff or die "Can't open $lncRNA_diff!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_lncRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		open EXO, $mRNA_diff or die "Can't open $mRNA_diff!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_mRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		##################################### report ###################################
		my $excel    = qq{$report/$base_out/LncRNA_with_MRNA_Interaction.xlsx};
		my $workbook = Excel::Writer::XLSX->new($excel);
		my %format   = package::format::run($workbook);

		# add worksheet "interaction"
		my $worksheet   = $workbook->add_worksheet(qq{Interaction});
		my %interaction = ();
		my $row =  0;
		my $cnt =  0;
		open DB, $db or die "Can't open $db!\n";
		while (<DB>) {
			chomp;
			my @arr = split /\t/;
			if ($cnt == 0 ) {
				$worksheet->write_row( $row, 0, \@arr, $format{'title'});
			} else {
				my ($lncRNA_id) = $arr[0] =~ /\((.+?)\)/;
				my ($mRNA_id)   = $arr[2] =~ /\((.+?)\)/;
				next if not exists $diff_lncRNA{$lncRNA_id};
				next if not exists $diff_mRNA{$mRNA_id};
				$interaction{qq{$lncRNA_id:$mRNA_id}} = "-";
				$worksheet->write_row( $row, 0, \@arr, $format{'normal'});

			}
			$cnt++;
			$row++;
		}
		close DB;

		# add worksheet "co_expression"
		my %co_expression = ();
		my $worksheet1    = $workbook->add_worksheet(qq{Coexpressionion});
		my $row = 0;
		open CORRELATION, qq{$result/$base_out/correlation.xls} or die "Can't open $result/$base_out/correlation.xls!\n";
		while (<CORRELATION>) {
			chomp;
			my @arr = split /\t/;
			if ($row == 0) {
				$worksheet1->write_row($row, 0, \@arr, $format{'title'});
			} else {
				#next if $arr[2] > 0;
				next if $arr[3] > 0.05;
				$co_expression{qq{$arr[0]:$arr[1]}} = [$arr[2], $arr[3]];
				$worksheet1->write_row($row, 0, \@arr, $format{'normal'});
			}
			$row++;

		}
		close CORRELATION;

		# add worksheet "Regulation"
		open EDGE, qq{>$report/$base_out/edge.txt} or die "Can't open $report/$base_out/edge.txt!\n";
		open NODE, qq{>$report/$base_out/node.txt} or die "Can't open $report/$base_out/node.txt!\n";
		my %unique = ();
		my $worksheet2  = $workbook->add_worksheet(qq{Regulation});
		my @head = qw/lncRNA mRNA lncRNA_log2(foldchange) lncRNA_pvalue lncRNA_regulation mRNA_log2(foldchange) mRNA_pvalue mRNA_regulation cor cor_pvalue/;
		$worksheet2->write_row(0, 0, \@head, $format{'title'});
		my $row = 1;
		foreach my $x (keys %interaction) {
			next if not exists $co_expression{$x};
			my ($lncRNA, $mRNA) = split /:/, $x;
			my @res = ();
			push @res, $lncRNA;
			push @res, $mRNA;
			push @res, $diff_lncRNA{$lncRNA}->[0];
			push @res, $diff_lncRNA{$lncRNA}->[1];
			push @res, $diff_lncRNA{$lncRNA}->[2];

			push @res, $diff_mRNA{$mRNA}->[0];
			push @res, $diff_mRNA{$mRNA}->[1];
			push @res, $diff_mRNA{$mRNA}->[2];

			print EDGE qq{$lncRNA\t$mRNA\n};
			print NODE qq{$lncRNA\tlncRNA\t$diff_lncRNA{$lncRNA}->[2]\n} if not exists $unique{$lncRNA};
			print NODE qq{$mRNA\tmRNA\t$diff_mRNA{$mRNA}->[2]\n}       if not exists $unique{$mRNA};

			$unique{$lncRNA} = "-";
			$unique{$mRNA}   = "-";

			push @res, $co_expression{qq{$lncRNA:$mRNA}}->[0];
			push @res, $co_expression{qq{$lncRNA:$mRNA}}->[1];

			$worksheet2->write_row($row, 0, \@res, $format{'normal'});
			$row++;
		}
		close EDGE;
		close NODE;
		
		$workbook->close();
		

		# cytoscape html
		system qq{perl $util/cytoscape_js.pl -node $report/$base_out/node.txt -edge $report/$base_out/edge.txt -o $report/$base_out};
		
		$pm->finish;
	}

	$pm->wait_all_children;
	print qq{lncRNA Cis 作用已经分析完成!\n};
}

sub res_check
{
	my $out    = shift;
	my $groups = shift;
	my @temp   = ();
	foreach my $x (@{$groups}) {
		my $control  = $x->[0];
		my $case     = $x->[1];
		my $base_out = qq{$case\_vs_$control};
		next if -e qq{$out/$base_out/LncRNA_with_MRNA_Interaction.xlsx};		
		push @temp, $x;
	}
	return @temp;
}


1;
