package package::association::mRNA_miRNA;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $organ    = $metadata->{'organ'};

	my $db       = qq{$base->{$organ}{mRNA_miRNA_db}};
	my $util     = $base->{'util'};
	my $rscript  = $base->{'rscript_bin'};

	my $miRNA_result = qq{$metadata->{miRNA_project}/Diff_Expression/result};
	my $mRNA_result  = qq{$metadata->{project}/mRNA/deseq/result};
	my $mRNA_count   = qq{$metadata->{project}/mRNA/gene_profile/result/gene_count_matrix.csv};
	my $report       = qq{$metadata->{report}/06_Association_Analysis/mRNA_miRNA};
	my $result       = qq{$metadata->{project}/association/mRNA_miRNA};
	
	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;

	my @res_groups = res_check(qq{$report}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{mRNA 和 miRNA 相互作用已经分析完成!\n};
		#return;
	}

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $pid = $pm->start and next;

		my $control  = $x->[0];
		my $case     = $x->[1];
		my $base_out = qq{$case\_vs_$control};
		
		system qq{mkdir -p $result/$base_out} if not -d qq{$result/$base_out};
		system qq{mkdir -p $report/$base_out} if not -d qq{$report/$base_out};

		my ($control_samples, $case_samples) = split /;/, $x->[2];
		my @control_sample = split /,/, $control_samples;
		my @case_sample    = split /,/, $case_samples;
		my $all_samples    = qq{$control_samples,$case_samples};

		my $mRNA_xls    = qq{$mRNA_result/$base_out/gene/diff.xls};
		my $miRNA_xls   = qq{$miRNA_result/$base_out/Known/sig_deseq_known.xls};
		my $miRNA_count = qq{$miRNA_result/$base_out/Known/known.count.xls};

		# extract count profile
		system qq{$rscript Rscript $util/extract_diff_count.R $control_samples $case_samples $mRNA_count $mRNA_xls $result/$base_out/mRNA.count.xls};
		system qq{perl $util/extract_profile_by_name.pl $miRNA_count $all_samples >$result/$base_out/miRNA.count.xls\n};
		#calculate co_expression
		system qq{perl $util/submit_co_expression.pl $result/$base_out/mRNA.count.xls $result/$base_out/miRNA.count.xls $result/$base_out/ mRNA miRNA $util/co_expression.R $rscript Rscript};

		my %diff_miRNA = ();
		my %diff_mRNA  = ();

		open EXO, $miRNA_xls or die "Can't open $miRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_miRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		open EXO, $mRNA_xls or die "Can't open $mRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_mRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		my %mRNA_with_miRNA = ();
		open EXO, qq{$result/$base_out/correlation.xls} or die "Can't open $result/$base_out/correlation.xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			$mRNA_with_miRNA{qq{$arr[0]:$arr[1]}} = [$arr[2], $arr[3]];

		}
		close EXO;

		##################################### report ###################################
		my $excel    = qq{$report/$base_out/MRNA_with_MiRNA_Interaction.xlsx};
		my $workbook = Excel::Writer::XLSX->new($excel);
		my %format   = package::format::run($workbook);

		# add worksheet "interaction"
		my @interaction_keys = ();
		my %interaction      = ();

		my $worksheet = $workbook->add_worksheet(qq{Interaction});
		my $row = 0;
		my $cnt = 0;
		open DB, $db or die "Can't open $db!\n";
		while (<DB>) {
			chomp;
			my @arr = split /\t/;
			if ($cnt == 0 ) {
				$worksheet->write_row( $row, 0, \@arr, $format{'title'});
			} else {
				next if not exists $diff_miRNA{$arr[2]};
				next if not exists $diff_mRNA{$arr[0]};
				$interaction{qq{$arr[0]:$arr[2]}} = $_;
				push @interaction_keys, qq{$arr[0]:$arr[2]};
				$worksheet->write_row( $row, 0, \@arr, $format{'normal'});
			}
			$cnt++;
			$row++;
		}
		close DB;

		# add worksheet "co_expression"
		my %co_expression = ();	
		my $worksheet1 = $workbook->add_worksheet(qq{Coexpressionion});
		my @head = qw/mRNA miRNA cor pvalue/;
		$worksheet1->write_row(0, 0, \@head, $format{'title'});

		my $row = 1;
		foreach my $x (@interaction_keys) {
			my ($mRNA, $miRNA) = split /:/, $x;
			my $cor  = $mRNA_with_miRNA{$x}->[0];
			my $pval = $mRNA_with_miRNA{$x}->[1];
			my @arr = ($mRNA, $miRNA, $cor, $pval);
			$co_expression{$x} = [$cor, $pval];
			$worksheet1->write_row($row, 0, \@arr, $format{'normal'});
			$row++;
		}

		# add worksheet "Regulation"
		my %unique = ();
		
		open SAVE, qq{>$result/$base_out/network.xls} or die "Can't open $result/$base_out/network.xls!\n";
		open EDEG, qq{>$report/$base_out/edge.txt} or die "Can't open $report/$base_out/edge.txt!\n";
		open NODE, qq{>$report/$base_out/node.txt} or die "Can't open $report/$base_out/node.txt!\n";
		
		my $worksheet2 = $workbook->add_worksheet(qq{Regulation});
		my @head = qw/mRNA miRNA mRNA_log2(foldchange) mRNA_pvalue mRNA_regulation miRNA_log2(foldchange) miRNA_pvalue miRNA_regulation cor cor_pvalue/;
		$worksheet2->write_row(0, 0, \@head, $format{'title'});
		
		my $row = 1;
		foreach my $x (keys %interaction) {
			next if not exists $co_expression{$x};
			my ($mRNA, $miRNA) = split /:/, $x;
			my @res = ();
			push @res, $mRNA;
			push @res, $miRNA;
			push @res, $diff_mRNA{$mRNA}->[0];
			push @res, $diff_mRNA{$mRNA}->[1];
			push @res, $diff_mRNA{$mRNA}->[2];

			push @res, $diff_miRNA{$miRNA}->[0];
			push @res, $diff_miRNA{$miRNA}->[1];
			push @res, $diff_miRNA{$miRNA}->[2];

			push @res, $co_expression{qq{$mRNA:$miRNA}}->[0];
			push @res, $co_expression{qq{$mRNA:$miRNA}}->[1];

			next if $co_expression{qq{$mRNA:$miRNA}}->[0] > 0;
			next if $co_expression{qq{$mRNA:$miRNA}}->[1] > 0.05;

			print SAVE qq{$interaction{qq{$mRNA:$miRNA}}\t$co_expression{qq{$mRNA:$miRNA}}->[0]\t$co_expression{qq{$mRNA:$miRNA}}->[1]\n};
			$worksheet2->write_row($row, 0, \@res, $format{'normal'});
			$row++;

			print EDEG qq{$mRNA\t$miRNA\n};
			print NODE qq{$mRNA\tmRNA\t$diff_mRNA{$mRNA}->[2]\n}    if not exists $unique{$mRNA};
			print NODE qq{$miRNA\tmiRNA\t$diff_miRNA{$miRNA}->[2]\n} if not exists $unique{$miRNA};

			$unique{$miRNA} = "-";
			$unique{$mRNA} = "-";
		}

		close EDEG;
		close NODE;
		$workbook->close();
		close SAVE;

		# cytoscape html
		system qq{perl $util/cytoscape_js.pl -node $report/$base_out/node.txt -edge $report/$base_out/edge.txt -o $report/$base_out};

		
		$pm->finish;

	}

	$pm->wait_all_children;
	print qq{mRNA 和 miRNA 相互作用已经分析完成!\n};

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
		next if -e qq{$out/$base_out/MRNA_with_MiRNA_Interaction.xlsx};	
		push @temp, $x;
	}
	return @temp;
}


1;
