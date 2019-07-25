package package::association::lncRNA_miRNA;
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

	my $rscript  = $base->{'rscript_bin'};
	my $util     = qq{$base->{util}};
	my $db       = $base->{$organ}{'lncRNA_miRNA_db'};

	my $miRNA_result  = qq{$metadata->{miRNA_project}/Diff_Expression/result};
	my $lncRNA_result = qq{$metadata->{project}/lncRNA/deseq/result};
	my $lncRNA_count  = qq{$metadata->{project}/lncRNA/gene_profile/result/gene_count_matrix.csv};
	my $report        = qq{$metadata->{report}/06_Association_Analysis/lncRNA_miRNA};
	my $result        = qq{$metadata->{project}/association/lncRNA_miRNA};

	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;

	my @res_groups = res_check( qq{$report}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{lncRNA 和 miRNA 相互作用已经分析完成!\n};
		return;
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
		my @control_sample  = split /,/, $control_samples;
		my @case_sample     = split /,/, $case_samples;
		my $all_samples     = qq{$control_samples,$case_samples};

		my $lncRNA_xls  = qq{$lncRNA_result/$base_out/gene/diff.xls};
		my $miRNA_xls   = qq{$miRNA_result/$base_out/Known/sig_deseq_known.xls};
		my $miRNA_count = qq{$miRNA_result/$base_out/Known/known.count.xls};

		# extract count profile
		system qq{$rscript Rscript $util/extract_diff_count.R $control_samples $case_samples $lncRNA_count $lncRNA_xls $result/$base_out/lncRNA.count.xls};
		system qq{perl $util/extract_profile_by_name.pl $miRNA_count $all_samples >$result/$base_out/miRNA.count.xls\n};
		#　calculate co_expression
		system qq{perl $util/submit_co_expression.pl $result/$base_out/lncRNA.count.xls $result/$base_out/miRNA.count.xls $result/$base_out/ lncRNA miRNA  $util/co_expression.R $rscript Rscript} if not -e qq{$result/$base_out/correlation.xls};
		
		my %diff_miRNA  = ();
		my %diff_lncRNA = ();

		open EXO, $miRNA_xls or die "Can't open $miRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_miRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		open EXO, $lncRNA_xls or die "Can't open $lncRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_lncRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		my %lncRNA_with_miRNA = ();
		open EXO, qq{$result/$base_out/correlation.xls} or die "Can't open $result/$base_out/correlation.xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			$lncRNA_with_miRNA{qq{$arr[0]:$arr[1]}} = [$arr[2], $arr[3]];

		}
		close EXO;
	
		##################################### report ###################################
		my $excel      = qq{$report/$base_out/LncRNA_with_MiRNA_Interaction.xlsx};
		my $workbook   = Excel::Writer::XLSX->new($excel);
		my %format     = package::format::run($workbook);

		my $worksheet  = $workbook->add_worksheet(qq{Interaction});
		my $worksheet1 = $workbook->add_worksheet(qq{Coexpression});
		my $worksheet2 = $workbook->add_worksheet(qq{Regulation});

		my $row     =  0;
		my $cnt     = 0;
		my $res_cnt = 1;
		my %unique  = ();

		open DB, $db or die "Can't open $db!\n";	
		open SAVE, qq{>$result/$base_out/network.xls} or die "Can't open $result/$base_out/network.xls\n";
		open EDEG, qq{>$report/$base_out/edge.txt} or die "Can't open $report/$base_out/edge.txt!\n";
		open NODE, qq{>$report/$base_out/node.txt} or die "Can't open $report/$base_out/node.txt!\n";

		while (<DB>) {
			chomp;
			my @arr = split /\t/;
			if ($cnt == 0 ) {
				$worksheet->write_row( $row, 0, \@arr, $format{'title'});
				my @head1 = qw/lncRNA miRNA cor pvalue/;
				my @head2 = qw/lncRNA miRNA lncRNA_log2(foldchange) lncRNA_pvalue lncRNA_regulation miRNA_log2(foldchange) miRNA_pvalue miRNA_regulation cor pvalue/;
				$worksheet1->write_row(  $row, 0, \@head1, $format{'title'});
				$worksheet2->write_row( $row, 0, \@head2, $format{'title'});

			} else {
				next if not exists $diff_miRNA{$arr[1]};
				my ($lncRNA_id) = $arr[0] =~ /\((.+?)\)/;
				next if not exists $diff_lncRNA{$lncRNA_id};

				$worksheet->write_row( $row, 0, \@arr, $format{'normal'});

				my @worksheet1_line = ($lncRNA_id, $arr[1], $lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[0], $lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[1]);
				$worksheet1->write_row( $row, 0, \@worksheet1_line, $format{'normal'});

				my @res = ();
				push @res, $lncRNA_id;
				push @res, $arr[1];
				push @res, $diff_lncRNA{$lncRNA_id}->[0];
				push @res, $diff_lncRNA{$lncRNA_id}->[1];
				push @res, $diff_lncRNA{$lncRNA_id}->[2];

				push @res, $diff_miRNA{$arr[1]}->[0];
				push @res, $diff_miRNA{$arr[1]}->[1];
				push @res, $diff_miRNA{$arr[1]}->[2];
				push @res, $lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[0];
				push @res, $lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[1];

				next if $lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[0] > 0;
				next if $lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[1] > 0.5;

				print SAVE $result_line = qq{$_\t$lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[0]\t$lncRNA_with_miRNA{qq{$lncRNA_id:$arr[1]}}->[1]\n};
				$worksheet2->write_row( $res_cnt, 0, \@res, $format{'normal'});
				$res_cnt++;

				print EDEG qq{$lncRNA_id\t$arr[1]\n};
				print NODE qq{$lncRNA_id\tlncRNA\t$diff_lncRNA{$lncRNA_id}->[2]\n} if not exists $unique{$lncRNA_id};
				print NODE qq{$arr[1]\tmiRNA\t$diff_miRNA{$arr[1]}->[2]\n}         if not exists $unique{$arr[1]};

				$unique{$lncRNA_id} = "-";
				$unique{$arr[1]}    = "-";
			}

			$cnt++;
			$row++;
		}

		close DB;
		close EDEG;
		close NODE;
		$workbook->close();
		close SAVE;


		# cytoscape html
		system qq{perl $util/cytoscape_js.pl -node $report/$base_out/node.txt -edge $report/$base_out/edge.txt -o $report/$base_out/};

		
		$pm->finish;



	}

	$pm->wait_all_children;
	print qq{lncRNA 和 miRNA 相互作用已经分析完成!\n};

}

sub res_check
{
	my $out    = shift;
	my $groups = shift;
	my @temp   = ();
	foreach my $x (@{$groups}) {
		my $control  = $x->[0];
		my $case     = $x->[1];
		my $base_out = scalar @{$groups} == 1 ? '' : qq{$case\_vs_$control};
		next if -e qq{$out/$base_out/LncRNA_with_MiRNA_Interaction.xlsx};	
		push @temp, $x;
	}
	return @temp;
}

1;
