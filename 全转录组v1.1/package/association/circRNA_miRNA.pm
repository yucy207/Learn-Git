package package::association::circRNA_miRNA;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $rscript  = $base->{'rscript_bin'};
	my $util     = qq{$base->{util}};

	my $circRNA_result = qq{$metadata->{project}/circRNA/different_analysis/result};
	my $miRNA_result   = qq{$metadata->{miRNA_project}/Diff_Expression/result};
	my $circRNA_target = qq{$metadata->{project}/circRNA/circrna_target/result};
	my $report         = qq{$metadata->{report}/06_Association_Analysis/circRNA_miRNA};
	my $result         = qq{$metadata->{project}/association/circRNA_miRNA};

	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;

	my @res_groups = res_check(qq{$report}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{circRNA 和 miRNA 相互作用已经分析完成!\n};
		return;
	}

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $control  = $x->[0];
		my $case     = $x->[1];
		my $base_out = qq{$case\_vs_$control};

		system qq{mkdir -p $result/$base_out} if not -d qq{$result/$base_out};
		system qq{mkdir -p $report/$base_out} if not -d qq{$report/$base_out};

		my ($control_samples, $case_samples) = split /;/, $x->[2];
		my @control_sample = split /,/, $control_samples;
		my @case_sample    = split /,/, $case_samples;
		my $all_samples    = qq{$control_samples,$case_samples};

		my $miRNA_xls      = qq{$miRNA_result/$base_out/Known/sig_deseq_known.xls};
		my $miRNA_count    = qq{$miRNA_result/$base_out/Known/known.count.xls};
		my $circRNA_count  = qq{$circRNA_result/$base_out/gene/circRNA.different.expression.xls};	
		my $db             = qq{$circRNA_target/$base_out/miranda.result.xls};

		# extract count profile
		system qq{perl $util/extract_profile_by_name.pl $circRNA_count $all_samples > $result/$base_out/circRNA.count.xls\n};
		system qq{perl $util/extract_profile_by_name.pl $miRNA_count   $all_samples > $result/$base_out/miRNA.count.xls\n};
		#　calculate co_expression
		system qq{perl $util/submit_co_expression.pl $result/$base_out/circRNA.count.xls $result/$base_out/miRNA.count.xls $result/$base_out/ circRNA miRNA  $util/co_expression.R $rscript Rscript};
		
		my %diff_miRNA = ();
		my %diff_circRNA = ();

		open EXO, $miRNA_xls or die "Can't open $miRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$arr[0] =~ s/\s+//g;	
			$diff_miRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		open EXO, $circRNA_count or die "Can't open $circRNA_count!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_circRNA{$arr[0]} = [$arr[$#arr - 3], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		my %circRNA_with_miRNA = ();
		open EXO, qq{$result/$base_out/correlation.xls} or die "Can't open $result/$base_out/correlation.xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			my $key = qq{$arr[0]:$arr[1]};
			$circRNA_with_miRNA{$key} = [$arr[2], $arr[3]];
		}
		close EXO;

		##################################### report ###################################
		my $excel      = qq{$report/$base_out/CircRNA_with_MiRNA_Interaction.xlsx};
		my $workbook   = Excel::Writer::XLSX->new($excel);
		my %format     = package::format::run($workbook);

		my $worksheet  = $workbook->add_worksheet(qq{Interaction});
		my $worksheet1 = $workbook->add_worksheet(qq{Coexpression});
		my $worksheet2 = $workbook->add_worksheet(qq{Regulation});

		my $row = 0;
		my $cnt = 0;
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
				my @head1 = qw/circRNA miRNA cor pvalue/;
				my @head2 = qw/circRNA miRNA circRNA_log2(foldchange) circRNA_pvalue circRNA_regulation miRNA_log2(foldchange) miRNA_pvalue miRNA_regulation cor pvalue/;
				$worksheet1->write_row(  $row, 0, \@head1, $format{'title'});
				$worksheet2->write_row( $row, 0, \@head2, $format{'title'});

			} else {

				$arr[0] =~ s/\|/-/;
				next if not exists $diff_miRNA{$arr[1]};
				next if not exists $diff_circRNA{$arr[0]};
				#print qq{$arr[0]\t$arr[1]\t$diff_miRNA{$arr[1]}->[0]\t$diff_circRNA{$arr[0]}->[0]\n};
				$worksheet->write_row( $row, 0, \@arr, $format{'normal'});
				
				my @worksheet1_line = ($arr[0], $arr[1], $circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[0], $circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[1]);
				$worksheet1->write_row( $row, 0, \@worksheet1_line, $format{'normal'});

				my @res = ();
				push @res, $arr[0];
				push @res, $arr[1];
				push @res, $diff_circRNA{$arr[0]}->[0];
				push @res, $diff_circRNA{$arr[0]}->[1];
				push @res, $diff_circRNA{$arr[0]}->[2];

				push @res, $diff_miRNA{$arr[1]}->[0];
				push @res, $diff_miRNA{$arr[1]}->[1];
				push @res, $diff_miRNA{$arr[1]}->[2];
				push @res, $circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[0];
				push @res, $circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[1];

				print EDEG qq{$arr[0]\t$arr[1]\n};
				print NODE qq{$arr[0]\tcircRNA\t$diff_circRNA{$arr[0]}->[2]\n} if not exists $unique{$arr[0]};
				print NODE qq{$arr[1]\tmiRNA\t$diff_miRNA{$arr[1]}->[2]\n}     if not exists $unique{$arr[1]};

				$unique{$arr[0]} = "-";
				$unique{$arr[1]}    = "-";
				next if $circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[0] > 0;
				#next if $circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[1] > 0.5;

				print SAVE $result_line = qq{$_\t$circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[0]\t$circRNA_with_miRNA{qq{$arr[0]:$arr[1]}}->[1]\n};
				$worksheet2->write_row( $res_cnt, 0, \@res, $format{'normal'});
				$res_cnt++;

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
	}

	print qq{circRNA 和 miRNA 相互作用已经分析完成!\n};

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
		next if -e qq{$out/$base_out/LncRNA_with_MiRNA_Interaction.xlsx};		
		push @temp, $x;
	}
	return @temp;
}

1;
