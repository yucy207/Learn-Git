package package::ceRNA::mRNA_miRNA_circRNA;
use package::format;
use Excel::Writer::XLSX;
use Parallel::ForkManager;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $organ    = qq{$metadata->{organ}};

	my $mRNA_with_miRNA_db = qq{$base->{$organ}{mRNA_miRNA_db}};

	my $util     = qq{$base->{util}};
	my $rscript  = qq{$base->{rscript_bin}};

	my $miRNA_result   = qq{$metadata->{miRNA_project}/Diff_Expression/result};
	my $mRNA_result    = qq{$metadata->{project}/mRNA/deseq/result};
	my $circRNA_result = qq{$metadata->{project}/circRNA/different_analysis/result};
	my $circRNA_target = qq{$metadata->{project}/circRNA/circrna_target/result};

	my $report = qq{$metadata->{report}/07_CeRNA_Analysis/mRNA_miRNA_circRNA};
	my $result = qq{$metadata->{project}/ceRNA/mRNA_miRNA_circRNA};

	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;

	my @res_groups = res_check( qq{$report}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{mRNA miRNA circRNA ceRNA已经分析完成!\n};
		return;
	}

	foreach my $x (@res_groups) {

		my $control  = $x->[0];
		my $case     = $x->[1];
		my $base_out = qq{$case\_vs_$control};

		system qq{mkdir -p $report/$base_out} if not -d qq{$report/$base_out/};
		system qq{mkdir -p $result/$base_out} if not -d qq{$result/$base_out/};

		my ($control_samples, $case_samples) = split /;/, $x->[2];
		my @control_sample = split /,/, $control_samples;
		my @case_sample    = split /,/, $case_samples;
		my $all_samples    = qq{$control_samples,$case_samples};

		my $mRNA_xls    = qq{$mRNA_result/$base_out/gene/diff.xls};
		my $miRNA_xls   = qq{$miRNA_result/$base_out/Known/sig_deseq_known.xls};
		my $circRNA_xls = qq{$circRNA_result/$base_out/gene/circRNA.different.expression.xls};
		my $mRNA_count  = qq{$metadata->{project}/mRNA/gene_profile/result/gene_count_matrix.csv};

		# run mRNA with circRNA co_expression
		system qq{$rscript Rscript  $util/extract_diff_count.R $control_samples $case_samples $mRNA_count $mRNA_xls $result/$base_out/mRNA.count.xls};
		system qq{perl $util/extract_profile_by_name.pl $circRNA_xls $all_samples >$result/$base_out/circRNA.count.xls};
		system qq{perl $util/submit_co_expression.pl $result/$base_out/circRNA.count.xls $result/$base_out/mRNA.count.xls $result/$base_out/ circRNA mRNA  $util/co_expression.R $rscript Rscript} if not -e qq{$result/$base_out/correlation.xls};
	
		my %diff_miRNA  = ();
		my %diff_circRNA = ();
		my %diff_mRNA   = ();

		open EXO, $mRNA_xls or die "Can't open $mRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_mRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		open EXO, $miRNA_xls or die "Can't open $miRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
                        next if $arr[$#arr] eq 'type';
			$diff_miRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		open EXO, $circRNA_xls or die "Can't open $circRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
			next if $arr[$#arr] eq 'type';
			$diff_circRNA{$arr[0]} = [$arr[$#arr - 3], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		################################## report ################################	
		my $excel    = qq{$report/$base_out/CircRNA_with_mRNA.ceRNA.Network.xlsx};
		my $workbook = Excel::Writer::XLSX->new($excel);
		my %format   = package::format::run($workbook);

		my %circRNA_with_miRNA = ();
		my %mRNA_with_miRNA    = ();
		my %circcRNA_with_mRNA = ();

		# add worksheet "CircRNA_with_miRNA_interaction"
		my $worksheet = $workbook->add_worksheet("CircRNA_with_miRNA_Interaction");
		my $row = 0;
		my $cnt = 0;
		open DB, qq{$circRNA_target/$base_out/miranda.result.xls} or die "Can't open $circRNA_target/$base_out/miranda.result.xls!\n";
		while (<DB>) {
			chomp;
			my @arr = split /\t/;
			if ($cnt == 0 ) {
				$worksheet->write_row( $row, 0, \@arr, $format{'title'});
			} else {			
				my $circRNA_id = $arr[0];
				$circRNA_id =~ s/\|/-/;
				next if not exists $diff_circRNA{$circRNA_id};			
				$circRNA_with_miRNA{$circRNA_id}{$arr[1]} = "-";
				$worksheet->write_row( $row, 0, \@arr, $format{'normal'});
			}
			$cnt++;
			$row++;
		}
		close DB;

		# add worksheet "mRNA_with_miRNA_interaction"
		my $worksheet1 = $workbook->add_worksheet("mRNA_with_miRNA_Interaction");
		my $row = 0;
		my $cnt = 0;
		open DB, qq{$mRNA_with_miRNA_db} or die "Can't open $mRNA_with_miRNA_db!\n";
		while (<DB>) {
			chomp;
			my @arr = split /\t/;
			if ($cnt == 0 ) {
				$worksheet1->write_row( $row, 0, \@arr, $format{'title'});
			} else {
				# next if not exists $diff_miRNA{$arr[2]};
				next if not exists $diff_mRNA{$arr[0]};
				$mRNA_with_miRNA{$arr[0]}{$arr[2]} = "-";
				$worksheet1->write_row( $row, 0, \@arr, $format{'normal'});
			}
			$cnt++;
			$row++;
		}
		close DB;

		# add worksheet "CircRNA_with_mRNA_Coexpresion"
		my $worksheet2 = $workbook->add_worksheet("CircRNA_with_mRNA_Coexpresion");
		my $co_expression_xls = qq{$result/$base_out/correlation.xls};
		my $row = 0;
		open CORRELATION, qq{$co_expression_xls} or die "Can't open $co_expression_xls!\n";
		while (<CORRELATION>) {
			chomp;
			my @arr = split /\t/;
			if ($row == 0) {
				$worksheet2->write_row($row, 0, \@arr, $format{'title'});
			} else {
				next if $arr[2] < 0;
				next if $arr[3] > 0.05;

				$circRNA_with_mRNA{$arr[1]}{$arr[0]} = [$arr[2], $arr[3]];
				$worksheet2->write_row($row, 0, \@arr, $format{'normal'});
			}
			$row++;
		}
		close CORRELATION;

		# 生成ceRNA.xls
		open CERNA, qq{>$result/$base_out/ceRNA.xls};

		my %common_miRNA = ();
		my %common_cnt   = ();

		foreach my $x (keys %mRNA_with_miRNA) {
			my $mRNA_target_cnt = scalar keys %{$mRNA_with_miRNA{$x}};
			foreach my $y (keys %{$circRNA_with_mRNA{$x}}) {
				next if not exists $circRNA_with_miRNA{$y};
				my $lncRNA_target_cnt = scalar keys  %{$circRNA_with_miRNA{$y}};
				my @mRNA_targets      = keys %{$mRNA_with_miRNA{$x}};
				my @lncRNA_targets    = keys %{$circRNA_with_miRNA{$y}};
				# cal overlap with two sets;
				my %cnt = ();
				map { $cnt{$_}++; } (@mRNA_targets, @lncRNA_targets);
				my @common_targets  = grep {$cnt{$_} == 2} keys %cnt;

				foreach my $miRNA (@common_targets) {
					push @{$common_miRNA{$miRNA}}, qq{$x\t$y};
					$common_cnt{$miRNA}++;
				}
				my $common  = scalar @common_targets;
				my @diff_common_miRNAs = ();

				foreach my $miRNA (@common_targets) {					
					push @diff_common_miRNAs, $miRNA if exists $diff_miRNA{$miRNA};
				}

				my $common_miRNAs = join ",", @diff_common_miRNAs;

				next if $common == 0;
				next if $common_miRNAs eq '';
				
				my $ceRNA_score = $common / ($mRNA_target_cnt + $lncRNA_target_cnt -  $common);
				my $cor         = $circRNA_with_mRNA{$x}{$y}->[0];
				my $cor_pvalue  = $circRNA_with_mRNA{$x}{$y}->[1];

				my @res  = ($x, qq{protein_coding}, $diff_mRNA{$x}->[2], $y, qq{circRNA}, $diff_circRNA{$y}->[2], $common, $mRNA_target_cnt, $lncRNA_target_cnt, $ceRNA_score, "NA", $cor, $cor_pvalue);
				my $line = join "\t", @res;

				print CERNA qq{$line\n};
			}
		}

		close CERNA;

		# add worksheet "ceRNA_Network"
		my $worksheet3  = $workbook->add_worksheet("ceRNA_Network");
		my @head = qw/geneName geneType geneRegulation cernaName cernaType cernaRegulation commonMirNum geneMirNum cernaMirNum ceRNA_score pvalue FDR cor cor_pvalue/;
		$worksheet3->write_row(0, 0, \@head, $format{'title'});
		system qq{$rscript Rscript $util/dphyer_test.R  -s $organ $result/$base_out/ceRNA.xls $result/$base_out/ceRNA.pvalue.xls};
		open RESULT, qq{$result/$base_out/ceRNA.pvalue.xls} or die "Can't open $result/$base_out/ceRNA.pvalue.xls!\n";
		my $row = 1;
		while (<RESULT>) {
			chomp;
			my @arr = split /\t/;
			$worksheet3->write_row($row, 0, \@arr, $format{'normal'});
			$row++;
		}
		close RESULT;

		$workbook->close();

		# 生成edge.tx和node.txt
		open EDGE, qq{>$report/$base_out/edge.txt} or die "Can't open $report/$base_out/edge.txt!\n";	
		my %nodes = ();
		my $cnt = 0;
		my %unique = ();
		foreach my $x (sort {$common_cnt{$b} <=> $common_cnt{$a}} keys %common_cnt) {
			next if not exists $diff_miRNA{$x};
			$cnt++;
			next if $cnt > 5;
			foreach my $y (@{$common_miRNA{$x}}) {
				my ($mRNA, $lncRNA) = split /\t/, $y;
				print EDGE qq{$lncRNA\t$x\n} if not exists $unique{qq{$lncRNA:$x}};
				print EDGE qq{$mRNA\t$x\n}   if not exists $unique{qq{$mRNA:$x}};

				$unique{qq{$lncRNA:$x}} = "-";
				$unique{qq{$mRNA:$x}}   = "-";			

				$nodes{$lncRNA} = "circRNA";
				$nodes{$mRNA}   = "mRNA";
				$nodes{$x}      = "miRNA";
				
			}
		}
		close EDGE;

		open NODE, qq{>$report/$base_out/node.txt} or die "Can't open $report/$base_out/node.txt!\n";
		foreach my $x (keys %nodes) {

			print NODE qq{$x\tcircRNA\t$diff_circRNA{$x}->[2]\n} if $nodes{$x} eq 'circRNA';
			print NODE qq{$x\tmRNA\t$diff_mRNA{$x}->[2]\n}     if $nodes{$x} eq 'mRNA';
			print NODE qq{$x\tmiRNA\t$diff_miRNA{$x}->[2]\n}   if $nodes{$x} eq 'miRNA';

		}
		close NODE;  

		# cytoscape html
		system qq{perl $util/cytoscape_js.pl -node $report/$base_out/node.txt -edge $report/$base_out/edge.txt -o $report/$base_out};

	}

	

	print qq{mRNA miRNA circRNA ceRNA已经分析完成!\n};

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
		next if -e qq{$out/$base_out/CircRNA_with_mRNA.ceRNA.Network.xlsx};	
		push @temp, $x;
	}
	return @temp;
}


1;

