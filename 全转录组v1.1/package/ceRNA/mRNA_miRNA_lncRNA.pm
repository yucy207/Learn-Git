package package::ceRNA::mRNA_miRNA_lncRNA;
use package::format;
use Excel::Writer::XLSX;
use Parallel::ForkManager;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $organ    = qq{$metadata->{organ}};

	my $util     = qq{$base->{util}};
	my $rscript  = qq{$base->{rscript_bin}};

	my $organ1   = qq{$base->{$organ}{common_name}};
	my $lncRNA_with_miRNA_db = qq{$base->{$organ}{lncRNA_miRNA_db}};
	my $mRNA_with_miRNA_db   = qq{$base->{$organ}{mRNA_miRNA_db}};

	my $miRNA_result  = qq{$metadata->{miRNA_project}/Diff_Expression/result};
	my $mRNA_result   = qq{$metadata->{project}/mRNA/deseq/result};
	my $lncRNA_result = qq{$metadata->{project}/lncRNA/deseq/result};
	my $lncRNA_mRNA_coexpression = qq{$metadata->{project}/association/lncRNA_mRNA/Cis_Analysis};

	my $report = qq{$metadata->{report}/07_CeRNA_Analysis/lncRNA_miRNA_mRNA};
	my $result = qq{$metadata->{project}/ceRNA/lncRNA_miRNA_mRNA};

	system qq{mkdir -p $report} if not -d $report;
	system qq{mkdir -p $result} if not -d $result;

	my @res_groups = res_check( qq{$report}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{mRNA miRNA circRNA ceRNA已经分析完成!\n};
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

		my $mRNA_xls   = qq{$mRNA_result/$base_out/gene/diff.xls};
		my $miRNA_xls  = qq{$miRNA_result/$base_out/Known/sig_deseq_known.xls};
		my $lncRNA_xls = qq{$lncRNA_result/$base_out/gene/diff.xls};

		my %diff_miRNA  = ();
		my %diff_lncRNA = ();
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

		open EXO, $lncRNA_xls or die "Can't open $lncRNA_xls!\n";
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if $arr[$#arr] eq 'Not DEG';
            next if $arr[$#arr] eq 'type';
			$diff_lncRNA{$arr[0]} = [$arr[$#arr - 5], $arr[$#arr - 2], $arr[$#arr]];
		}
		close EXO;

		################################## report ################################		
		my $excel     = qq{$report/$base_out/LncRNA_with_mRNA.ceRNA.Network.xlsx};
		my $workbook  = Excel::Writer::XLSX->new($excel);
		my %format    = package::format::run($workbook);

		my %lncRNA_with_miRNA = ();
		my %mRNA_with_miRNA   = ();
		my %lncRNA_with_mRNA  = ();

		# add worksheet "LncRNA_with_miRNA_interaction"
		my $worksheet  = $workbook->add_worksheet("LncRNA_with_miRNA_Interaction");
		my $row = 0;
		my $cnt = 0;
		open DB, qq{$lncRNA_with_miRNA_db} or die "Can't open $lncRNA_with_miRNA_db!\n";
		while (<DB>) {
			chomp;
			my @arr = split /\t/;
			if ($cnt == 0 ) {
				$worksheet->write_row( $row, 0, \@arr, $format{'title'});

			} else {
				# next if not exists $diff_miRNA{$arr[1]};
				my ($lncRNA_id) = $arr[0] =~ /\((.+?)\)/;
				next if not exists $diff_lncRNA{$lncRNA_id};

				$lncRNA_with_miRNA{$lncRNA_id}{$arr[1]} = "-";

				$worksheet->write_row( $row, 0, \@arr, $format{'normal'});
			}
			$cnt++;
			$row++;
		}
		close DB;

		# add worksheet "mRNA_with_miRNA_interaction"
		my $worksheet1  = $workbook->add_worksheet("mRNA_with_miRNA_Interaction");
		my $row =  0;
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

		# add worksheet "LncRNA_with_mRNA_Coexpresion"
		my $worksheet2  = $workbook->add_worksheet("LncRNA_with_mRNA_Coexpresion");
		my $co_expression_xls = qq{$lncRNA_mRNA_coexpression/$base_out/correlation.xls};
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

				$lncRNA_with_mRNA{$arr[1]}{$arr[0]} = [$arr[2], $arr[3]];
				$worksheet2->write_row($row, 0, \@arr, $format{'normal'});
			}
			$row++;

		}
		close CORRELATION;

		# 生成ceRNA.xls
		my %common_miRNA = ();
		my %common_cnt   = ();

		open CERNA, qq{>$result/$base_out/ceRNA.xls};

		foreach my $x (keys %mRNA_with_miRNA) {
			my $mRNA_target_cnt = scalar keys %{$mRNA_with_miRNA{$x}};
			foreach my $y (keys %{$lncRNA_with_mRNA{$x}}) {
				next if not exists $lncRNA_with_miRNA{$y};
				my $lncRNA_target_cnt = scalar keys  %{$lncRNA_with_miRNA{$y}};
				my @mRNA_targets      = keys %{$mRNA_with_miRNA{$x}};
				my @lncRNA_targets    = keys %{$lncRNA_with_miRNA{$y}};
				# cal overlap with two sets;
				my %cnt = ();
				map { $cnt{$_}++; } (@mRNA_targets, @lncRNA_targets);
				my @common_targets  = grep {$cnt{$_} == 2} keys %cnt;

				foreach my $miRNA (@common_targets) {
					push @{$common_miRNA{$miRNA}}, qq{$x:$y};
					$common_cnt{$miRNA}++;
				}

				my $common = scalar @common_targets;

				my @diff_common_miRNAs = ();

				foreach my $miRNA (@common_targets) {					
					push @diff_common_miRNAs, $miRNA if exists $diff_miRNA{$miRNA};
				}
				
				my $common_miRNAs = join ",", @diff_common_miRNAs;

				next if $common == 0;
				next if $common_miRNAs eq '';
				
				my $ceRNA_score = $common / ($mRNA_target_cnt + $lncRNA_target_cnt -  $common);
				my $cor         = $lncRNA_with_mRNA{$x}{$y}->[0];
				my $cor_pvalue  = $lncRNA_with_mRNA{$x}{$y}->[1];

				my @res  = ($x, qq{protein_coding}, $diff_mRNA{$x}->[2], $y, qq{lncRNA}, $diff_lncRNA{$y}->[2], $common, $mRNA_target_cnt, $lncRNA_target_cnt, $ceRNA_score, "NA", $cor, $cor_pvalue);
				my $line = join "\t", @res;

				print CERNA  qq{$line\n};

			}
		}

		close CERNA;

		# add worksheet "ceRNA_Network"
		my $worksheet3 = $workbook->add_worksheet("ceRNA_Network");
		my @head = qw/geneName geneType geneRegulation cernaName cernaType cernaRegulation commonMirNum geneMirNum cernaMirNum ceRNA_score pvalue FDR cor cor_pvalue/;
		$worksheet3->write_row(0, 0, \@head, $format{'title'});
		system qq{$rscript Rscript $util/dphyer_test.R -s $organ1 $result/$base_out/ceRNA.xls $result/$base_out/ceRNA.pvalue.xls};
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
				my ($mRNA, $lncRNA) = split /:/, $y;
				print EDGE qq{$lncRNA\t$x\n} if not exists $unique{qq{$lncRNA:$x}};
				print EDGE qq{$mRNA\t$x\n}   if not exists $unique{qq{$mRNA:$x}};

				$unique{qq{$lncRNA:$x}} = "-";
				$unique{qq{$mRNA:$x}}   = "-";

				$nodes{$lncRNA} = "lncRNA";
				$nodes{$mRNA}   = "mRNA";
				$nodes{$x}      = "miRNA";			
			}
		}
		close EDGE;

		open NODE, qq{>$report/$base_out/node.txt} or die "Can't open $report/$base_out/node.txt!\n";
		foreach my $x (keys %nodes) {

			print NODE qq{$x\tlncRNA\t$diff_lncRNA{$x}->[2]\n} if $nodes{$x} eq 'lncRNA';
			print NODE qq{$x\tmRNA\t$diff_mRNA{$x}->[2]\n}     if $nodes{$x} eq 'mRNA';
			print NODE qq{$x\tmiRNA\t$diff_miRNA{$x}->[2]\n}   if $nodes{$x} eq 'miRNA';

		}
		close NODE;

		# cytoscape html
		system qq{perl $util/cytoscape_js.pl -node $report/$base_out/node.txt -edge $report/$base_out/edge.txt -o $report/$base_out};


		$pm->finish;
	}

	$pm->wait_all_children;

	print qq{mRNA miRNA lncRNA ceRNA 分析已经运行完成!\n};

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
		next if -e qq{$out/$base_out/LncRNA_with_mRNA.ceRNA.Network.xlsx};	
		push @temp, $x;
	}
	return @temp;
}


1;
