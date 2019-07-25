package package::circRNA::circRNA_diff;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $organ1   = qq{$metadata->{organ}};
	my $organ    = qq{$base->{$organ1}{common_name}};

	my $rscript  = qq{$base->{rscript_bin}};
	my $util     = qq{$base->{util}};
	my $log2fc   = qq{$base->{log2fc_circ}};

	my $profile  = qq{$metadata->{project}/circRNA/circrna_profile/result/circRNA.expression.xls};
	my $profile_anno  = qq{$metadata->{project}/circRNA/circrna_profile/result/circRNA.expression.add.annotation.xls};
	my $Summary  = qq{$metadata->{project}/circRNA/circRNA_predict/result/Summary.xls};
	my $result   = qq{$metadata->{project}/circRNA/different_analysis};
	my $report   = qq{$metadata->{report}/05_circRNA_Analysis/03_Differential_Expression_Analysis};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;

	my @res_groups = res_check(qq{$result/result},  \@groups);

	if (scalar @res_groups == 0) {
		print qq{circRNA 差异分析已经运行完成!\n};
		return;
	}

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $pid = $pm->start and next;

		my $control = $x->[0];
		my $case    = $x->[1];
		my @control_samples = split /,/, (split /;/, $x->[2])[0];
		my @case_samples    = split /,/, (split /;/, $x->[2])[1];
		my $control_name    = join "," , @control_samples;
		my $case_name       = join "," , @case_samples;

		my $base_out        = qq{$case\_vs_$control};
		my $de_out_dir      = qq{$result/result/$base_out/gene};
		my $de_report_dir   = qq{$report/$case\_vs_$control};

		system qq{mkdir -p $de_out_dir}    if not -d qq{$de_out_dir};
		system qq{mkdir -p $de_report_dir} if not -d qq{$de_report_dir};
		# DEseq
		my $cmd           = qq{$rscript Rscript $util/Deseq_diff.R $log2fc $case_name $control_name $case $control $profile $de_out_dir/circRNA.different.expression.xls};
		# MA
		my $ma            = qq{$rscript Rscript $util/circRNA_MA.R $de_out_dir/circRNA.different.expression.xls $de_report_dir/MA.pdf};
		# valcano
		my $valcano       = qq{$rscript Rscript $util/circRNA_valcano.R $de_out_dir/circRNA.different.expression.xls $de_report_dir/valcano.pdf};
		# correlation
		my $correlation   = qq{$rscript Rscript $util/circRNA_corrleation.R $de_out_dir/circRNA.different.expression.xls $de_report_dir/correlation.pdf};	
		# heatmap
		my $heatmap       = qq{$rscript Rscript $util/circRNA_heatmap.R $de_out_dir/circRNA.different.expression.xls $de_report_dir/heatmap.pdf};
		# Top50.heatmap
		my $Top50_heatmap = qq{$rscript Rscript $util/Top50.circRNA_heatmap.R $de_out_dir/circRNA.different.expression.xls $de_report_dir/Top50.heatmap.pdf};	
		#diff_ano
		my $add           = qq{perl $util/add_circRNA_diff_annotation.pl $Summary $profile_anno $de_out_dir/circRNA.different.expression.xls $organ $de_out_dir/circRNA.different.expression.add.annotation.xls};
		# finish
		my $finish        = qq{touch  $result/result/$base_out/gene/$case\_vs_$control.finish};

		open SAVE, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $result/run/$case\_vs_$control.sh!\n";
		print SAVE qq{$cmd\n$ma\n$valcano\n$correlation\n$heatmap\n$Top50_heatmap\n$add\n$finish};		
		close SAVE;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log\n};

		$pm->finish;
	}

	$pm->wait_all_children;

	############################### report #################################
	
	my $excel     = qq{$report/Differential_ Expression_circRNA_Summary.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);
	my $worksheet = $workbook->add_worksheet("Summary");

	my $title  = qq{Grpup\tAll\tUp\tDown};
	my @head   = split /\t/, $title;
	$worksheet->write_row( 0, 0, \@head, $format{'title'});

	my $row = 1;
	foreach my $x (@groups) {

		my $control = $x->[0];
		my $case    = $x->[1];
		my $sample  = $x->[2];
		my $base_out   = qq{$case\_vs_$control};
		my $de_out_dir = qq{$result/result/$base_out/gene};

		open EXO, qq{$de_out_dir/circRNA.different.expression.xls};
		my ($line_cnt, $up_cnt, $down_cnt) = (0, 0, 0, 0);
		while (<EXO>) {
			chomp;
			$line_cnt++;
			next if $line_cnt == 1;
			my @arr = split /\t/;
			$up_cnt++   if $arr[$#arr] eq 'Up';
			$down_cnt++ if $arr[$#arr] eq 'Down';	
		}
		close EXO;
		my $all = $up_cnt + $down_cnt;
		my $line = qq{$case\_vs_$control\t$all\t$up_cnt\t$down_cnt};
		my @result = split /\t/, $line;
		$worksheet->write_row( $row, 0, \@result, $format{'normal'});
		$row++;
	}

	foreach my $x (@groups) {

		my $control = $x->[0];
		my $case    = $x->[1];
		my $sample  = $x->[2];
		my $base_out   = qq{$case\_vs_$control};
		my $de_out_dir = qq{$result/result/$base_out/gene};
		
		my $worksheet = $workbook->add_worksheet(qq{$case\_vs_$control});

		my $row = 0;
		open EXO, qq{$de_out_dir/circRNA.different.expression.add.annotation.xls};
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			if ($row == 0) {
				$worksheet->write_row( $row, 0, \@arr, $format{'title'});	
			} else {
				$worksheet->write_row( $row, 0, \@arr, $format{'normal'});
			}
			$row++;
		}
		close EXO;

	}

	my $worksheet1 = $workbook->add_worksheet("README");
	my $A1 = decode("utf8", "标题"); my $A2 = decode("utf8", "说明");
	my $B1 = decode("utf8", "gene"); my $B2 = decode("utf8", "基因");
	my $C1 = decode("utf8", "basemean"); my $C2 = decode("utf8", "总体表达量均值");
	my $D1 = decode("utf8", "basemeanA"); my $D2 = decode("utf8", "第一组表达量均值");
	my $E1 = decode("utf8", "basemeanB"); my $E2 = decode("utf8", "第二组表达量均值");
    my $F1 = decode("utf8", "fold_change"); my $F2 = decode("utf8", "group_2/group1的表达量比值");
	my $G1 = decode("utf8", "log2FoldChange"); my $G2 = decode("utf8", "group_2/group1的表达量比值取log2的结果（inf 表示group_1的value为0，-inf表示group_2的value为0）");
	my $H1 = decode("utf8", "pval"); my $H2 = decode("utf8", "P值");
	my $I1 = decode("utf8", "pdaj"); my $I2 = decode("utf8", "P值校正（BH方法校正后的P值）");
	my $J1 = decode("utf8", "type"); my $J2 = decode("utf8", "P值小于0.05，并且根据log2(fold_change)给出的差异类型，大于0为Up，小于0为Down");
	my $K1 = decode("utf8", "circBase_ID"); my $K2 = decode("utf8", "circBase数据库ID");
	my $L1 = decode("utf8", "chr"); my $L2 = decode("utf8", "染色体");
	my $M1 = decode("utf8", "start"); my $M2 = decode("utf8", "起始位置");
	my $N1 = decode("utf8", "end"); my $N2 = decode("utf8", "终止位置");
	my $O1 = decode("utf8", "strand"); my $O2 = decode("utf8", "链的正负");
	my $P1 = decode("utf8", "genomic length"); my $P2 = decode("utf8", "基因组上的长度");
	my $Q1 = decode("utf8", "spliced length"); my $Q2 = decode("utf8", "剪切后转录本长度");
	my $R1 = decode("utf8", "circRNA_type"); my $R2 = decode("utf8", "circRNA类型");
	my $S1 = decode("utf8", "best_transcript"); my $S2 = decode("utf8", "最佳转录本");
	my $T1 = decode("utf8", "IRESfinder_Index"); my $T2 = decode("utf8", "IRESfinder软件里是否存在");
	my $U1 = decode("utf8", "IRESfinder_Score"); my $U2 = decode("utf8", "IRESfinder软件打分");

	my @Organ = ("human","mouse","fly","zebrafish");
	my ($V1, $V2, $W1, $W2, $X1, $X2, $Y1, $Y2);
	if(grep {$_ eq $organ} @Organ){
		$V1 = decode("utf8", "ORF_size"); $V2 = decode("utf8", "开放阅读框大小");
		$W1 = decode("utf8", "Fickett_score"); $W2 = decode("utf8", "Fickett评分");
		$X1 = decode("utf8", "Hexamer_score"); $X2 = decode("utf8", "Hexamer评分");
		$Y1 = decode("utf8", "coding_prob"); $Y2 = decode("utf8", "编码可能性");
	}

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
	$worksheet1->write_row( 11, 0, [$L1, $L2], $format{'normal'});
	$worksheet1->write_row( 12, 0, [$M1, $M2], $format{'normal'});
	$worksheet1->write_row( 13, 0, [$N1, $N2], $format{'normal'});
	$worksheet1->write_row( 14, 0, [$O1, $O2], $format{'normal'});
	$worksheet1->write_row( 15, 0, [$P1, $P2], $format{'normal'});
	$worksheet1->write_row( 16, 0, [$Q1, $Q2], $format{'normal'});
	$worksheet1->write_row( 17, 0, [$R1, $R2], $format{'normal'});
	$worksheet1->write_row( 18, 0, [$S1, $S2], $format{'normal'});
	$worksheet1->write_row( 19, 0, [$T1, $T2], $format{'normal'});
	$worksheet1->write_row( 20, 0, [$U1, $U2], $format{'normal'});

	if(grep {$_ eq $organ} @Organ){

		$worksheet1->write_row( 21, 0, [$V1, $V2], $format{'normal'});
		$worksheet1->write_row( 22, 0, [$W1, $W2], $format{'normal'});
		$worksheet1->write_row( 23, 0, [$X1, $X2], $format{'normal'});
		$worksheet1->write_row( 24, 0, [$Y1, $Y2], $format{'normal'});
	}

	$workbook->close();

	print qq{circRNA 差异分析运行完成!\n};

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
		next if  -e qq{$out/$base_out/gene/$case\_vs_$control.finish};
		push @temp, $x;
	}
	return @temp;
}

1;
