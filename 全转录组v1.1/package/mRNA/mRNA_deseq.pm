package package::mRNA::mRNA_deseq;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use Encode qw/decode/;
use package::format;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my @samples  = split /,/, $metadata->{'samples'};
	my $organ    = qq{$metadata->{organ}};

	my $expression   = qq{$metadata->{project}/mRNA/gene_profile/result};
	my $result       = qq{$metadata->{project}/mRNA/deseq};
	my $diff_report  = qq{$metadata->{report}/03_mRNA_Analysis/02_Differential_Expression_Analysis};

	my $log2fc       = qq{$base->{log2fc}};
	my $rscript      = qq{$base->{rscript_bin}};
	my $util         = qq{$base->{util}};
	my $pca_plot     = qq{$base->{pca_bin}};

	my $ref_mRNA_gtf = qq{$base->{$organ}{genome_mRNA_gtf}};	
	my $go_db        = qq{$base->{$organ}{go_info_db}};
	my $kegg_db      = qq{$base->{$organ}{kegg_info_db}};

	if (not exists $metadata->{'group'}) {

		print "[Warning] : You must set group in the config.txt!\n";
		exit;
	}

	pre_check($metadata, $base);
	
	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $diff_report}   if not -d $diff_report;

	my @res_groups = res_check(qq{$result/result}, \@groups);

	if (scalar @res_groups == 0) {

		print qq{mRNA 差异分析已经运行完成!\n};
		#return 0;
	}  

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $pid = $pm->start and next;

		my $control = $x->[0];
		my $case    = $x->[1];
		my @control_samples  = split /,/, (split /;/, $x->[2])[0];
		my @case_samples     = split /,/, (split /;/, $x->[2])[1];
		my $heatmap_case     = join ",",  @case_samples;
		my $heatmap_control  = join ",",  @control_samples;

		my $base_out         = qq{$case\_vs_$control};
		my $de_temp_dir      = qq{$result/result/$base_out};
		my $de_report_dir    = qq{$diff_report/$case\_vs_$control\_Differential_Expression_mRNA_Figures};

		system qq{mkdir -p $de_report_dir}          if not -d qq{$de_report_dir};
		system qq{mkdir -p $de_temp_dir/gene}       if not -d qq{$de_temp_dir/gene}; 
		system qq{mkdir -p $de_temp_dir/transcript} if not -d qq{$de_temp_dir/transcript};

		open OUT, qq{>$de_temp_dir/gene/sample.groups} or die "Can't make $de_temp_dir/gene/sample.groups!\n";
		foreach my $x (@control_samples) {
			print OUT qq{$x\t$control\n};
		}
		foreach my $y (@case_samples) {
			print OUT qq{$y\t$case\n};
		}

		# DESeq2
		my $gene_cmd           = qq{$rscript Rscript $util/Deseq2_diff.R $expression/gene_count_matrix.csv $case $control $heatmap_case $heatmap_control $log2fc $de_temp_dir/gene};
		my $extract_diff_count = qq{perl $util/extract_diff_count.pl $de_temp_dir/gene/diff.xls $expression/count.rlog.xls $de_temp_dir/gene/de_count.rlog.xls};
		my $diff_pca_cmd       = qq{$pca_plot Rscript $util/pca_plot.R $de_temp_dir/gene/de_count.rlog.xls $de_temp_dir/gene/sample.groups $de_temp_dir/gene/ &> $de_temp_dir/gene/pca_plot.log};
		my $trans_cmd          = qq{$rscript Rscript $util/Deseq2_diff.R $expression/transcript_count_matrix.csv $case $control $heatmap_case $heatmap_control $log2fc  $de_temp_dir/transcript/};
		my $trans_fmt          = qq{perl $util/deseq_trans_fmt.pl $ref_mRNA_gtf $de_temp_dir/transcript/diff.xls > $de_temp_dir/transcript/trans.diff.xls};
		# report pdf
		my $cp = qq{cp $de_temp_dir/gene/MA.pdf $de_temp_dir/gene/heatmap.pdf $de_temp_dir/gene/Top50.heatmap.pdf $de_temp_dir/gene/valcano.pdf $de_temp_dir/gene/correlation.pdf $de_temp_dir/gene/pca.pdf $de_report_dir};	
		# finish
		my $touch              = qq{touch $de_temp_dir/$case\_vs_$control.finish};

		open SAVE, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $result/run/$case\_vs_$control.sh!\n";
		print SAVE qq{$gene_cmd\n$extract_diff_count\n$diff_pca_cmd\n$trans_cmd\n$cp\n$trans_fmt\n};
		print SAVE qq{$touch};
		close SAVE;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log\n};
		$pm->finish;

	}

	$pm->wait_all_children;

	##################### read go annotaion info #################
	my %go_info = ();
	open GODB, $go_db;
	while (<GODB>) {
		chomp;
		my @arr = split /\t/;
		$go_info{$arr[0]}{'bp'} = $arr[1];
		$go_info{$arr[0]}{'mf'} = $arr[2];
		$go_info{$arr[0]}{'cc'} = $arr[3];
	}
	close GODB;
	#################### read kegg annotation info ################
	my %kegg_info = ();
	open KEGGDB, $kegg_db;
	while (<KEGGDB>) {
		chomp;
		my @arr = split /\t/;
		$kegg_info{$arr[1]} = $arr[2];
	}
	close KEGGDB;

	#################### report:Sample.groups.xlsx ################
	my $excel     = qq{$diff_report/Sample.groups.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);
	my $worksheet = $workbook->add_worksheet("Summary");

	my %unique = ();
	foreach my $x (@groups) {
		my $control = $x->[0];
		my $case    = $x->[1];

		my $control_samples = (split /;/, $x->[2])[0];
		my $case_samples    = (split /;/, $x->[2])[1];

		$unique{$control} = $control_samples;
		$unique{$case}    = $case_samples;
	}

	my @head   = ("Group", "Samples");
	$worksheet->write_row( 0, 0, \@head, $format{'title'});
	my $row = 1;
	foreach my $x (keys %unique) {
		my @line = ($x,  $unique{$x});
		$worksheet->write_row( $row, 0, \@line, $format{'normal'});
		$row++;
	}

	$workbook->close();
	
	################### report:different gene summary xlsx ##############
	my $excel     = qq{$diff_report/Differential_Expression_Genes_Summary.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);
	report_xlsx($workbook, \@groups, "gene", \%go_info, \%kegg_info, \%format, $result);
	$workbook->close();


	############### report:different transcript summary xlsx ############
	my $excel     = qq{$diff_report/Differential_Expression_Transcripts_Summary.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);
	report_xlsx($workbook, \@groups, "transcript", \%go_info, \%kegg_info, \%format, $result);
	$workbook->close();

	print qq{mRNA 差异分析已经运行完成!\n};


}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;

	my $organ    = qq{$metadata->{organ}};
	my $ref_mRNA_gtf  = qq{$base->{$organ}{genome_mRNA_gtf}};
	# reference file is exists
	die "reference mRNA gtf file is not exists!\n" if not -e $ref_mRNA_gtf;	
	# group config is exists 
	my @groups = @{$metadata->{'group'}};
	die "group config info is not exists!\n" if scalar @groups == 0;

	foreach my $x (@groups) {
		my $control = $x->[0];
		my $case    = $x->[1];

		my @control_samples = split /,/, (split /;/, $x->[2])[0];
		my @case_samples    = split /,/, (split /;/, $x->[2])[1];

		die qq{$control must have at least one sample!\n} if scalar @control_samples == 0;
		die qq{$case    must have at least one sample!\n} if scalar @case_samples == 0;
	}	

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

		next if  -e qq{$out/$base_out/$case\_vs_$control.finish};
		push @temp, $x;
	}
	return @temp;
}

sub report_xlsx
{
	my ($workbook, $groups, $type, $go_info, $kegg_info, $format, $result) = @_;

	my $worksheet = $workbook->add_worksheet("Summary");	
	my $title     = qq{Grpup\tAll\tUp\tDown};
	my @head      = split /\t/, $title;
	$worksheet->write_row( 0, 0, \@head, $format->{'title'});

	my $row = 1;
	foreach my $x (@$groups) {
		my $control = $x->[0];
		my $case    = $x->[1];

		my $base_out    = qq{$case\_vs_$control};
		my $de_temp_dir = qq{$result/result/$base_out};

		open EXO, qq{$de_temp_dir/$type/diff.xls};
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

		$worksheet->write_row( $row, 0, \@result, $format->{'normal'});
		$row++;
	}

	foreach my $x (@$groups) {
		my $control = $x->[0];
		my $case    = $x->[1];

		my $base_out    = qq{$case\_vs_$control};
		my $de_temp_dir = qq{$result/result/$base_out};
		my $worksheet   = $workbook->add_worksheet(qq{$case\_vs_$control});

		my $row = 0;
		if($type eq "gene"){open EXO, qq{$de_temp_dir/$type/diff.xls};}
		if($type eq "transcript"){open EXO, qq{$de_temp_dir/$type/trans.diff.xls};}
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			my ($go_bp, $go_mf, $go_cc, $kegg);
			if($type eq "gene"){

				$go_bp = $go_info->{$arr[0]}{'bp'};
				$go_mf = $go_info->{$arr[0]}{'mf'};
				$go_cc = $go_info->{$arr[0]}{'cc'};
				$kegg  = $kegg_info->{$arr[0]};

			}
			if($type eq "transcript"){

				$go_bp = $go_info->{$arr[1]}{'bp'};
				$go_mf = $go_info->{$arr[1]}{'mf'};
				$go_cc = $go_info->{$arr[1]}{'cc'};
				$kegg  = $kegg_info->{$arr[1]};

			}
			push @arr, $go_bp;
			push @arr, $go_mf;
			push @arr, $go_cc;
			push @arr, $kegg;

			if ($row == 0) {

				if($type eq "gene"){$arr[0] = "gene_id";}
				if($type eq "transcript"){$arr[0] = "transcript_id";$arr[1] = "gene_id";}
				$arr[$#arr - 3 ] = "GO_BP";
				$arr[$#arr - 2] = "GO_MF";
				$arr[$#arr - 1] = "GO_CC";
				$arr[$#arr] = "KEGG_PATHWAY";	

				$worksheet->write_row( $row, 0, \@arr, $format->{'title'});	
			} else {

				$worksheet->write_row( $row, 0, \@arr, $format->{'normal'});
			}
			$row++;
		}
		close EXO;
	}

	my $worksheet = $workbook->add_worksheet("README");
	$worksheet->write_row( 0, 0, [decode("utf8", "标题"), decode("utf8", "说明")], $format->{'title'});

	my @readme = ("sample_name",       "归一化之后的表达量",
                  "baseMeanA",         "control组中的平均表达量",
                  "baseMeanB",         "case组中的平均表达量",
                  "baseMean",          "在所有样本中的平均表达量",
                  "log2(fold_change)", "case组的表达量比上control组中的表达量",
                  "lfcSE",             "标准误的对数值",
                  "stat",              "假设检验的统计量值",
                  "pvalue",            "P值",
                  "padj",              "校正之后的P值（BH方法校正后的P值）",
                  "type",              "P值小于0.05，并且根据log2(fold_change)给出的差异类型，大于1为Up，小于-1为Down",
                  "GO_BP",             "基因对应的GO数据库BP类别的注释信息",
                  "GO_MF",             "基因对应的GO数据库MF类别的注释信息",
                  "GO_CC",             "基因对应的GO数据库CC类别的注释信息",
                  "KEGG_PATHWAY",      "基因对应的kegg数据库的pathway信息"
				);

	if($type eq "gene"){
	
		my $B1 = decode("utf8", "gene_id"); my $B2 = decode("utf8", "基因ID");	
		$worksheet->write_row( 1, 0, [$B1, $B2], $format->{'normal'});
		my $i = 0;
		foreach my $row(2..15){
			$worksheet->write_row($row, 0, [decode("utf8", $readme[$i]), decode("utf8", $readme[$i + 1])], $format->{'normal'});
			$i = $i + 2;
		}
	}

	if($type eq "transcript"){

		my $B1 = decode("utf8", "transcript_id"); my $B2 = decode("utf8", "转录本ID");
		my $C1 = decode("utf8", "gene_id");       my $C2 = decode("utf8", "转录本对应的基因ID");
		$worksheet->write_row( 1, 0, [$B1, $B2], $format->{'normal'});
		$worksheet->write_row( 2, 0, [$C1, $C2], $format->{'normal'});
		my $i = 0;
		foreach my $row(3..16){
			$worksheet->write_row($row, 0, [decode("utf8", $readme[$i]), decode("utf8", $readme[$i + 1])], $format->{'normal'});
			$i = $i + 2;
		}
	}

}


1;
