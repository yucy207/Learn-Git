package package::lncRNA::lncRNA_target_deseq;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $organ    = qq{$metadata->{organ}};

	my $result   = qq{$metadata->{project}/lncRNA/target_predict};
	my $report   = qq{$metadata->{report}/04_lncRNA_Analysis/03_Target_Gene_Prediction};
	my $deseq    = qq{$metadata->{project}/lncRNA/deseq/result};
	
	my $util           = qq{$base->{util}};
	my $gtf_to_fasta   = qq{$base->{lncRNA_gtf_to_fasta_bin}};
	my $ref_lncRNA_gtf = qq{$base->{$organ}{genome_lncRNA_one_gtf}};
	my $ref_mRNA_gtf   = qq{$base->{$organ}{genome_mRNA_one_gtf}};
	my $ref_fasta      = qq{$base->{$organ}{genome_fasta}};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;

	my @res_groups = res_check( qq{$result/result}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{lncRNA 靶基因预测已经分析完成!\n};
		#return;
	}

	pre_check($metadata, $base);

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $pid = $pm->start and next;

		my $control  = $x->[0];
		my $case     = $x->[1];
		my $base_out = qq{$case\_vs_$control};

		my $de_temp_dir  = qq{$deseq/$base_out};
		my $report_dir   = qq{$report/$base_out};
		my $lncRNA_deseq = qq{$de_temp_dir/gene/diff.xls};

		my $cmd    = qq{perl $util/lncTar_run_subimit_by_gene.pl $lncRNA_deseq $gtf_to_fasta $ref_lncRNA_gtf $ref_fasta $ref_mRNA_gtf $result/result/$base_out};
		my $cat    = qq{cat $result/result/$base_out/*/*target.txt | grep -v "Query"> $result/result/$base_out/lncRNA.target.xls};
		my $title  = qq{(echo -e  "Query\\tLength_Query\\tTarget\\tLength_Target\\tdG\\tndG\\tStart_Position_Query\\tEnd_Position_Query\\tStart_Position_Target\\tEnd_Position_Target";awk '\$6 < -0.1' $result/result/$base_out/lncRNA.target.xls) > $result/result/$base_out/lncRNA.target.fmt.xls};
		my $finish = qq{touch $result/result/$base_out/$case\_vs_$control.finish};

		open  SAVE, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $result/run/$case\_vs_$control.sh!\n";
		print SAVE qq{$cmd\n};
		print SAVE qq{$cat\n};
		print SAVE qq{$title\n$finish\n};
		close SAVE;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log\n};

		$pm->finish;
	}

	$pm->wait_all_children;

	#################### report:lncRNA_Targets.xlsx ################
	my $excel    = qq{$report/Differential_Expression_lncRNA_Targets.xlsx};
	my $workbook = Excel::Writer::XLSX->new($excel);
	my %format   = package::format::run($workbook);

	foreach my $x (@groups) {

		my $control = $x->[0];
		my $case    = $x->[1];

		my $base_out  = qq{$case\_vs_$control};
		my $worksheet = $workbook->add_worksheet(qq{$case\_vs_$control});

		my $row = 0;
		open EXO, qq{$result/result/$base_out/lncRNA.target.fmt.xls};
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

	my $worksheet = $workbook->add_worksheet("README");
	$worksheet->write_row( 0, 0, [decode("utf8", "标题"), decode("utf8", "说明")], $format{'title'});

	my @readme = ("Query",                 "查询序列（lncRNA）",
                  "Length_Query",          "查询序列长度",
                  "Target",                "目标序列（基因）",
                  "Length_Target",         "目标序列长度",
                  "dG",                    "结合自由能",
                  "ndG",                   "结合自由能（归一化，除以序列长度）",
                  "Start_Position_Query",  "查询序列比对起始位置",
                  "End_Position_Query",    "查询序列比对结束位置",
                  "Start_Position_Target", "目标序列比对起始位置",
                  "End_Position_Target",   "目标序列比对结束位置"
				);
	my $i = 0;
	foreach my $row(1..10){

		$worksheet->write_row($row, 0, [decode("utf8", $readme[$i]), decode("utf8", $readme[$i + 1])], $format{'normal'});
		$i = $i + 2;

	}
	$workbook->close();

	print qq{lncRNA 靶基因预测分析完成!\n};

}

sub pre_check
{

	my $metadata = shift;
	my $base     = shift;

	my $organ           = qq{$metadata->{organ}};
	my $ref_lncRNA_gtf  = qq{$base->{$organ}{genome_lncRNA_gtf}};
	my $ref_mRNA_gtf    = qq{$base->{$organ}{genome_mRNA_gtf}};
	my $ref_fasta       = qq{$base->{$organ}{genome_fasta}};
	# reference genome fiel is exists 
	die "reference genome fasta file is not exists!\n"      if not -e $ref_fasta;
	die "reference genome mRNA gtf file is not exists!\n"   if not -e $ref_mRNA_gtf;
	die "reference genome lncRNA gtf file is not exists!\n" if not -e $ref_lncRNA_gtf;

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

1;
