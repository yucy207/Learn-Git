package package::lncRNA::lncRNA_overlap_target_deseq;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'overlap'}};
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
		print qq{lncRNA 差异overlap靶基因预测已经分析完成!\n};
		#return;
	}

	pre_check($metadata, $base);

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_groups) {

		my $pid = $pm->start and next;

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @line  = ();
		my @names = ();

		foreach my $a (0..$#diffs) {

			my ($control, $case) = split /,/, $diffs[$a];
			my $diff = qq{$deseq/$case\_vs_$control/gene/diff.xls};

			system qq{grep -v "Not DEG" $diff | grep -v "pvalue" | cut -f1 | sort | uniq > $deseq/$case\_vs_$control/gene/de_genes.list};
			system qq{grep "Up" $diff | grep -v "pvalue" | cut -f1 | sort | uniq > $deseq/$case\_vs_$control/gene/up.list};
			system qq{grep "Down" $diff | grep -v "pvalue" | cut -f1 | sort | uniq > $deseq/$case\_vs_$control/gene/down.list};
			
			my $type      = $genes[$a];
			my %type2list = ('all' => 'de_genes.list',  'up' => 'up.list', 'down' => 'down.list');
			my $file      = qq{$deseq/$case\_vs_$control/gene/$type2list{$type}};
			push @names, qq{$case\_vs_$control\_$type};
			push @line , $file;
				
		}

		my $list = join " ", @line;
		my $name_list      = join "_with_", @names;
		my $target_tmp_dir = qq{$result/result/$name_list};
		system qq{mkdir -p $target_tmp_dir} if not -d $de_temp_dir;

		my $venn   = qq{$venndiagram Rscript $util/venn_plot.r $list $target_tmp_dir};
		my $cmd    = qq{perl $util/lncTar_run_subimit_by_overlap_gene.pl $target_tmp_dir/common_gene.tsv $gtf_to_fasta $ref_lncRNA_gtf $ref_fasta $ref_mRNA_gtf $target_tmp_dir};
		my $cat    = qq{cat $target_tmp_dir/*/*target.txt | grep -v "Query"> $target_tmp_dir/lncRNA.target.xls};
		my $title  = qq{(echo -e  "Query\\tLength_Query\\tTarget\\tLength_Target\\tdG\\tndG\\tStart_Position_Query\\tEnd_Position_Query\\tStart_Position_Target\\tEnd_Position_Target";awk '\$6 < -0.1' $target_tmp_dir/lncRNA.target.xls) > $target_tmp_dir/lncRNA.target.fmt.xls};
		my $finish = qq{touch $target_tmp_dir/$name_list.finish};

		open  SAVE, qq{>$result/run/$name_list.sh} or die "Can't open $result/run/$name_list.sh!\n";
		print SAVE qq{$venn\n};
		print SAVE qq{$cmd\n};
		print SAVE qq{$cat\n};
		print SAVE qq{$title\n$finish\n};
		close SAVE;

		system qq{bash $result/run/$name_list.sh &> $result/log/$name_list.log\n};

		$pm->finish;
	}

	$pm->wait_all_children;

	#################### report:overlap_lncRNA_Targets.xlsx ################
	my $excel    = qq{$report/Overlap_Differential_Expression_lncRNA_Targets.xlsx};
	my $workbook = Excel::Writer::XLSX->new($excel);
	my %format   = package::format::run($workbook);

	my $worksheet1 = $workbook->add_worksheet(qq{overlap_info});
	$worksheet1->write_row( 0, 0, [decode("utf8", "编号"), decode("utf8", "overlap组别信息")], $format{'title'});
	
	my $row1 = 1;
	foreach my $x (@groups) {

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @names = ();

		foreach my $a (0..$#diffs) {
			my ($control, $case) = split /,/, $diffs[$a];
			push @names, qq{$case\_vs_$control\_$genes[$a]};			
		}

		my $name_list = join "_with_", @names;
		$worksheet1->write_row( $row1, 0, [decode("utf8", "overlap\_$row1"), decode("utf8",  $name_list)], $format{'normal'});
	
		my $worksheet = $workbook->add_worksheet(qq{overlap\_$row1});
		$row1++;

		my $row = 0;
		open EXO, qq{$result/result/$name_list/lncRNA.target.fmt.xls};
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

	print qq{lncRNA 差异overlap靶基因预测分析完成!\n};

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

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @names = ();

		foreach my $a (0..$#diffs) {
			my ($control, $case) = split /,/, $diffs[$a];		
			push @names, qq{$case\_vs_$control\_$genes[$a]};			
		}

		my $name_list = join "_with_", @names;
		next if -e qq{$out/$name_list.finish};
		push @temp, $x;

	}
	return @temp;
}

1;
