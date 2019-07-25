package package::lncRNA::lncRNA_maSigPro_target_deseq;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @design   = @{$metadata->{'time'}};
	my $organ    = qq{$metadata->{organ}};

	my $result   = qq{$metadata->{project}/lncRNA/maSigPro_target_predict};
	my $report   = qq{$metadata->{report}/04_lncRNA_Analysis/03_Target_Gene_Prediction};
	my $maSigPro = qq{$metadata->{project}/lncRNA/maSigPro/result/gene};
	
	my $util           = qq{$base->{util}};
	my $gtf_to_fasta   = qq{$base->{lncRNA_gtf_to_fasta_bin}};
	my $ref_lncRNA_gtf = qq{$base->{$organ}{genome_lncRNA_one_gtf}};
	my $ref_mRNA_gtf   = qq{$base->{$organ}{genome_mRNA_one_gtf}};
	my $ref_fasta      = qq{$base->{$organ}{genome_fasta}};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;

	my @res = res_check(qq{$result/result}, \@design);

	if (scalar @res == 0) {
		print qq{lncRNA 时间序列靶基因预测已经分析完成!\n};
	}

	pre_check($metadata, $base);

	my $num = 1;
	foreach my $x (sort @res) {

		my $base_out = scalar @design == 1 ? '' : qq{time\_$num};
		my $tmp_dir  = qq{$result/result/$base_out};
		system qq{mkdir -p $tmp_dir} if not -d $tmp_dir;

		my $diff   = qq{$maSigPro/$base_out/all.diff.xls};		
		my $cmd    = qq{perl $util/lncTar_run_subimit_by_gene.pl $diff $gtf_to_fasta $ref_lncRNA_gtf $ref_fasta $ref_mRNA_gtf $tmp_dir};
		my $cat    = qq{cat $tmp_dir/*/*target.txt | grep -v "Query"> $tmp_dir/lncRNA.target.xls};
		my $title  = qq{(echo -e  "Query\\tLength_Query\\tTarget\\tLength_Target\\tdG\\tndG\\tStart_Position_Query\\tEnd_Position_Query\\tStart_Position_Target\\tEnd_Position_Target";awk '\$6 < -0.1' $tmp_dir/lncRNA.target.xls) > $tmp_dir/lncRNA.target.fmt.xls};
		my $finish = qq{touch $tmp_dir/target.finish};

		open  SAVE, qq{>$result/run/target$base_out.sh} or die "Can't open $result/run/target$base_out.sh!\n";
		print SAVE qq{$cmd\n};
		print SAVE qq{$cat\n};
		print SAVE qq{$title\n$finish\n};
		close SAVE;

		system qq{bash $result/run/target$base_out.sh &> $result/log/target$base_out.log\n};

		$num++;
	}

	####################### report:maSigPro_lncRNA_Targets.xlsx #####################
	my $excel    = qq{$report/MaSigPro_Differential_Expression_lncRNA_Targets.xlsx};
	my $workbook = Excel::Writer::XLSX->new($excel);
	my %format   = package::format::run($workbook);

	if(scalar @design != 1){

		my $worksheet1 = $workbook->add_worksheet(qq{MaSigPro_info});
		$worksheet1->write_row( 0, 0, [decode("utf8", "编号"), decode("utf8", "时间信息")], $format{'title'});
		
		my $num = 1;
		foreach my $x (sort @design) {
			my $base_out = qq{time\_$num};
			my %hash     = ();
			my @arr;
			open(IN,$x) or die "cannot open $x";
			while(<IN>){
				my($sample,$time,$replicate,$group) = split/\t/,$_;
				next if $time eq 'Time';
				if(not exists $hash{$time}){
					$hash{$time} = "-";
					push @arr,$time;
				}					
			}
			close IN;
			my $info = join(",",@arr);
			$worksheet1->write_row($num, 0, [decode("utf8", "time\_$num"), decode("utf8", $info)], $format{'normal'});
			
			my $row = 0;
			my $worksheet = $workbook->add_worksheet(qq{time\_$num});
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
			$num++;
		}

	}else{

		my $row = 0;
		my $worksheet = $workbook->add_worksheet(qq{lncRNA_target_result});
		open EXO, qq{$result/result/lncRNA.target.fmt.xls};
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

	print qq{lncRNA 时间序列靶基因预测分析完成!\n};

}

sub pre_check
{

	my $metadata = shift;
	my $base     = shift;

	my $organ          = qq{$metadata->{organ}};
	my $ref_lncRNA_gtf = qq{$base->{$organ}{genome_lncRNA_gtf}};
	my $ref_mRNA_gtf   = qq{$base->{$organ}{genome_mRNA_gtf}};
	my $ref_fasta      = qq{$base->{$organ}{genome_fasta}};
	# reference genome fiel is exists 
	die "reference genome fasta file is not exists!\n"      if not -e $ref_fasta;
	die "reference genome mRNA gtf file is not exists!\n"   if not -e $ref_mRNA_gtf;
	die "reference genome lncRNA gtf file is not exists!\n" if not -e $ref_lncRNA_gtf;

}

sub res_check
{
	my $out    = shift;
	my $design = shift;
	my @temp   = ();

	my $num    = 1;
	foreach my $x (sort @{$design}) {
		my $base_out = scalar @{$design} == 1 ? '' : qq{time\_$num};
		next if -e qq{$out/$base_out/target.finish};
		push @temp, $x;
		$num++;
	}
	return @temp;

}

1;
