package package::circRNA::circRNA_target;
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
	
	my $util     = qq{$base->{util}};
	my $miranda  = qq{$base->{miranda_bin}};
	my $miRNA    = qq{$base->{$organ}{miRNA_fasta}};

	my $deseq    = qq{$metadata->{project}/circRNA/different_analysis/result};
	my $circrna  = qq{$metadata->{project}/circRNA/circRNA_predict/result};
	my $result   = qq{$metadata->{project}/circRNA/circrna_target};
	my $report   = qq{$metadata->{report}/05_circRNA_Analysis/05_Target_Prediction};

	system qq{mkdir -p $report}        if not -d $report;
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};

	my @res_groups = res_check(qq{$result/result},  \@groups);

	if (scalar @res_groups == 0) {
		print qq{circRNA 靶标预测分析已经完成!\n};
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
		
	    system qq{mkdir -p $result/result/$base_out} if not -d qq{$result/result/$base_out};	
		# different transcript
		my $fasta  = qq{perl $util/extract_diff_circRNA_seq.pl -input $deseq/$base_out/gene/circRNA.different.expression.xls -fasta $circrna/circrnas_splice_transcript.fasta -output $result/result/$base_out/circRNA.fmt.fasta};
		my $cmd    = qq{$miranda miranda $miRNA $result/result/$base_out/circRNA.fmt.fasta -out $result/result/$base_out/miranda_targets_out};
		my $res    = qq{perl $util/fmt_miranda_out.pl $result/result/$base_out/miranda_targets_out > $result/result/$base_out/miranda.result.xls};    
		my $finish = qq{touch $result/result/$base_out/$case\_vs_$control.finish};

		open SAVE, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $result/run/$case\_vs_$control.sh\n";
		print SAVE qq{$fasta\n$cmd\n};
		print SAVE qq{$res\n$finish\n};
		close SAVE;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log};

		$pm->finish;
	}

	$pm->wait_all_children;

	################################### report ###################################
	my $excel     = qq{$report/Differential_Expression_circRNA_Targets.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);

	foreach my $x (@groups) {

		my $control   = $x->[0];
		my $case      = $x->[1];
		my $base_out  = qq{$case\_vs_$control};
		my $worksheet = $workbook->add_worksheet(qq{$case\_vs_$control});

		my %miranda_hash = miranda("$result/result/$base_out/miranda_targets_out");

		my $row = 0;
		open EXO, qq{$result/result/$base_out/miranda.result.xls};
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			if ($row == 0) {
				push @arr,"Result";
				$worksheet->write_row( $row, 0, \@arr, $format{'title'});	
			} else {
				if(exists $miranda_hash{$arr[0]}{$arr[1]}){
					push @arr,$miranda_hash{$arr[0]}{$arr[1]};
					$worksheet->write_row( $row, 0, \@arr);
				}			
			}
			$row++;
		}
		close EXO;
	}

	my $worksheet1 = $workbook->add_worksheet("README");
	my $A1 = decode("utf8", "标题"); my $A2 = decode("utf8", "说明");
	my $B1 = decode("utf8", "Seq1"); my $B2 = decode("utf8", "查询的序列1");
	my $C1 = decode("utf8", "Seq2"); my $C2 = decode("utf8", "查询的序列2");
	my $D1 = decode("utf8", "Tot Score"); my $D2 = decode("utf8", "两条序列比对情况总打分");
	my $E1 = decode("utf8", "Tot Energy"); my $E2 = decode("utf8", "两条序列结合结构的总的自由能");
    my $F1 = decode("utf8", "Max Score"); my $F2 = decode("utf8", "两条序列比对情况的最高打分");
	my $G1 = decode("utf8", "Max Energy"); my $G2 = decode("utf8", "两条序列结合结构的最大的自由能");
	my $H1 = decode("utf8", "Len1"); my $H2 = decode("utf8", "序列1的长度");
	my $I1 = decode("utf8", "Len2"); my $I2 = decode("utf8", "序列2的长度");
	my $J1 = decode("utf8", "Positions"); my $J2 = decode("utf8", "两条序列结合的位置，如果有多个位置，用|分隔");
	my $K1 = decode("utf8", "Result"); my $K2 = decode("utf8", "具体比对结果");

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

	$workbook->close();
	print qq{circRNA 靶标预测分析已经完成!\n};

}


sub pre_check
{

	my $metadata = shift;
	my $base     = shift;
	my $organ    = qq{$metadata->{organ}};
	my $miRNA    = qq{$base->{$organ}{miRNA_fasta}};
	die "genome miRNA fasta is not exists!\n" if not -e $miRNA; 

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

sub miranda
{
	my $in = shift;
	$/ = "Complete";
	my %hash = ();
	open(IN,$in) or die "cannot open $in!\n";
	while(<IN>){

		next if $_ =~/Scan Complete/;
		next if $_ =~/No Hits Found above Threshold/;
		next if $_ !~/\S+/;
		my ($miRNA)= $_ =~/Performing Scan:\s+(.*)\s+vs/;
		my ($circ) = $_ =~/vs\s+(\S+)/;	
		my @arr    = split/\n/,$_;
		my @query  = grep/Query:/, @arr;
		my ($index)= grep {$arr[$_] eq $query[0]} 0..$#arr;
		my $txt    = $arr[$index+1];
		my @ref    = grep/Ref:/, @arr;
		my $info   = "$query[0]\n$txt\n$ref[0]";
		$hash{$circ}{$miRNA} = $info;

	}
	close IN;
	$/ = "\n";
	return %hash;
}

1;
