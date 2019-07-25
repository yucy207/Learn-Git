package package::lncRNA::lncRNA_anova;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use Encode qw/decode/;
use package::format;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'anova'}};
	my @samples  = split /,/, $metadata->{'samples'};

	my $expression   = qq{$metadata->{project}/lncRNA/gene_profile/result};
	my $result       = qq{$metadata->{project}/lncRNA/deseq};
	my $diff_report  = qq{$metadata->{report}/04_lncRNA_Analysis/02_Different_Expression_Analysis};

	my $rscript      = qq{$base->{rscript_bin}};
	my $util         = qq{$base->{util}};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $diff_report}   if not -d $diff_report;

	my @res_groups = res_check(qq{$result/result}, \@groups);

	if (scalar @res_groups == 0) {

		print qq{lncRNA anova差异分析已经运行完成!\n};
		return 0;
	}  

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);
			
	foreach my $x (@res_groups) {

		my @ano_group     = split /\,/, $x->[0];
		my $base          = join("_vs_",@ano_group);
		my $anova_tmp_dir = qq{$result/result/$base};		
		my $anova_out_dir = qq{$diff_report/$base\_ANOVA_Test_Result};	

		system qq{mkdir -p $anova_tmp_dir} if not -d qq{$anova_tmp_dir};
		system qq{mkdir -p $anova_out_dir} if not -d qq{$anova_out_dir};

		my $i = 0;
		open OUT, ">$anova_tmp_dir/sample_groups.xls" or die "Can't make sample_groups.xls!\n";
		foreach my $y (@ano_group) {
			my $ano_sample = (split /;/, $x->[1])[$i];
			$i++;
			my @arr = split/\,/,$ano_sample;
			foreach my $z(@arr){
				print OUT "$z\t$y\n";
			}	
		}

		# ANOVA (gene and transcript level)
		system qq{$rscript Rscript $util/Anova.r $expression/gene_count_matrix.csv $anova_tmp_dir/sample_groups.xls $anova_tmp_dir &> $anova_tmp_dir/gene.anova.log};
		system qq{cp $anova_tmp_dir/heatmap.pdf $anova_tmp_dir/Top50.heatmap.pdf $anova_out_dir};
		system qq{$rscript Rscript $util/trans.Anova.r $expression/transcript_count_matrix.fmt.csv $anova_tmp_dir/sample_groups.xls $anova_tmp_dir &> $anova_tmp_dir/transcript.anova.log};
		system qq{touch $anova_tmp_dir/$base.anova.finish};

		# Report
		my $excel     = qq{$anova_out_dir/gene_ANOVA_Test_Result.xlsx};	
		my $workbook  = Excel::Writer::XLSX->new($excel);		
		my %format    = package::format::run($workbook);
		my $worksheet = $workbook->add_worksheet("gene_ANOVA_Test_Result");

		report_anova("$anova_tmp_dir/gene_ANOVA_Test_Result.xls", $worksheet, \%format);
		$workbook->close();

		my $excel     = qq{$anova_out_dir/transcript_ANOVA_Test_Result.xlsx};
		my $workbook  = Excel::Writer::XLSX->new($excel);
		my %format    = package::format::run($workbook);
		my $worksheet = $workbook->add_worksheet("transcript_ANOVA_Test_Result");

		report_anova("$anova_tmp_dir/transcript_ANOVA_Test_Result.xls", $worksheet, \%format);
		$workbook->close();
		
		print qq{lncRNA anova差异分析已经运行完成!\n};

	}

}


sub res_check
{
	my $out    = shift;
	my $groups = shift;
	my @temp   = ();

	foreach my $x (@{$groups}) {

		my @ano_group = split /\s+/,$x->[0];
		my $base      = join("_vs_",@ano_group);
		next if  -e qq{$out/$base/$base.anova.finish};
		push @temp, $x;
	}
	return @temp;
	
}

sub report_anova
{
	my $txt       = shift;
	my $worksheet = shift;
	my $format	  = shift;

	my $row = 0;
	open IN, $txt or die "Can't open $txt!\n";
	while (<IN>){
		chomp;
		my @arr = split /\t/;
		if ($row == 0) {
			$worksheet->write_row( $row, 0, \@arr, $format->{'title'});
		} else {
			$worksheet->write_row( $row, 0, \@arr, $format->{'normal'});
		}
		$row++;
	}
	close IN;

}

1;
