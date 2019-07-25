package package::lncRNA::lncRNA_enrichment;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $organ    = qq{$metadata->{organ}};

	my $target   = qq{$metadata->{project}/lncRNA/target_predict/result};
	my $result   = qq{$metadata->{project}/lncRNA/enrich_analysis};
	my $report   = qq{$metadata->{report}/04_lncRNA_Analysis/04_Target_Gene_Set_Enrichment_Analysis};

	my $util     = qq{$base->{util}};
	my $rscript  = qq{$base->{rscript_bin}};
	my $three    = qq{$base->{$organ}{three_letter}};
	my $clusterprofiler = qq{$base->{clusterprofiler_bin}};	
	my $pathway_db      = qq{$base->{$organ}{pathway_position_db}};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;

	my @res_groups = res_check(qq{$result/result},  \@groups);

	if (scalar @res_groups == 0) {
		print qq{lncRNA 靶基因富集分析已经运行完成!\n};
		#return;
	}

	pre_check($metadata, $base);

	foreach my $x (@res_groups) {

		my $control = $x->[0];
		my $case    = $x->[1];

		my $base_out   = qq{$case\_vs_$control};
		my $enrich_dir = qq{$report/$case\_vs_$control\_Figures};

		system qq{mkdir -p $enrich_dir}                  if not -d $enrich_dir;
		system qq{mkdir -p $enrich_dir/Custom_Figures}   if not -d qq{$enrich_dir/Custom_Figures};
		system qq{mkdir -p $enrich_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_dir/KEGG_Pathway_Illustrations};

		# 提取对应的信息
		my $cmd = qq{perl $util/lncRNA_target_gene_fmt.pl $target/$base_out/lncRNA.target.fmt.xls $result/result/$base_out/};

		# ClusterProfiler的GO/KEGG分析
		my $enrich = qq{$clusterprofiler $util/clusterProfiler.R -s $three $result/result/$base_out/de_genes.list $result/result/$base_out/};
		my $png    = qq{perl $util/get_circRNA_kegg.pl $result/result/$base_out/kegg_enrichment.xls > $result/result/$base_out/kegg.gene.list};
		my $url    = qq{perl $util/kegg_pathway_png.pl $pathway_db $result/result/$base_out/kegg.gene.list > $result/result/$base_out/kegg.urls.txt};
		# KEGG 通路图下载
		my $download = qq{perl $util/download_kegg_png.pl -i $result/result/$base_out/kegg.urls.txt -o $result/result/$base_out/};		
		# 图片cp到输出目录
		my $cp_png   = qq{cp $result/result/$base_out/png/*png $enrich_dir/KEGG_Pathway_Illustrations};
		my $cp       = qq{cp $result/result/$base_out/go.bp.pdf $result/result/$base_out/go.cc.pdf $result/result/$base_out/go.mf.pdf $result/result/$base_out/GO_barplot.pdf $result/result/$base_out/kegg_dotplot.pdf $enrich_dir/Custom_Figures};
		# KEGG 热图绘制
		my $fmt      = qq{perl $util/enrich_xls_to_heatmap.pl $result/result/$base_out/kegg_enrichment.xls > $result/result/$base_out/data.txt};
		my $heatmap  = qq{$rscript Rscript $util/enrich_heatmap.R $result/result/$base_out/data.txt $enrich_dir/Custom_Figures/kegg_heatmap.pdf};
		# finish
		my $touch    = qq{touch $result/result/$base_out/$case\_vs_$control.finish};

		open  SAVE, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $result/run/$case\_vs_$control.sh!\n";
		print SAVE qq{$cmd\n};
		print SAVE qq{$enrich\n};
		print SAVE qq{$fmt\n$heatmap\n};
		print SAVE qq{$png\n$url\n$download\n$cp_png\n};
		print SAVE qq{$cp\n};
		print SAVE qq{$touch\n};
		close SAVE;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log\n};
	}

	# report:GO summary xlsx
	my $excel_go    = qq{$report/GO_Enrichment_Summary.xlsx};
	my $workbook_go = Excel::Writer::XLSX->new($excel_go);
	my %format_go   = package::format::run($workbook_go);

	report_xlsx($result, $workbook_go, \%format_go , \@groups, "go");
	$workbook_go->close();

	# report:KEGG summary xlsx
	my $excel_kegg    = qq{$report/KEGG_Enrichment_Summary.xlsx};
	my $workbook_kegg = Excel::Writer::XLSX->new($excel_kegg);
	my %format_kegg   = package::format::run($workbook_kegg);

	report_xlsx($result, $workbook_kegg, \%format_kegg , \@groups, "kegg");
	$workbook_kegg->close();

	print qq{lncRNA 靶基因富集分析运行完成!\n};

}


sub pre_check
{
	my $metadata = shift;
	my $base     = shift;	

	# group config is exists 
	my @groups   = @{$metadata->{'group'}};
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
		my $control = $x->[0];
		my $case    = $x->[1];
		my $base_out         = scalar @{$groups} == 1 ? '' : qq{$case\_vs_$control};
		next if  -e qq{$out/$base_out/$case\_vs_$control.finish};
		push @temp, $x;
	}
	return @temp;
}

sub report_xlsx
{
	my ($result, $workbook, $format, $groups, $type) = @_;

	foreach my $x (@$groups) {

		my $control    = $x->[0];
		my $case       = $x->[1];
		my $base_out   = qq{$case\_vs_$control};
		my $worksheet  = $workbook->add_worksheet(qq{$case\_vs_$control});

		open EXO, qq{$result/result/$base_out/$type\_enrichment.xls};
		my $row = 0;
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			if ($row == 0) {
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

	my $id;
	if($type eq "go"){$id = "GO";}
	if($type eq "kegg"){$id = "KO";}
	my @readme = ("ID",          "$id号",
                  "Description", "$id对应的描述信息",
                  "GeneRatio",   "富集基因在该$id上的比例",
                  "BgRatio",     "背景基因在该$id上的比例",
                  "pvalue",      "P值",
                  "p.adjust",    "校正P值（BH方法校正后的P值）",
                  "qvalue",      "Q值（Q方法校正后的P值）",
                  "geneID",      "该$id上所有富集基因的基因名，用“/”分割",
                  "Count",       "该$id上所有富集基因的数量"                 
                );

	my $i = 0;
	foreach my $row(1..9){

		$worksheet->write_row($row, 0, [decode("utf8", $readme[$i]), decode("utf8", $readme[$i + 1])], $format->{'normal'});
		$i = $i + 2;

	}

}

1;
