package package::circRNA::circRNA_enrich;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'group'}};
	my $organ    = qq{$metadata->{organ}};

	my $rscript  = qq{$base->{rscript_bin}};
	my $util     = qq{$base->{util}};
	my $three    = qq{$base->{$organ}{three_letter}};

	my $clusterprofiler     =  qq{$base->{clusterprofiler_bin}};
	my $pathway_position_db = qq{$base->{$organ}{pathway_position_db}};

	my $deseq    = qq{$metadata->{project}/circRNA/different_analysis/result};
	my $result   = qq{$metadata->{project}/circRNA/enrich_analysis};
	my $report   = qq{$metadata->{report}/05_circRNA_Analysis/04_Parental_Gene_Set_Enrichment_Analysis};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;


	my @res_groups = res_check(qq{$result/result}, \@groups);

	if (scalar @res_groups == 0) {
		print qq{circRNA 来源基因富集分析已经完成!\n};
		#return;
	}

	foreach my $x (@res_groups) {

		my $control  = $x->[0];
		my $case     = $x->[1];

		my $base_out          = qq{$case\_vs_$control};
		my $de_reprot_dir     = qq{$deseq/$base_out/gene};
		my $enrich_out_dir    = qq{$result/result/$base_out};	
		my $enrich_report_dir = qq{$report/$case\_vs_$control\_Figures};

		system qq{mkdir -p $enrich_out_dir}    if not -d $enrich_out_dir;
		system qq{mkdir -p $enrich_report_dir} if not -d $enrich_report_dir;
		system qq{mkdir -p $enrich_report_dir/Custom_Figures} if not -d qq{$enrich_report_dir/Custom_Figures};
		system qq{mkdir -p $enrich_report_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_report_dir/KEGG_Pathway_Illustrations};
		
		# ClusterProfiler的GO/KEGG分析	
		my $extract  = qq{perl $util/extract_circRNA_gene.pl $de_reprot_dir/circRNA.different.expression.xls $de_reprot_dir};
		my $cmd      = qq{$clusterprofiler $util/clusterProfiler.R -s $three $de_reprot_dir/de_genes.list $enrich_out_dir};			
		my $png      = qq{perl $util/get_circRNA_kegg.pl $enrich_out_dir/kegg_enrichment.xls > $enrich_out_dir/kegg_gene.list};
		my $url      = qq{perl $util/kegg_pathway_png.pl $pathway_position_db $enrich_out_dir/kegg_gene.list > $enrich_out_dir/kegg.urls.txt};
		# KEGG 通路图下载
		my $download = qq{perl $util/download_kegg_png.pl -i $enrich_out_dir/kegg.urls.txt -o $enrich_out_dir};
		# 图片cp到输出目录
		my $cp_png   = qq{cp $enrich_out_dir/png/*png $enrich_report_dir/KEGG_Pathway_Illustrations};
		my $cp       = qq{cp $enrich_out_dir/go.bp.pdf $enrich_out_dir/go.cc.pdf $enrich_out_dir/go.mf.pdf $enrich_out_dir/GO_barplot.pdf $enrich_out_dir/kegg_dotplot.pdf $enrich_report_dir/Custom_Figures};
		# KEGG 热图绘制
		my $fmt     = qq{perl $util/enrich_xls_to_heatmap.pl $enrich_out_dir/kegg_enrichment.xls > $result/result/$base_out/data.txt};
		my $heatmap = qq{$rscript Rscript $util/enrich_heatmap.R $result/result/$base_out/data.txt $enrich_report_dir/Custom_Figures/kegg_heatmap.pdf};
		# finish
		my $finish  = qq{touch $result/result/$base_out/$case\_vs_$control.finish};

		open SAVE, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $result/run/$case\_vs_$control.sh!\n";
		print SAVE qq{$extract\n};
		print SAVE qq{$cmd\n};
		print SAVE qq{$png\n$url\n$download\n$cp_png\n$cp\n};		
		print SAVE qq{$fmt\n$heatmap\n$finish};
		close SAVE;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log\n};

	}

	# report:GO summary xlsx
	my $excel_go     = qq{$report/GO_Enrichment_Summary.xlsx};
	my $workbook_go  = Excel::Writer::XLSX->new($excel_go);
	my %format_go    = package::format::run($workbook_go);

	report_xlsx($result, $workbook_go, \%format_go, \@groups, "go");
	$workbook_go->close();

	# report:KEGG summary xlsx
	my $excel_kegg     = qq{$report/KEGG_Enrichment_Summary.xlsx};
	my $workbook_kegg  = Excel::Writer::XLSX->new($excel_kegg);
	my %format_kegg    = package::format::run($workbook_kegg);

	report_xlsx($result, $workbook_kegg, \%format_kegg, \@groups,"kegg");
	$workbook_kegg->close();

	print qq{circRNA 来源基因富集分析运行完成!\n};
	
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
