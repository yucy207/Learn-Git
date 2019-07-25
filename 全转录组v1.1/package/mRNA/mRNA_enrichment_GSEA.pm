package package::mRNA::mRNA_enrichment_GSEA;
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
	my $three    = exists $base->{$organ}{"three_letter"} ? qq{$base->{$organ}{three_letter}} : "custom";

	my $rscript  = qq{$base->{rscript_bin}};
	my $util     = qq{$base->{util}};
	my $clusterprofiler = qq{$base->{clusterprofiler_bin}};
	my $pathway_db      = qq{$base->{$organ}{pathway_position_db}};

	if (not exists $metadata->{'group'}) {

		print "[Warning] : You must set group in the config.txt!\n";
		exit;
	}
	
	my $deseq           = qq{$metadata->{project}/mRNA/deseq};
	my $result          = qq{$metadata->{project}/mRNA/GSEA};
	my $enrich_report   = qq{$metadata->{report}/03_mRNA_Analysis/03_Gene_Set_Enrichment_Analysis/GSEA};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $enrich_report} if not -d $enrich_report;

	my @res_groups = res_check(qq{$result/result}, \@groups);
	if (scalar @res_groups == 0) {

		print qq{mRNA GSEA富集分析已经运行完成!\n};
		#return 0;
	} 

	pre_check($metadata, $base);

	foreach my $x (@res_groups) {

		my $control = $x->[0];
		my $case    = $x->[1];
		my @control_samples = split /,/, (split /;/, $x->[2])[0];
		my @case_samples    = split /,/, (split /;/, $x->[2])[1];
		my $heatmap_case    = join ",",  @case_samples;
		my $heatmap_control = join ",",  @control_samples;

		my $base_out        = qq{$case\_vs_$control};
		my $deseq_out_dir   = qq{$deseq/result/$base_out/gene};
		my $GSEA_out_dir    = qq{$result/result/$base_out};
		my $enrich_dir      = qq{$enrich_report/$case\_vs_$control\_Figures};

		system qq{mkdir -p $GSEA_out_dir}     if not -d qq{$GSEA_out_dir};
		system qq{mkdir -p $enrich_dir}       if not -d qq{$enrich_dir};
		system qq{mkdir -p $enrich_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_dir/KEGG_Pathway_Illustrations};
		system qq{mkdir -p $enrich_dir/Custom_Figures}             if not -d qq{$enrich_dir/Custom_Figures};
		
		my $mkdir  = qq{mkdir -p $enrich_dir};
		# ClusterProfiler的GO/KEGG分析	
		my $enrich = qq{$clusterprofiler $util/clusterProfiler_GSEA.R -s $three $deseq_out_dir/diff.xls $GSEA_out_dir};
		my $png    = qq{perl $util/get_GSEAkegg_png.pl $deseq_out_dir/diff.xls $GSEA_out_dir/kegg_enrichment.xls >$GSEA_out_dir/kegg.gene.list};
		my $url    = qq{perl $util/kegg_pathway_png.pl $pathway_db $GSEA_out_dir/kegg.gene.list >$GSEA_out_dir/kegg.urls.txt};
		# KEGG 通路图下载
		my $download     = qq{perl $util/download_kegg_png.pl -i $GSEA_out_dir/kegg.urls.txt -o $GSEA_out_dir/};
		my $cp_png       = qq{cp $GSEA_out_dir/png/*png $enrich_dir/KEGG_Pathway_Illustrations};
		# KEGG 热图绘制
		my $fmt_kegg     = qq{perl $util/GSEAenrich_xls_to_heatmap.pl $GSEA_out_dir/kegg_enrichment.xls > $GSEA_out_dir/data.txt};
		my $heatmap_kegg = qq{$rscript Rscript $util/enrich_heatmap.R $GSEA_out_dir/data.txt $enrich_dir/Custom_Figures/kegg_heatmap.pdf};
		# 图片cp到输出目录
		my $cp_go        = qq{cp $GSEA_out_dir/go.bp.pdf $GSEA_out_dir/go.cc.pdf $GSEA_out_dir/go.mf.pdf $enrich_dir/Custom_Figures};	
		my $cp_kegg      = qq{cp $GSEA_out_dir/kegg_dotplot.pdf $enrich_dir/Custom_Figures};
		my $cp_gsea      = qq{cp -r $GSEA_out_dir/GSEA_Figures $enrich_dir/};
		# finish
		my $touch        = qq{touch $GSEA_out_dir/$case\_vs_$control.enrich.finish};

		open SAVE, qq{>$result/run/$case\_vs_$control.sh} or die "Can't open $result/run/$case\_vs_$control.sh!\n";
		print SAVE qq{$mkdir\n};
		print SAVE qq{$enrich\n};
		print SAVE qq{$png\n};
		print SAVE qq{$url\n};
		print SAVE qq{$download\n$cp_png\n};
		print SAVE qq{$fmt_kegg\n$heatmap_kegg\n};
		print SAVE qq{$cp_kegg\n$cp_go\n$cp_gsea\n};
		print SAVE qq{$touch};
		close SAVE;

		system qq{bash $result/run/$case\_vs_$control.sh &> $result/log/$case\_vs_$control.log\n};

	}
	
	# report:GO summary xlsx
	my $excel_go      = qq{$enrich_report/GO_Enrichment_Summary.xlsx};
	my $workbook_go   = Excel::Writer::XLSX->new($excel_go);
	my %format        = package::format::run($workbook_go);
	
	report_xlsx($result, $workbook_go, \%format, \@groups, "go");
	$workbook_go->close();
		
	# report:KEGG summary xlsx
	my $excel_kegg    = qq{$enrich_report/KEGG_Enrichment_Summary.xlsx};
	my $workbook_kegg = Excel::Writer::XLSX->new($excel_kegg);
	my %format        = package::format::run($workbook_kegg);

	report_xlsx($result, $workbook_kegg, \%format, \@groups, "kegg");
	$workbook_kegg->close();

	print qq{mRNA GSEA富集分析运行完成!\n};

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	
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
		next if  -e qq{$out/$base_out/$case\_vs_$control.enrich.finish};	
		push @temp, $x;

	}
	return @temp;
}

sub report_xlsx
{
	my ($result, $workbook, $format, $groups, $type) = @_;

	foreach my $x (@$groups) {

		my $control   = $x->[0];
		my $case      = $x->[1];
		my $base_out  = qq{$case\_vs_$control};
		my $worksheet = $workbook->add_worksheet(qq{$case\_vs_$control});

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
	my @readme = ("ID",              "$id号",
                  "Description",     "$id对应的描述信息",
                  "setSize",         "每次置换检验的数据集大小",
                  "enrichmentScore", "富集分数",
                  "NES",             "归一化之后的富集分数",
                  "pvalue",          "P值",
                  "p.adjust",        "校正P值（BH方法校正后的P值）",
                  "qvalue",          "Q值（Q方法校正后的P值）",
                  "rank",            "核心基因的最大排序值",
                  "leading_edge",    "tags表示核心基因占该基因集基因总数的比例，list表示核心基因占所有基因总数的比例，singal通过公式tag*(1-list)*(N/(N-Nh))计算得到,N代表该基因集下的基因总数，Nh代表核心基因数",                
                  "core_enrichment", "该GO上富集的核心基因的基因名，用“/”分割"
                );

	my $i = 0;
	foreach my $row(1..11){
		$worksheet->write_row($row, 0, [decode("utf8", $readme[$i]), decode("utf8", $readme[$i + 1])], $format->{'normal'});
		$i = $i + 2;
	}
	if($type eq "go"){
		$worksheet->write_row(12, 0, [decode("utf8", "Type"), decode("utf8", "BP/CC/MF")], $format->{'normal'});
	}

}

1;
