package package::mRNA::mRNA_anova_enrichment;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'anova'}};
	my $result   = qq{$metadata->{project}/mRNA/deseq};

	my @res_groups = res_check(qq{$result/result}, \@groups, "gene");
	if (scalar @res_groups == 0) {
		print qq{mRNA anova后续富集分析已经运行完成!\n};
		#return 0;
	}  
	
	main($metadata, $base, \@res_groups, "gene");

	print qq{mRNA anova后续富集分析运行完成!\n};


}

sub main
{

	my ($metadata, $base, $res_groups, $level) = @_;

	my $organ    = qq{$metadata->{organ}};
	my @groups   = @{$metadata->{'anova'}};

	my $result          = qq{$metadata->{project}/mRNA/deseq};
	my $enrich_report   = qq{$metadata->{report}/03_mRNA_Analysis/03_Gene_Set_Enrichment_Analysis};

	my $rscript         = qq{$base->{rscript_bin}};
	my $util            = qq{$base->{util}};
	my $clusterprofiler = qq{$base->{clusterprofiler_bin}};
	my $three           = exists $base->{$organ}{"three_letter"} ? qq{$base->{$organ}{three_letter}} : "custom";

	my $pathway_db      = qq{$base->{$organ}{pathway_position_db}};
	my $go_enrich_db    = qq{$base->{$organ}{go_enrich_db}};
	my $kegg_enrich_db  = qq{$base->{$organ}{kegg_enrich_db}};
	
	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $enrich_report} if not -d $enrich_report;


	foreach my $x (@$res_groups) {

		my @ano_group = split /\,/, $x->[0];
		my $base      = join("_vs_",@ano_group);
		my $anova_tmp_dir = "$result/result/$base/$level";		
		my $anova_out_dir = "$enrich_report/$level\_level/$base\_Figures";	
		system qq{mkdir -p $anova_tmp_dir} if not -d qq{$anova_tmp_dir};
		system qq{mkdir -p $anova_out_dir} if not -d qq{$anova_out_dir};
		system qq{mkdir -p $anova_out_dir/KEGG_Pathway_Illustrations} if not -d qq{$anova_out_dir/KEGG_Pathway_Illustrations};
		system qq{mkdir -p $anova_out_dir/Custom_Figures}             if not -d qq{$anova_out_dir/Custom_Figures};
		
		my $de_gene_all = qq{perl $util/get_anova_gene_list.pl  $result/result/$base/gene_ANOVA_Test_Result.xls  $anova_tmp_dir};
	
		# ClusterProfiler的GO/KEGG分析
			# 有数据库的物种	
		my $enrich = qq{$clusterprofiler $util/clusterProfiler.R -s $three $anova_tmp_dir/de_genes.list $anova_tmp_dir};
		my $png    = qq{perl $util/get_circRNA_kegg.pl $anova_tmp_dir/kegg_enrichment.xls >$anova_tmp_dir/kegg.gene.list};
		my $url    = qq{perl $util/kegg_pathway_png.pl $pathway_db $anova_tmp_dir/kegg.gene.list >$anova_tmp_dir/kegg.urls.txt};
			# 自定义数据库的物种
		my $enrich_go    = qq{$clusterprofiler Rscript $util/clusterProfiler_go.R $anova_tmp_dir/de_genes.list $go_enrich_db $anova_tmp_dir};
		my $enrich_kegg  = qq{$clusterprofiler Rscript $util/clusterProfiler_kegg.R $anova_tmp_dir/de_genes.list $kegg_enrich_db $anova_tmp_dir};	
		my $url_custom   = qq{perl $util/kegg_pathway_png_custom.pl $kegg_enrich_db $anova_tmp_dir/kegg.gene.list >$anova_tmp_dir/kegg.urls.txt};
		
		# KEGG 通路图下载
		my $download     = qq{perl $util/download_kegg_png.pl -i $anova_tmp_dir/kegg.urls.txt -o $anova_tmp_dir/};
		my $cp_png       = qq{cp $anova_tmp_dir/png/*png $anova_out_dir/KEGG_Pathway_Illustrations};

		# KEGG 热图绘制
		my $fmt_kegg     = qq{perl $util/enrich_xls_to_heatmap.pl $anova_tmp_dir/kegg_enrichment.xls > $anova_tmp_dir/data.txt};
		my $heatmap_kegg = qq{$rscript Rscript $util/enrich_heatmap.R $anova_tmp_dir/data.txt $anova_out_dir/Custom_Figures/kegg_heatmap.pdf};

		# 图片cp到输出目录
		my $cp_go        = qq{cp $anova_tmp_dir/GO_barplot.pdf $anova_tmp_dir/go.bp.pdf $anova_tmp_dir/go.cc.pdf $anova_tmp_dir/go.mf.pdf $anova_out_dir/Custom_Figures};	
		my $cp_kegg      = qq{cp $anova_tmp_dir/kegg_dotplot.pdf $anova_out_dir/Custom_Figures};

		# finish
		my $touch        = qq{touch $anova_tmp_dir/$base.$level.enrich.finish};

		open SAVE, qq{>$result/run/$base.$level.sh} or die "Can't open $result/run/$base.$level.sh!\n";
		print SAVE qq{$de_gene_all\n};
		print SAVE qq{$mkdir\n};

		if ($three ne "custom") {
			print SAVE qq{$enrich\n};
			print SAVE qq{$png\n};
			print SAVE qq{$url\n};

		} else {
			print SAVE qq{$enrich_go\n$enrich_kegg\n};
			print SAVE qq{$png\n};
			print SAVE qq{$url_custom\n};
		}

		print SAVE qq{$download\n$cp_png\n};
		print SAVE qq{$fmt_kegg\n$heatmap_kegg\n};
		print SAVE qq{$cp_go\n$cp_kegg\n};
		print SAVE qq{$touch};
		close SAVE;

		system qq{bash $result/run/$base.$level.sh &> $result/log/$base.$level.log\n};

	}
	
	# report:GO summary xlsx
	my $excel_go      = qq{$enrich_report/$level\_level/ANOVA_GO_Enrichment_Summary.xlsx};
	my $workbook_go   = Excel::Writer::XLSX->new($excel_go);
	my %format        = package::format::run($workbook_go);

	report_xlsx($result, $workbook_go, \%format, \@groups, $level,"go");
	$workbook_go->close();

	# report:KEGG summary xlsx
	my $excel_kegg    = qq{$enrich_report/$level\_level/ANOVA_KEGG_Enrichment_Summary.xlsx};
	my $workbook_kegg = Excel::Writer::XLSX->new($excel_kegg);
	my %format        = package::format::run($workbook_kegg);

	report_xlsx($result, $workbook_kegg, \%format, \@groups, $level,"kegg");
	$workbook_kegg->close();

}

sub res_check
{
	my ($out, $groups, $level) = @_;

	my @temp   = ();
	foreach my $x (@{$groups}) {

		my @ano_group = split /\,/, $x->[0];
		my $base      = join("_vs_",@ano_group);
		next if -e qq{$out/$base/$level/$base.$level.enrich.finish};
		push @temp, $x;

	}
	return @temp;
}

sub report_xlsx
{
	my ($result, $workbook, $format, $groups, $level, $type) = @_;

	foreach my $x (@$groups) {

		my @ano_group = split /\,/, $x->[0];
		my $base      = join("_vs_",@ano_group);
		my $worksheet = $workbook->add_worksheet($base);

		open EXO, qq{$result/result/$base/$level/$type\_enrichment.xls};
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
