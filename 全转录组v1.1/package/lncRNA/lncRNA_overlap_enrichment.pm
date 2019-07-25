package package::lncRNA::lncRNA_overlap_enrichment;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my $organ    = qq{$metadata->{organ}};
	my @groups   = @{$metadata->{'overlap'}};

	my $target   = qq{$metadata->{project}/lncRNA/target_predict/result};
	my $result   = qq{$metadata->{project}/lncRNA/enrich_analysis};
	my $report   = qq{$metadata->{report}/04_lncRNA_Analysis/04_Target_Gene_Set_Enrichment_Analysis/Venn};

	my $util     = qq{$base->{util}};
	my $rscript  = qq{$base->{rscript_bin}};
	my $three    = qq{$base->{$organ}{three_letter}};
	my $clusterprofiler = qq{$base->{clusterprofiler_bin}};	
	my $pathway_db      = qq{$base->{$organ}{pathway_position_db}};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;

	my @res_groups = res_check(qq{$result/result}, \@groups);

	if (scalar @res_groups == 0) {

		print qq{lncRNA overlap靶基因富集分析已经运行完成!\n};
		#return 0;
	}  
	
	foreach my $x (@res_groups) {

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @names = ();

		foreach my $a (0..$#diffs) {
			my ($control, $case) = split /,/, $diffs[$a];
			push @names, qq{$case\_vs_$control\_$genes[$a]};		
		}

		my $name_list   = join "_with_", @names;
		my $temp_dir    = qq{$result/result/$name_list};
		my $enrich_dir  = qq{$report/$name_list\_Figures};

		system qq{mkdir -p $temp_dir}    if not -d $temp_dir;
		system qq{mkdir -p $enrich_dir}  if not -d $enrich_dir;
		system qq{mkdir -p $enrich_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_dir/KEGG_Pathway_Illustrations};
		system qq{mkdir -p $enrich_dir/Custom_Figures} if not -d qq{$enrich_dir/Custom_Figures};

		# 提取对应的信息
		my $cmd    = qq{perl $util/lncRNA_target_gene_fmt.pl $target/$name_list/lncRNA.target.fmt.xls $temp_dir};
		# ClusterProfiler的GO/KEGG分析
		my $enrich = qq{$clusterprofiler $util/clusterProfiler.R -s $three $temp_dir/de_genes.list $temp_dir};	
		my $png    = qq{perl $util/get_circRNA_kegg.pl $temp_dir/kegg_enrichment.xls >$temp_dir/kegg.gene.list};
		my $url    = qq{perl $util/kegg_pathway_png.pl $pathway_db $temp_dir/kegg.gene.list >$temp_dir/kegg.urls.txt};		
		# KEGG 通路图下载
		my $download     = qq{perl $util/download_kegg_png.pl -i $temp_dir/kegg.urls.txt -o $temp_dir/};
		my $cp_png       = qq{cp $temp_dir/png/*png $enrich_dir/KEGG_Pathway_Illustrations};
		# KEGG 热图绘制
		my $fmt_kegg     = qq{perl $util/enrich_xls_to_heatmap.pl $temp_dir/kegg_enrichment.xls > $temp_dir/data.txt};
		my $heatmap_kegg = qq{$rscript Rscript $util/enrich_heatmap.R $temp_dir/data.txt $enrich_dir/Custom_Figures/kegg_heatmap.pdf};
		# 图片cp到输出目录
		my $cp_go        = qq{cp $temp_dir/venn.pdf $temp_dir/GO_barplot.pdf $temp_dir/go.bp.pdf $temp_dir/go.cc.pdf $temp_dir/go.mf.pdf $enrich_dir/Custom_Figures};	
		my $cp_kegg      = qq{cp $temp_dir/kegg_dotplot.pdf $enrich_dir/Custom_Figures};
		# finish
		my $touch        = qq{touch $temp_dir/$name_list.enrich.finish};

		open SAVE, qq{>$result/run/$name_list.sh} or die "Can't open $result/run/$name_list.sh!\n";
		print SAVE qq{$cmd\n};
		print SAVE qq{$enrich\n};
		print SAVE qq{$png\n};
		print SAVE qq{$url\n};
		print SAVE qq{$download\n$cp_png\n};
		print SAVE qq{$fmt_kegg\n$heatmap_kegg\n};
		print SAVE qq{$cp_go\n$cp_kegg\n};
		print SAVE qq{$touch};
		close SAVE;

		system qq{bash $result/run/$name_list.sh &> $result/log/$name_list.log\n};

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

	print qq{lncRNA overlap靶基因富集分析运行完成!\n};


}


sub res_check
{
	my ($out, $groups) = @_;
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
		next if -e qq{$out/$name_list/$name_list.enrich.finish};
		push @temp, $x;

	}
	return @temp;
}

sub report_xlsx
{
	my ($result, $workbook, $format, $groups, $type) = @_;

	my $worksheet1 = $workbook->add_worksheet(qq{overlap_info});
	$worksheet1->write_row( 0, 0, [decode("utf8", "编号"), decode("utf8", "overlap组别信息")], $format->{'title'});
	
	my $row1 = 1;
	foreach my $x (@$groups) {

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @names = ();

		foreach my $a (0..$#diffs) {
			my ($control, $case) = split /,/, $diffs[$a];
			push @names, qq{$case\_vs_$control\_$genes[$a]};			
		}

		my $name_list = join "_with_", @names;
		$worksheet1->write_row( $row1, 0, [decode("utf8", "overlap\_$row1"), decode("utf8",  $name_list)], $format->{'normal'});
	
		my $worksheet = $workbook->add_worksheet(qq{overlap\_$row1});
		$row1++;

		open EXO, qq{$result/result/$name_list/$type\_enrichment.xls};
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
