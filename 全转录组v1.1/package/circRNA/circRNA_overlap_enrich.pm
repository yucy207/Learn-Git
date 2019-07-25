package package::circRNA::circRNA_overlap_enrich;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{

	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'overlap'}};
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
		print qq{circRNA overlap基因富集分析已经完成!\n};
		#return;
	}

	foreach my $x (@res_groups) {

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @line  = ();
		my @names = ();

		foreach my $a (0..$#diffs) {

			my ($control, $case) = split /,/, $diffs[$a];		
			my $type = $genes[$a];
			my %type2list = ('all' => 'de_genes.list',  'up' => 'up.list', 'down' => 'down.list' );
			my $file = qq{$deseq/$case\_vs_$control/gene/$type2list{$type}};
			push @names, qq{$case\_vs_$control\_$type};
			push @line , $file;
				
		}

		my $list = join " ", @line;
		my $name_list         = join "_with_", @names;
		my $enrich_out_dir    = qq{$result/result/$name_list};	
		my $enrich_report_dir = qq{$report/Venn/$name_list\_Figures};
	
		system qq{mkdir -p $enrich_out_dir}    if not -d $enrich_out_dir;
		system qq{mkdir -p $enrich_report_dir} if not -d $enrich_report_dir;
		system qq{mkdir -p $enrich_report_dir/Custom_Figures} if not -d qq{$enrich_report_dir/Custom_Figures};
		system qq{mkdir -p $enrich_report_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_report_dir/KEGG_Pathway_Illustrations};

		my $venn     = qq{$venndiagram Rscript $util/venn_plot.r $list $enrich_out_dir};
		# ClusterProfiler的GO/KEGG分析	
		my $cmd      = qq{$clusterprofiler $util/clusterProfiler.R -s $three $enrich_out_dir/common_gene.tsv $enrich_out_dir};			
		my $png      = qq{perl $util/get_circRNA_kegg.pl $enrich_out_dir/kegg_enrichment.xls > $enrich_out_dir/kegg_gene.list};
		my $url      = qq{perl $util/kegg_pathway_png.pl $pathway_position_db $enrich_out_dir/kegg_gene.list > $enrich_out_dir/kegg.urls.txt};
		# KEGG 通路图下载
		my $download = qq{perl $util/download_kegg_png.pl -i $enrich_out_dir/kegg.urls.txt -o $enrich_out_dir};
		# 图片cp到输出目录
		my $cp_png   = qq{cp $enrich_out_dir/png/*png $enrich_report_dir/KEGG_Pathway_Illustrations};
		my $cp       = qq{cp $enrich_out_dir/go.bp.pdf $enrich_out_dir/go.cc.pdf $enrich_out_dir/go.mf.pdf $enrich_out_dir/GO_barplot.pdf $enrich_out_dir/kegg_dotplot.pdf $enrich_report_dir/Custom_Figures};
		# KEGG 热图绘制
		my $fmt     = qq{perl $util/enrich_xls_to_heatmap.pl $enrich_out_dir/kegg_enrichment.xls > $result/result/$name_list/data.txt};
		my $heatmap = qq{$rscript Rscript $util/enrich_heatmap.R $result/result/$name_list/data.txt $enrich_report_dir/Custom_Figures/kegg_heatmap.pdf};
		# finish
		my $finish  = qq{touch $enrich_out_dir/$name_list.finish};

		open SAVE, qq{>$result/run/$name_list.sh} or die "Can't open $result/run/$name_list.sh!\n";
		print SAVE qq{$venn\n};
		print SAVE qq{$cmd\n};
		print SAVE qq{$png\n$url\n$download\n$cp_png\n$cp\n};		
		print SAVE qq{$fmt\n$heatmap\n$finish};
		close SAVE;

		system qq{bash $result/run/$name_list.sh &> $result/log/$name_list.log\n};

	}
	
	# report:GO summary xlsx
	my $excel_go     = qq{$report/Venn/GO_Enrichment_Summary.xlsx};
	my $workbook_go  = Excel::Writer::XLSX->new($excel_go);
	my %format_go    = package::format::run($workbook_go);

	report_xlsx($result, $workbook_go, \%format_go, \@groups, "go");
	$workbook_go->close();

	# report:KEGG summary xlsx
	my $excel_kegg     = qq{$report/Venn/KEGG_Enrichment_Summary.xlsx};
	my $workbook_kegg  = Excel::Writer::XLSX->new($excel_kegg);
	my %format_kegg    = package::format::run($workbook_kegg);

	report_xlsx($result, $workbook_kegg, \%format_kegg, \@groups,"kegg");
	$workbook_kegg->close();

	print qq{circRNA overlap基因富集分析运行完成!\n};

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
		my $name_list   = join "_with_", @names;
		next if -e qq{$out/$name_list/$name_list.finish};
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
