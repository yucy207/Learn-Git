package package::mRNA::mRNA_overlap_enrichment;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @groups   = @{$metadata->{'overlap'}};
	my $result   = qq{$metadata->{project}/mRNA/deseq};

	############################### gene level ##############################
	my @res_groups = res_check(qq{$result/result}, \@groups, "gene");
	if (scalar @res_groups == 0) {
		print qq{mRNA 基因水平overlap富集分析已经运行完成!\n};
		#return 0;
	}  
	
	main($metadata, $base, \@res_groups, "gene");

	print qq{mRNA 基因水平overlap富集分析运行完成!\n};

	############################# transcript level ##########################
	my @res_groups = res_check(qq{$result/result}, \@groups,"transcript");
	if (scalar @res_groups == 0) {
		print qq{mRNA 转录本水平overlap富集分析已经运行完成!\n};
		#return 0;
	}  

	main($metadata, $base, \@res_groups,"transcript");

	print qq{mRNA 转录本水平overlap富集分析运行完成!\n};

}

sub main
{

	my ($metadata, $base, $res_groups, $level) = @_;

	my $organ    = qq{$metadata->{organ}};
	my @groups   = @{$metadata->{'overlap'}};

	my $result          = qq{$metadata->{project}/mRNA/deseq};
	my $enrich_report   = qq{$metadata->{report}/03_mRNA_Analysis/03_Gene_Set_Enrichment_Analysis};

	my $rscript         = qq{$base->{rscript_bin}};
	my $util            = qq{$base->{util}};
	my $venndiagram     = qq{$base->{venndiagram_bin}};
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

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @line  = ();
		my @names = ();

		foreach my $a (0..$#diffs) {

			my ($control, $case) = split /,/, $diffs[$a];
			my $type = $genes[$a];
			my %type2list = ('all' => 'de_genes.list',  'up' => 'up.list', 'down' => 'down.list' );
			my $file = qq{$result/result/$case\_vs_$control/$level/$type2list{$type}};
			push @names, qq{$case\_vs_$control\_$type};
			push @line , $file;
				
		}

		my $list = join " ", @line;
		my $name_list   = join "_with_", @names;
		my $de_temp_dir = qq{$result/result/$name_list/$level};
		my $enrich_dir  = qq{$enrich_report/Venn/$level\_level/$name_list\_Figures};
		system qq{mkdir -p $de_temp_dir} if not -d $de_temp_dir;
		system qq{mkdir -p $enrich_dir}  if not -d $enrich_dir;
		system qq{mkdir -p $enrich_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_dir/KEGG_Pathway_Illustrations};
		system qq{mkdir -p $enrich_dir/Custom_Figures} if not -d qq{$enrich_dir/Custom_Figures};

		# ClusterProfiler的GO/KEGG分析
		my $venn = qq{$venndiagram Rscript $util/venn_plot.r $list $de_temp_dir/};
			# 有数据库的物种	
		my $enrich = qq{$clusterprofiler $util/clusterProfiler.R -s $three $de_temp_dir/common_gene.tsv $de_temp_dir};	
		my $png    = qq{perl $util/get_circRNA_kegg.pl $de_temp_dir/kegg_enrichment.xls >$de_temp_dir/kegg.gene.list};
		my $url    = qq{perl $util/kegg_pathway_png.pl $pathway_db $de_temp_dir/kegg.gene.list >$de_temp_dir/kegg.urls.txt};
			# 自定义数据库的物种
		my $enrich_go    = qq{$clusterprofiler Rscript $util/clusterProfiler_go.R $de_temp_dir/common_gene.tsv $go_enrich_db $de_temp_dir};
		my $enrich_kegg  = qq{$clusterprofiler Rscript $util/clusterProfiler_kegg.R $de_temp_dir/common_gene.tsv $kegg_enrich_db $de_temp_dir};	
		my $url_custom   = qq{perl $util/kegg_pathway_png_custom.pl $kegg_enrich_db  $de_temp_dir/kegg.gene.list >$de_temp_dir/kegg.urls.txt};
		
		# KEGG 通路图下载
		my $download     = qq{perl $util/download_kegg_png.pl -i $de_temp_dir/kegg.urls.txt -o $de_temp_dir/};
		my $cp_png       = qq{cp $de_temp_dir/png/*png $enrich_dir/KEGG_Pathway_Illustrations};

		# KEGG 热图绘制
		my $fmt_kegg     = qq{perl $util/enrich_xls_to_heatmap.pl $de_temp_dir/kegg_enrichment.xls > $de_temp_dir/data.txt};
		my $heatmap_kegg = qq{$rscript Rscript $util/enrich_heatmap.R $de_temp_dir/data.txt $enrich_dir/Custom_Figures/kegg_heatmap.pdf};

		# 图片cp到输出目录
		my $cp_go        = qq{cp $de_temp_dir/venn.pdf $de_temp_dir/GO_barplot.pdf $de_temp_dir/go.bp.pdf $de_temp_dir/go.cc.pdf $de_temp_dir/go.mf.pdf $enrich_dir/Custom_Figures};	
		my $cp_kegg      = qq{cp $de_temp_dir/kegg_dotplot.pdf $enrich_dir/Custom_Figures};

		# finish
		my $touch        = qq{touch $de_temp_dir/$name_list.$level.enrich.finish};

		open SAVE, qq{>$result/run/$name_list.$level.sh} or die "Can't open $result/run/$name_list.$level.sh!\n";
		print SAVE qq{$venn\n};
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

		system qq{bash $result/run/$name_list.$level.sh &> $result/log/$name_list.$level.log\n};

	}
	
	# report:GO summary xlsx
	my $excel_go      = qq{$enrich_report/Venn/$level\_level/GO_Enrichment_Summary.xlsx};
	my $workbook_go   = Excel::Writer::XLSX->new($excel_go);
	my %format        = package::format::run($workbook_go);

	report_xlsx($result, $workbook_go, \%format, \@groups, $level,"go");
	$workbook_go->close();

	# report:KEGG summary xlsx
	my $excel_kegg    = qq{$enrich_report/Venn/$level\_level/KEGG_Enrichment_Summary.xlsx};
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

		my @diffs = split /;/, $x->[0];
		my @genes = split /;/, $x->[1];
		my @names = ();
		foreach my $a (0..$#diffs) {
			my ($control, $case) = split /,/, $diffs[$a];
			push @names, qq{$case\_vs_$control\_$genes[$a]};		
		}
		my $name_list   = join "_with_", @names;
		next if -e qq{$out/$name_list/$level/$name_list.$level.enrich.finish};
		push @temp, $x;

	}
	return @temp;
}


sub report_xlsx
{
	my ($result, $workbook, $format, $groups, $level, $type) = @_;

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

		open EXO, qq{$result/result/$name_list/$level/$type\_enrichment.xls};
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
