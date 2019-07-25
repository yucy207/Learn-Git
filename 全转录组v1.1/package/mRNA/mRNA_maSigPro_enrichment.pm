package package::mRNA::mRNA_maSigPro_enrichment;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{

	my $metadata = shift;
	my $base     = shift;
	my @design   = @{$metadata->{'time'}};
	my $result   = qq{$metadata->{project}/mRNA/maSigPro};

	# gene level
	my @res = res_check(qq{$result/result/gene}, \@design);

	if (scalar @res == 0) {
		print qq{mRNA maSigPro 基因水平时间序列富集分析已经运行完成!\n};
	}

	main($metadata, $base, @res, "gene");
	print qq{mRNA maSigPro 基因水平时间序列富集分析运行完成!\n};

	# transcript level
	my @res = res_check(qq{$result/result/transcript}, \@design);

	if (scalar @res == 0) {
		print qq{mRNA maSigPro 转录本水平时间序列富集分析已经运行完成!\n};
	}

	main($metadata, $base, @res,"transcript");
	print qq{mRNA maSigPro 转录本水平时间序列富集分析运行完成!\n};

}

sub main
{

	my ($metadata, $base, $res, $level) = @_;

	my $organ           = qq{$metadata->{organ}};
	my $result          = qq{$metadata->{project}/mRNA/maSigPro};
	my $enrich_report   = qq{$metadata->{report}/03_mRNA_Analysis/03_Gene_Set_Enrichment_Analysis/maSigPro};
	my @design          = @{$metadata->{'time'}};

	my $rscript         = qq{$base->{rscript_bin}};
	my $util            = qq{$base->{util}};
	my $clusterprofiler = qq{$base->{clusterprofiler_bin}};
	my $three           = exists $base->{$organ}{"three_letter"} ? qq{$base->{$organ}{three_letter}} : "custom";
	my $pathway_db      = qq{$base->{$organ}{pathway_position_db}};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $enrich_report} if not -d $enrich_report;

	my $num = 1;
	foreach my $x (sort @$res) {

		my $base_out   = scalar @design == 1 ? '' : qq{time\_$num};
		my $tmp_dir    = qq{$result/result/$level/$base_out};
		my $enrich_dir = qq{$enrich_report/$level\_level/$base_out};
		system qq{mkdir -p $tmp_dir}    if not -d $tmp_dir;
		system qq{mkdir -p $enrich_dir} if not -d $enrich_dir;
		system qq{mkdir -p $enrich_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_dir/KEGG_Pathway_Illustrations};
		system qq{mkdir -p $enrich_dir/Custom_Figures}             if not -d qq{$enrich_dir/Custom_Figures};

		my $de_gene_all;
		if ($level eq "gene")      {$de_gene_all = qq{grep -v "Not DEG" $tmp_dir/all.diff.xls | grep -v "p-value" | cut -f1 | sort | uniq >$tmp_dir/de_genes.list};}
		if ($level eq "transcript"){$de_gene_all = qq{grep -v "Not DEG" $tmp_dir/all.diff.xls | grep -v "p-value" | cut -f2 | sort | uniq >$tmp_dir/de_genes.list};}
		# ClusterProfiler的GO/KEGG分析
			# 有数据库的物种	
		my $enrich = qq{$clusterprofiler $util/clusterProfiler.R -s $three $tmp_dir/de_genes.list $tmp_dir};
		my $png    = qq{perl $util/get_circRNA_kegg.pl $tmp_dir/kegg_enrichment.xls >$tmp_dir/kegg.gene.list};
		my $url    = qq{perl $util/kegg_pathway_png.pl $pathway_db $tmp_dir/kegg.gene.list >$tmp_dir/kegg.urls.txt};
		
		# KEGG 通路图下载
		my $download     = qq{perl $util/download_kegg_png.pl -i $tmp_dir/kegg.urls.txt -o $tmp_dir/};
		my $cp_png       = qq{cp $tmp_dir/png/*png $enrich_dir/KEGG_Pathway_Illustrations};

		# KEGG 热图绘制
		my $fmt_kegg     = qq{perl $util/enrich_xls_to_heatmap.pl $tmp_dir/kegg_enrichment.xls > $tmp_dir/data.txt};
		my $heatmap_kegg = qq{$rscript Rscript $util/enrich_heatmap.R $tmp_dir/data.txt $enrich_dir/Custom_Figures/kegg_heatmap.pdf};

		# 图片cp到输出目录
		my $cp_go        = qq{cp $tmp_dir/GO_barplot.pdf $tmp_dir/go.bp.pdf  $tmp_dir/go.cc.pdf $tmp_dir/go.mf.pdf $enrich_dir/Custom_Figures};	
		my $cp_kegg      = qq{cp $tmp_dir/kegg_dotplot.pdf $enrich_dir/Custom_Figures};

		# finish
		my $touch        = qq{touch $tmp_dir/enrich.finish};

		open SAVE, qq{>$result/run/$level$base_out.enrich.sh} or die "Can't open $result/run/$level$base_out.enrich.sh!\n";
		print SAVE qq{$de_gene_all\n};
		print SAVE qq{$enrich\n};
		print SAVE qq{$png\n};
		print SAVE qq{$url\n};
		print SAVE qq{$download\n$cp_png\n};
		print SAVE qq{$fmt_kegg\n$heatmap_kegg\n};
		print SAVE qq{$cp_go\n$cp_kegg\n};
		print SAVE qq{$touch};
		close SAVE;

		system qq{bash $result/run/$level$base_out.enrich.sh &> $result/log/$level$base_out.enrich.log\n};
		$num++;
	}

	# report:GO summary xlsx
	my $excel_go      = qq{$enrich_report/$level\_level/GO_Enrichment_Summary.xlsx};
	my $workbook_go   = Excel::Writer::XLSX->new($excel_go);
	my %format        = package::format::run($workbook_go);

	report_xlsx($result, $workbook_go, \%format,\@design, $level, "go");
	$workbook_go->close();

	# report:KEGG summary xlsx
	my $excel_kegg    = qq{$enrich_report/$level\_level/KEGG_Enrichment_Summary.xlsx};
	my $workbook_kegg = Excel::Writer::XLSX->new($excel_kegg);
	my %format        = package::format::run($workbook_kegg);

	report_xlsx($result, $workbook_kegg, \%format, \@design, $level,"kegg");
	$workbook_kegg->close();
}

sub res_check
{
	my $out    = shift;
	my $design = shift;
	my @temp   = ();

	my $num    = 1;
	foreach my $x (sort @{$design}) {
		my $base_out = scalar @{$design} == 1 ? '' : qq{time\_$num};
		next if  -e qq{$out/$base_out/enrich.finish};
		push @temp, $x;
		$num++;
	}
	return @temp;
}

sub report_xlsx
{
	my ($result, $workbook, $format, $design, $level, $type) = @_;

	if(scalar @$design != 1){

		my $worksheet1 = $workbook->add_worksheet(qq{MaSigPro_info});
		$worksheet1->write_row( 0, 0, [decode("utf8", "编号"), decode("utf8", "时间信息")], $format{'title'});
		
		my $num = 1;
		foreach my $x (sort @$design) {

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
			my $worksheet  = $workbook->add_worksheet(qq{time\_$num});
			open EXO, qq{$result/result/$level/$base_out/$type\_enrichment.xls};	
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
			$num++;
		}
		
	}else{

		my $worksheet = $workbook->add_worksheet(qq{$type\_enrichment});
		open EXO, qq{$result/result/$level/$type\_enrichment.xls};
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
