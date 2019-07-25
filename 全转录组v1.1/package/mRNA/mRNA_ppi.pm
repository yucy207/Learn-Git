package package::mRNA::mRNA_ppi;
use Excel::Writer::XLSX;
use package::format;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;

	my @groups   = @{$metadata->{'group'}};
	my $organ    = $metadata->{'organ'};
	my $organ_id = qq{$base->{$organ}{official_name}};

	my $diff     = qq{$metadata->{project}/mRNA/deseq/result};
	my $result   = qq{$metadata->{project}/mRNA/ppi};
	my $report   = qq{$metadata->{report}/03_mRNA_Analysis/07_Ppi_Analysis};

	my $util     = qq{$base->{util}};
	my $rscript  = qq{$base->{rscript_bin}};
	my $stringdb = qq{$base->{stringdb_bin}};
	my $clusterprofiler = qq{$base->{clusterprofiler_bin}};
	my $three           = exists $base->{$organ}{"three_letter"} ?  qq{$base->{$organ}{three_letter}} : "custom";
	
	my $pathway_db      = qq{$base->{$organ}{pathway_position_db}};
	my $go_enrich_db    = qq{$base->{$organ}{go_enrich_db}};
	my $kegg_enrich_db  = qq{$base->{$organ}{kegg_enrich_db}};;

	pre_check($metadata, $base);

	my @res_groups = res_check(qq{$result/result}, \@groups);
	if (scalar @res_groups == 0) {

		print qq{mRNA ppi分析&富集分析已经运行完成!\n};

	}  
	
	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;

	foreach my $x (@res_groups) {

		my $control     = $x->[0];
		my $case        = $x->[1];
		my $base_out    = qq{$case\_vs_$control};

		my $ppi_tmp_dir = qq{$result/result/$base_out};
		my $ppi_out_dir = qq{$report/$base_out};
		my $bash_dir    = qq{$result/run/$base_out};
		my $log_dir     = qq{$result/log/$base_out};

		# 差异基因小于10的组别跳过该分析
		my $gene_num    = `grep -v "Not DEG" $diff/$base_out/gene/diff.xls | grep -v "pvalue" | cut -f1 | sort | uniq | wc -l`;	
		if ($gene_num < 10) {
			print "分组$base_out的差异基因不足10个，因此跳过ppi分析!\n";
			next;
		}

		system qq{mkdir -p $ppi_tmp_dir} if not -d qq{$ppi_tmp_dir};
		system qq{mkdir -p $ppi_out_dir} if not -d qq{$ppi_out_dir}; 
		system qq{mkdir -p $bash_dir}    if not -d qq{$bash_dir};
		system qq{mkdir -p $log_dir}     if not -d qq{$log_dir};

		#蛋白互作分析(ppi)
		my $ppi  = qq{$stringdb Rscript $util/stringdb.r $diff/$base_out/gene/diff.xls $organ_id $ppi_tmp_dir $ppi_out_dir};
		my $cp   = qq{cp $ppi_tmp_dir/interactions.xls $ppi_out_dir};

		open SAVE, qq{>$bash_dir/$case\_vs_$control.ppi.sh} or die "Can't open $bash_dir/$case\_vs_$control.ppi.sh!\n";
		print SAVE qq{$ppi\n};
		print SAVE qq{$cp\n};

		close SAVE;

		system qq{bash $bash_dir/$case\_vs_$control.ppi.sh &> $log_dir/$case\_vs_$control.ppi.log\n};

		#富集分析
		my %cluster = ();
		open IN, qq{$ppi_out_dir/cluster.xls} or die "Can't open $ppi_out_dir/cluster.xls!\n";
		while(<IN>){
			next if $_ =~/^cluster/;
			$_ =~s/[\r\n]//g;
			my ($id,$gene) = split/\t/,$_;
			push @{$cluster{$id}},$gene;
		}
		close IN;

		my $excel_go      = qq{$ppi_out_dir/GO_Enrichment_Summary.xlsx};
		my $workbook_go   = Excel::Writer::XLSX->new($excel_go);
		my %format_go     = package::format::run($workbook_go);

		my $excel_kegg    = qq{$ppi_out_dir/KEGG_Enrichment_Summary.xlsx};
		my $workbook_kegg = Excel::Writer::XLSX->new($excel_kegg);
		my %format_kegg   = package::format::run($workbook_kegg);

		foreach my $k(sort keys %cluster){

			my $num = @{$cluster{$k}};
			next if $num < 10;
			my $enrich_tmp_dir = qq{$ppi_tmp_dir/cluster\_$k};
			my $enrich_out_dir = qq{$ppi_out_dir/cluster\_$k};
			system qq{mkdir -p $enrich_tmp_dir} if not -d qq{$enrich_tmp_dir};
			system qq{mkdir -p $enrich_out_dir} if not -d qq{$enrich_out_dir};
			system qq{mkdir -p $enrich_out_dir/KEGG_Pathway_Illustrations} if not -d qq{$enrich_out_dir/KEGG_Pathway_Illustrations};
			system qq{mkdir -p $enrich_out_dir/Custom_Figures}             if not -d qq{$enrich_out_dir/Custom_Figures};

			next if -e qq{$enrich_tmp_dir/cluster\_$k.enrich.finish};

			open DE, qq{>$enrich_tmp_dir/de_genes.list} or die "Can't make $enrich_tmp_dir/de_genes.list!\n";
			my $list = join("\n",@{$cluster{$k}});
			print DE qq{$list};
			close DE;

			# 常规物种ClusterProfiler的GO/KEGG分析
			my $enrich   = qq{$clusterprofiler $util/clusterProfiler.R -s $three $enrich_tmp_dir/de_genes.list $enrich_tmp_dir};
			my $png      = qq{perl $util/get_kegg_png.pl $diff/$base_out/gene/diff.xls $enrich_tmp_dir/kegg_enrichment.xls > $enrich_tmp_dir/kegg.gene.list};
			my $url      = qq{perl $util/kegg_pathway_png.pl $pathway_db $enrich_tmp_dir/kegg.gene.list > $enrich_tmp_dir/kegg.urls.txt};	
			my $download = qq{perl $util/download_kegg_png.pl -i $enrich_tmp_dir/kegg.urls.txt -o $enrich_tmp_dir/};
		    my $cp_png   = qq{cp $enrich_tmp_dir/png/*png $enrich_out_dir/KEGG_Pathway_Illustrations};
		    # KEGG 热图绘制
			my $fmt_kegg      = qq{perl $util/enrich_xls_to_heatmap.pl $enrich_tmp_dir/kegg_enrichment.xls > $enrich_tmp_dir/data.txt};
			my $heatmap_kegg  = qq{$rscript Rscript $util/enrich_heatmap.R $enrich_tmp_dir/data.txt $enrich_out_dir/Custom_Figures/kegg_heatmap.pdf};
			# 图片cp到输出目录
			my $cp_go    = qq{cp $enrich_tmp_dir/GO_barplot.pdf $enrich_tmp_dir/go.bp.pdf  $enrich_tmp_dir/go.cc.pdf $enrich_tmp_dir/go.mf.pdf $enrich_out_dir/Custom_Figures};
			my $cp_kegg  = qq{cp $enrich_tmp_dir/kegg_dotplot.pdf $enrich_out_dir/Custom_Figures};		
			# finish
			my $touch    = qq{touch $enrich_tmp_dir/cluster\_$k.enrich.finish};

			open SAVE, qq{>$bash_dir/cluster\_$k.enrich.sh} or die "Can't open $bash_dir/cluster\_$k.enrich.sh!\n";
			print SAVE qq{$enrich\n};
			print SAVE qq{$png\n};
			print SAVE qq{$url\n};
			print SAVE qq{$download\n$cp_png\n};
			print SAVE qq{$fmt_kegg\n$heatmap_kegg\n};
			print SAVE qq{$cp_kegg\n$cp_go\n};
			print SAVE qq{$touch};
			close SAVE;

			system qq{bash $bash_dir/cluster\_$k.enrich.sh &> $log_dir/cluster\_$k.enrich.log\n};

			# report:GO summary xlsx
			report_xlsx($enrich_tmp_dir, $workbook_go, \%format_go, $k, "go");		
			# report:KEGG summary xlsx
			report_xlsx($enrich_tmp_dir, $workbook_kegg, \%format_kegg, $k, "kegg");

		}

		# report:GO summary readme
		report_readme($workbook_go, \%format_go,"go");	
		$workbook_go->close();	
		# report:KEGG summary readme
		report_readme($workbook_kegg, \%format_kegg,"kegg");
		$workbook_kegg->close();

		system qq{touch $ppi_tmp_dir/$case\_vs_$control.finish};

	}

	print qq{mRNA ppi分析&富集分析已经运行完成!\n};

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my $organ    = $metadata->{'organ'};
	
	my $organ_id = qq{$base->{$organ}{official_name}};
	die "organ official name is not exists!\n" if not defined $organ_id;
	my @groups   = @{$metadata->{'group'}};
	die "group config info is not exists!\n" if scalar @groups == 0;

	foreach my $x (@groups) {
		my $control = $x->[0];
		my $case    = $x->[1];
		my $sample  = $x->[2];

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
		my $control   = $x->[0];
		my $case      = $x->[1];
		my $base_out  = qq{$case\_vs_$control};
		next if  -e qq{$out/$base_out/$case\_vs_$control.finish};	
		push @temp, $x;
	}
	return @temp;
}

sub report_xlsx
{
	my ($enrich_tmp_dir, $workbook, $format, $k, $type) = @_;

	my $worksheet = $workbook->add_worksheet(qq{cluster\_$k});
	open IN, qq{$enrich_tmp_dir/go_enrichment.xls};
	my $row = 0;
	while (<IN>) {
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

sub report_readme
{
	my ($workbook, $format, $type) = @_;

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
