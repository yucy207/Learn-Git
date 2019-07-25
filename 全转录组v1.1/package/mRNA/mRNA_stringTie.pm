package package::mRNA::mRNA_stringTie;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::format;

sub run
{
	my $metadata  = shift;
	my $base      = shift;
	my @samples   = split /,/, $metadata->{'samples'};
	my $organ     = qq{$metadata->{organ}};

	my $map       = qq{$metadata->{project}/mapping/result};
	my $result    = qq{$metadata->{project}/mRNA/gene_profile}; 
	my $report    = qq{$metadata->{report}/03_mRNA_Analysis/01_mRNA_Expression_Profile};

	my $util      = qq{$base->{util}};
	my $rscript   = qq{$base->{rscript_bin}};
	my $stringtie = qq{$base->{stringtie_bin}};
	my $ref_gtf   = qq{$base->{$organ}{genome_mRNA_gtf}};
	my $ballgown  = qq{$base->{ballgown_bin}};
	my $pca_plot  = qq{$base->{pca_bin}};
	my $tsne_plot = qq{$base->{Rtsne_bin}};
	my $cor_plot  = qq{$base->{cor_bin}};

	if (not exists $metadata->{'organ'}) {

		print "[Warning] : You must set <organ> argument in the config.txt!\n";
		exit;
	}
		
	system qq{mkdir -p $report}         if not -d qq{$report};
	system qq{mkdir -p $result/run}     if not -d qq{$result/run};
	system qq{mkdir -p $result/result}  if not -d qq{$result/result};
	system qq{mkdir -p $result/log}     if not -d qq{$result/log};
	system qq{mkdir -p $report/mRNA_Expression_Figures}  if not -d qq{$report/mRNA_Expression_Figures};

	my @res_samples = res_check(qq{$result/result}, \@samples);

	if (exists $base->{'force_sample'}) {
		foreach my $x (split /,/, $base->{'force_sample'}) {
			push @res_samples, $x;
		}
	}
	if (exists $base->{'force_step'}) {
		my @steps = split /,/, $base->{'force_step'};
		@res_samples = @samples if 6 ~~ @steps;
	}

	if (scalar @res_samples == 0) {
	 	print qq{mRNA 表达定量已经运行完成!\n};
	 	return 0;
	}  
	
	pre_check($metadata, $base);

	my $max_threads = 10;
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_samples) {

		system qq{mkdir -p $result/result/$x} if not -d qq{$result/result/$x};
		my $pid     = $pm->start and next;
		my $cmd     = qq{$stringtie stringtie -p 10 -G $ref_gtf -o $result/result/$x/$x.gtf -b $result/result/$x -e $map/$x/accepted_hits.bam};
		my $finish  = qq{touch $result/result/$x/$x.stringtie.finish};

		open SAVE, qq{>$result/run/$x.sh} or die "Can't open $run/$x.sh\n";
		print SAVE qq{$cmd\n$finish\n};
		close SAVE;

		system qq{bash $result/run/$x.sh &> $result/log/$x.log};
		$pm->finish;
	}

	$pm->wait_all_children;
	
	open SAVE, qq{>$result/result/gtf.list} or die "Can't open $result/result/gtf.list!\n";
	foreach my $x (@samples) {
		print SAVE qq{$x\t$result/result/$x/$x.gtf\n};
	}	
	close SAVE;

	
	system qq{$stringtie prepDE.py -i $result/result/gtf.list -g $result/result/gene_count_matrix.csv -t $result/result/transcript_count_matrix.csv};
	system qq{sed -i 's/,/\\t/g' $result/result/gene_count_matrix.csv};
	system qq{sed -i 's/\r//g' $result/result/gene_count_matrix.csv};
	system qq{sed -i 's/,/\\t/g' $result/result/transcript_count_matrix.csv};
	system qq{sed -i 's/\r//g' $result/result/transcript_count_matrix.csv};

	system qq{perl $util/trans_count_fmt.pl -g $ref_gtf -c $result/result/transcript_count_matrix.csv > $result/result/transcript_count_matrix.fmt.csv\n};
	my $files = join ",", map {qq{$result/result/$_}} @samples;
	system qq{$ballgown Rscript $util/ballgown_fpkm.R $files $result/result &> $result/log/ballgown.log\n};
	system qq{sed -i "s/FPKM.//g" $result/result/genes.fpkm.xls};
	system qq{sed -i "s/FPKM.//g" $result/result/trans.fpkm.xls};

	# plot common figures
	system qq{$rscript Rscript $util/fpkm_density.R $result/result/genes.fpkm.xls $report/mRNA_Expression_Figures/fpkm.density.pdf &> $result/log/fpkm_density.log};
	system qq{$rscript Rscript $util/Deseq2_PCA.R $result/result/gene_count_matrix.csv $result/result/ &> $result/log/deseq2.log};
	system qq{perl $util/trans.tpm.stat.pl $result/result/gtf.list $result/result/trans.tpm.xls};
	system qq{$cor_plot Rscript $util/cor.r $result/result/count.norm.xls $report/mRNA_Expression_Figures &> $result/log/cor.log};	

	# plot pca and tsne
	if (exists $metadata->{'group'}) {

		my @groups = @{$metadata->{'group'}};
		my %unique = ();
		foreach my $x (@groups) {

			my $control  = $x->[0];
			my $case     = $x->[1];
			my $control_samples = (split /;/, $x->[2])[0];
			my $case_samples    = (split /;/, $x->[2])[1];
			$unique{$control}   = $control_samples;
			$unique{$case}      = $case_samples;

		}

		open OUT, qq{>$result/result/sample.groups} or die "Can't open $result/result/sample.groups!\n";
		foreach my $x (keys %unique) {
			my @arr = split /\,/,$unique{$x};
			foreach my $y(@arr){
				print OUT qq{$y\t$x\n};
			}	
		}

		system qq{$pca_plot Rscript  $util/pca_plot.R $result/result/count.rlog.xls $result/result/sample.groups $report/mRNA_Expression_Figures &> $result/log/pca.log};
		system qq{$tsne_plot Rscript $util/tsne_plot.R $result/result/count.rlog.xls $result/result/sample.groups $report/mRNA_Expression_Figures &> $result/log/tsne.log};  

	}      

	#xlsx表格输出
	my $excel      = qq{$report/mRNA_Expression_Summary.xlsx};
	my $workbook   = Excel::Writer::XLSX->new($excel);
	my %format     = package::format::run($workbook);
	
	my $txt1       = qq{$result/result/gene_count_matrix.csv};
	my $worksheet1 = $workbook->add_worksheet("Gene_Expression");
	report_xlsx($txt1, $worksheet1, \%format);

	my $txt2       = qq{$result/result/transcript_count_matrix.fmt.csv};
	my $worksheet2 = $workbook->add_worksheet("Transcript_Expression");
	report_xlsx($txt2, $worksheet2, \%format);

	my $txt3       = qq{$result/result/genes.fpkm.xls};
	my $worksheet3 = $workbook->add_worksheet("Gene_FPKM");
	report_xlsx($txt3, $worksheet3, \%format);

	my $txt4       = qq{$result/result/trans.fpkm.xls};
	my $worksheet4 = $workbook->add_worksheet("Transcript_FPKM");
	report_xlsx($txt4, $worksheet4, \%format);

	my $txt5       = qq{$result/result/trans.tpm.xls};
	my $worksheet5 = $workbook->add_worksheet("Transcript_TPM");
	report_xlsx($txt5, $worksheet5, \%format);

	$workbook->close();
	print qq{mRNA 表达定量已经运行完成!\n};

}


sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};

	my $map      = qq{$metadata->{project}/mapping/result};
	# map bam is exist
	foreach my $x (@samples) {
		my $bam  = qq{$map/$x/accepted_hits.bam};
		die "$x align bam file  isn't exist!\n" if not -e $bam;
	}
	# reference genome gtf is exist
	my $organ    = qq{$metadata->{organ}};
	my $ref_gtf  = qq{$base->{$organ}{genome_mRNA_gtf}};
	die "reference mRNA gtf file isn't exist!\n" if not -e $ref_gtf;
                                            
}

sub res_check
{
	my $output = shift;
	my $sample = shift;

	my @temp_sample = ();
	foreach my $x (@{$sample}) {
		next if -e qq{$output/$x/$x.stringtie.finish};
		push @temp_sample, $x;
	}
	return @temp_sample;

}

sub report_xlsx
{
	my $txt       = shift;
	my $worksheet = shift;
	my $format	  = shift;

	my $row = 0;
	open TXT, $txt or die "Can't open $txt!\n";
	while (<TXT>) {
		chomp;
		my @arr = split /\t/;
		if($row == 0){
			$worksheet->write_row( 0, 0, \@arr, $format->{'title'});
		}else{
			$worksheet->write_row( $row, 0, \@arr, $format->{'normal'});
		}
		$row++;
	}
	close TXT;
}

1;


