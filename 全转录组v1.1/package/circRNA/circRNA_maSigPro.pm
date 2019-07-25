package package::circRNA::circRNA_maSigPro;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use Encode qw/decode/;
use package::format;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @design   = @{$metadata->{'time'}};
	my $result   = qq{$metadata->{project}/circRNA/maSigPro};

	# gene level
	my @res = res_check(qq{$result/result}, \@design);

	if (scalar @res == 0) {
		print qq{circRNA maSigPro 时间序列分析已经运行完成!\n};
	}

	maSigPro($metadata, $base, \@res);
	print qq{circRNA maSigPro 时间序列分析运行完成!\n};

}
	
sub maSigPro
{		
	my ($metadata, $base, $res) = @_;
	
	my @design     = @{$metadata->{'time'}};
	my $rscript    = qq{$base->{rscript_bin}};
	my $util       = qq{$base->{util}};
	my $masigpro   = qq{$base->{masigpro}};

	my $expression = qq{$metadata->{project}/circRNA/circrna_profile/result};
	my $result     = qq{$metadata->{project}/circRNA/maSigPro};
	my $report     = qq{$metadata->{report}/05_circRNA_Analysis/03_Differential_Expression_Analysis/maSigPro_Analysis};

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report}        if not -d $report;

	my $num = 1;
	foreach my $x (sort @$res) {

		my $base_out = scalar @design == 1 ? '' : qq{time\_$num};
		my $tmp_dir  = qq{$result/result/$base_out};
		my $out_dir  = qq{$report/$base_out};
		system qq{mkdir -p $tmp_dir} if not -d $tmp_dir;
		system qq{mkdir -p $out_dir} if not -d $out_dir;

		my $time      = qq{$masigpro Rscript $util/maSigPro.R $expression/circRNA.expression.xls $x $tmp_dir transcript};
		my $heatmap   = qq{$rscript Rscript $util/maSigPro_circRNA_heatmap.R $tmp_dir/all.diff.xls $out_dir};
		my $cp_design = qq{cp $x $out_dir};

		open SAVE, qq{>$result/run/maSigPro.sh};
		print SAVE qq{$time\n};
		print SAVE qq{$heatmap\n$cp_design\n};
		close SAVE;
		system qq{bash $result/run/maSigPro.sh &> $result/log/maSigPro.log};

		my $excel     = qq{$out_dir/maSigPro_analysis.xlsx};	
		my $workbook  = Excel::Writer::XLSX->new($excel);		
		my %format    = package::format::run($workbook);
		my $worksheet = $workbook->add_worksheet("maSigPro_result");
		report("$tmp_dir/all.diff.xls", $worksheet, \%format);
		$workbook->close();

		$num++;

	}

}

sub res_check
{
	my $out    = shift;
	my $design = shift;
	my @temp   = ();

	my $num    = 1;
	foreach my $x (sort @{$design}) {
		my $base_out = scalar @{$design} == 1 ? '' : qq{time\_$num};
		next if  -e qq{$out/$base_out/all.diff.xls};
		push @temp, $x;
		$num++;
	}
	return @temp;
}

sub report
{
	my $txt       = shift;
	my $worksheet = shift;
	my $format	  = shift;

	my $row = 0;
	open IN, $txt or die "Can't open $txt!\n";
	while (<IN>){
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

1;
