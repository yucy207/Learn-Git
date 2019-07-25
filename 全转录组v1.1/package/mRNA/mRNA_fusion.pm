package package::mRNA::mRNA_fusion;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use Encode qw/decode/;
use package::format;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my @groups   = @{$metadata->{'group'}};
	my $organ    = qq{$metadata->{organ}};

	my $mapping  = qq{$metadata->{project}/mapping/result};
	my $result   = qq{$metadata->{project}/mRNA/fusion_gene};
	my $report   = qq{$metadata->{report}/03_mRNA_Analysis/06_Fusion_gene}; 

	my $util         = qq{$base->{util}};
	my $mono         = qq{$base->{mono_bin}};
	my $fusionmap    = qq{$base->{fusionmap_bin}};
	my $fusionmap_db = qq{$base->{$organ}{fusionmap_db}};
	my $release      = qq{$base->{$organ}{fusionmap_db_release}};
	my $database     = qq{$base->{$organ}{fusionmap_db_model}};
	
	system qq{mkdir -p $result/run}       if not -d qq{$result/run};
	system qq{mkdir -p $result/result}    if not -d qq{$result/result};
	system qq{mkdir -p $result/log}       if not -d qq{$result/log};
	system qq{mkdir -p $report}           if not -d qq{$report};


	my @res_samples = res_check(qq{$result/result}, \@samples);

	if (scalar @res_samples == 0) {
		print qq{mRNA  融合基因分析已经运行完成!\n};
		#return;
	}

	pre_check($metadata, $base);

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_samples) {

		my $pid    = $pm->start and next;
		my $bam    = qq{$mapping/$x/accepted_hits.bam};
		#生成config.txt
		my $config =<<EOF;
<Files>
$bam
<Options>
MonoPath=/usr/local/bin/mono
FileFormat=BAM					// possible values: SAM, BAM. Default value=SAM 
ThreadNumber=10
RnaMode=True					//Possible values: True, False. Default value=True 
MinimalHit=2					//Possible values: 1-5, Default value =2, minimal read pairs
ReportUnannotatedFusion=False	//Possible values: True, False. Default value = False
FusionReportCutoff=1			//Possible values: 1-1000. Default value=1, require unique mapping of each read
OutputFusionReads=True			//Possible values: True, False. Default value = True  

<Output>
TempPath=$result/result/$x
OutputPath=$result/result/$x
OutputName=$x
EOF
		system qq{mkdir -p $result/result/$x} if not -d qq{$result/result/$x};
		my $config_txt = qq{>$result/result/$x/$x.config};
		open CFG, $config_txt or die "Can't make config.txt\n";
		print CFG $config;
		close CFG;
		#融合基因运行
		my $cmd    = qq{$mono mono $fusionmap --semap $fusionmap_db $release $database $result/result/$x/$x.config > $result/log/$x.log};
		my $finish = qq{touch $result/result/$x/$x.finish};

		open SAVE, qq{>$result/run/$x.sh} or die "Can't open $run/$x.sh\n";
		print SAVE qq{$cmd\n};
		print SAVE qq{$finish\n};
		close SAVE;

		system qq{bash $result/run/$x.sh &> $result/log/$x.log};

		$pm->finish;
	}

	$pm->wait_all_children;

	##################################### report #####################################
	# sample fusion gene
	my $excel     = qq{$report/Sample_Fusion_Gene.xlsx};	
	my $workbook  = Excel::Writer::XLSX->new($excel);		
	my %format    = package::format::run($workbook);
	my $worksheet = $workbook->add_worksheet("Sample_Fusion_Gene");
	
	my @head      = qw/Sample FusionID	accepted_hits.UniqueCuttingPositionCount	accepted_hits.SeedCount	accepted_hits.RescuedCount	Strand	Chromosome1	Position1	Chromosome2	Position2	KnownGene1	KnownTranscript1	KnownExonNumber1	KnownTranscriptStrand1	KnownGene2	KnownTranscript2	KnownExonNumber2	KnownTranscriptStrand2	FusionJunctionSequence	FusionGene	SplicePattern	SplicePatternClass	FrameShift	FrameShiftClass	Distance	OnExonBoundary	Filter/;
	$worksheet->write_row( 0, 0, \@head, $format{'title'});

	my %line = ();
	my $row  = 1;
	foreach my $x (@samples){	
		open IN, qq{$result/result/$x/$x.FusionReport.txt} or die "Can't open $txt!\n";
		while (<IN>){
			chomp;
			next if $_ =~/FusionID/;
			my @arr = split /\t/;
			unshift(@arr, $x);
			$line{$row} = \@arr;
			$row++;
		}
		close IN;
	}

	foreach my $k(keys %line){
		$worksheet->write_row( $k, 0, $line{$k}, $format{'normal'});
	}
	
	# extract case vs control unique fusion gene
	my $excel     = qq{$report/Group_Unique_Fusion_gene.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);

	my $format    = $workbook->add_format();
	$format->set_bold();
	$format->set_align('center');
	$format->set_bg_color('yellow');
	$format->set_border();
	$format->set_border_color('black');

	my $cell_format = $workbook->add_format();
	$cell_format->set_align('center');
	$cell_format->set_border();
	$cell_format->set_border_color('black');
	
	foreach my $x (@groups) {

		my $control = $x->[0];
		my $case    = $x->[1];
		my @control_samples1 = split /,/, (split /;/, $x->[2])[0];
		my @case_samples1    = split /,/, (split /;/, $x->[2])[1];
		my @control_samples;
		my @case_samples;

		for my $i (0..$#control_samples1){
			next if not -e "$result/result/$control_samples1[$i]/$control_samples1[$i]\.FusionReport\.txt";
			push @control_samples, $control_samples1[$i];
		}
		for my $i (0..$#case_samples1){
			next if not -e "$result/result/$case_samples1[$i]/$case_samples1[$i]\.FusionReport\.txt";
			push @case_samples, $case_samples1[$i];
		}

		my %control_id = read_count(qq{$result/result}, \@control_samples);
		my %case_id    = read_count(qq{$result/result}, \@case_samples);

		my $worksheet  = $workbook->add_worksheet(qq{$case\_vs_$control});
		my @head       = qw/Sample FusionID	accepted_hits.UniqueCuttingPositionCount	accepted_hits.SeedCount	accepted_hits.RescuedCount	Strand	Chromosome1	Position1	Chromosome2	Position2	KnownGene1	KnownTranscript1	KnownExonNumber1	KnownTranscriptStrand1	KnownGene2	KnownTranscript2	KnownExonNumber2	KnownTranscriptStrand2	FusionJunctionSequence	FusionGene	SplicePattern	SplicePatternClass	FrameShift	FrameShiftClass	Distance	OnExonBoundary	Filter/;
		$worksheet->write_row( 0, 0, \@head, $format);
		my $row = 1;

		foreach my $x (keys %case_id) {

			my @control_val = ();
			foreach my $y (@control_samples) {
				my $val = exists $control_id{$x}{$y} ? $control_id{$x}{$y}->[0] : 0;
				push @control_val, $val;
			}
			my @case_val   = ();
			foreach my $z (@case_samples) {
				my $val = exists $case_id{$x}{$z} ? $case_id{$x}{$z}->[0] : 0;
				push @case_val, $val;
			}

			my $sum_control = cal_sum(@control_val);
			my $sum_case    = cal_sum(@case_val);
			next if $sum_control > 0; 
			next if $sum_case == 0;
			my $cnt_control = join  "\t", @control_val;
			my $cnt_case   = join  "\t", @case_val;

			my ($gene1, $gene2) = split /:/, $x;
			foreach my $sample (keys %{$case_id{$x}}) {
				my @line = ();
				push @line, qq{$sample};
				foreach my $y (split /\t/, qq{$case_id{$x}{$sample}->[1]})  {
					push @line, $y;
				}
				push @line, qq{} if scalar @line != 27;
				next if $line[3] < 3;
				my @pattern = ();
				push @pattern, qq{CanonicalPattern[Major]};
				push @pattern, qq{CanonicalPattern[Minor]};
				next if not  $line[21] ~~ @pattern;
				next if $line[26] eq 'InBlackList';
				my $res = \@line;
				$worksheet->write_row( $row, 0, $res, $cell_format);
				$row++;
			}
		
		}

	}

	my $worksheet  = $workbook->add_worksheet("README");
	my $A1 = decode("utf8", "标题");     my $A2 = decode("utf8", "说明");
    my $B1 = decode("utf8", "FusionID"); my $B2 = decode("utf8", "融合基因ID");
    my $C1 = decode("utf8", "accepted_hits.UniqueCuttingPositionCount"); my $C2 = decode("utf8", "断点覆盖的Unique Reads数（去掉PCR重复）");
    my $D1 = decode("utf8", "accepted_hits.SeedCount"); my $D2 = decode("utf8", "能够比对到2个基因的Reads数（可信度高）");
    my $E1 = decode("utf8", "accepted_hits.RescuedCount"); my $E2 = decode("utf8", "能够比对到2个基因的Reads数（可信度一般）");
    my $F1 = decode("utf8", "Strand"); my $F2 = decode("utf8", "Reads比对到2个基因的方向");
    my $G1 = decode("utf8", "Chromosome1"); my $G2 = decode("utf8", "比对结果1所在染色体");
    my $H1 = decode("utf8", "Position1"); my $H2 = decode("utf8", "比对结果1所在位置");
    my $I1 = decode("utf8", "Chromosome2"); my $I2 = decode("utf8", "比对结果2所在染色体");
    my $J1 = decode("utf8", "Position2"); my $J2 = decode("utf8", "比对结果2所在位置");
    my $K1 = decode("utf8", "KnownGene1"); my $K2 = decode("utf8", "基因1（已知的）");
    my $L1 = decode("utf8", "KnownTranscript1"); my $L2 = decode("utf8", "基因1的转录本");
    my $M1 = decode("utf8", "KnownExonNumber1"); my $M2 = decode("utf8", "基因1的外显子个数");
    my $N1 = decode("utf8", "KnownTranscriptStrand1"); my $N2 = decode("utf8", "基因1的方向");
    my $O1 = decode("utf8", "KnownGene2"); my $O2 = decode("utf8", "基因2（已知的）");
    my $P1 = decode("utf8", "KnownTranscript2"); my $P2 = decode("utf8", "基因2的转录本");
    my $Q1 = decode("utf8", "KnownExonNumber2"); my $Q2 = decode("utf8", "基因2的外显子个数");
    my $R1 = decode("utf8", "KnownTranscriptStrand2"); my $R2 = decode("utf8", "基因2的方向");
    my $S1 = decode("utf8", "FusionJunctionSequence"); my $S2 = decode("utf8", "融合基因序列");
    my $T1 = decode("utf8", "FusionGene"); my $T2 = decode("utf8", "融合基因");
    my $U1 = decode("utf8", "SplicePattern"); my $U2 = decode("utf8", "剪切模式");

	my $V1 = decode("utf8", "SplicePatternClass"); my $V2 = decode("utf8", "剪切模式类型");
	my $W1 = decode("utf8", "FrameShift"); my $W2 = decode("utf8", "滑移");
	my $X1 = decode("utf8", "FrameShiftClass"); my $X2 = decode("utf8", "滑移类型");
	my $Y1 = decode("utf8", "Distance"); my $Y2 = decode("utf8", "融合基因距离");
	my $Z1 = decode("utf8", "Filter"); my $Z2 = decode("utf8", "过滤条件");
	my $AB1 = decode("utf8", "Sample"); my $AB2 = decode("utf8", "检测到该融合基因的样本");

    $worksheet->write_row(0, 0, [$A1, $A2], $format);
    $worksheet->write_row(1, 0, [$AB1, $AB2], $cell_format);
    $worksheet->write_row(2, 0, [$B1, $B2], $cell_format);
    $worksheet->write_row(3, 0, [$C1, $C2], $cell_format);
    $worksheet->write_row(4, 0, [$D1, $D2], $cell_format);
    $worksheet->write_row(5, 0, [$E1, $E2], $cell_format);
    $worksheet->write_row(6, 0, [$F1, $F2], $cell_format);
    $worksheet->write_row(7, 0, [$G1, $G2], $cell_format);
    $worksheet->write_row(8, 0, [$H1, $H2], $cell_format);
    $worksheet->write_row(9, 0, [$I1, $I2], $cell_format);
    $worksheet->write_row(10, 0, [$J1, $J2], $cell_format);
    $worksheet->write_row(11, 0, [$K1, $K2], $cell_format);
    $worksheet->write_row(12, 0, [$L1, $L2], $cell_format);
    $worksheet->write_row(13, 0, [$M1, $M2], $cell_format);
    $worksheet->write_row(14, 0, [$N1, $N2], $cell_format);
    $worksheet->write_row(15, 0, [$O1, $O2], $cell_format);
    $worksheet->write_row(16, 0, [$P1, $P2], $cell_format);
    $worksheet->write_row(17, 0, [$Q1, $Q2], $cell_format);
    $worksheet->write_row(18, 0, [$R1, $R2], $cell_format);
    $worksheet->write_row(19, 0, [$S1, $S2], $cell_format);
    $worksheet->write_row(20, 0, [$T1, $T2], $cell_format);
    $worksheet->write_row(21, 0, [$U1, $U2], $cell_format);
    $worksheet->write_row(22, 0, [$V1, $V2], $cell_format);
    $worksheet->write_row(23, 0, [$W1, $W2], $cell_format);
    $worksheet->write_row(24, 0, [$X1, $X2], $cell_format);
    $worksheet->write_row(25, 0, [$Y1, $Y2], $cell_format);
    $worksheet->write_row(26, 0, [$Z1, $Z2], $cell_format);
      
	$workbook->close();
	print qq{mRNA  融合基因分析已经运行完成!\n};
}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;

	my @samples  = split /,/, $metadata->{'samples'};
	my $mapping  = qq{$metadata->{project}/mapping/result};
	# map bam is exist
	foreach my $x (@samples) {
		my $bam = qq{$mapping/$x/accepted_hits.bam};
		die "$x mapping bam file  isn't exist!\n" if not -e $bam;
	}

}

sub res_check
{
	my $output = shift;
	my $sample = shift;
	my @temp_sample = ();
	foreach my $x (@{$sample}) {
		next if -e qq{$output/$x/$x.finish};
		push @temp_sample, $x;
	}
	return @temp_sample;
}

sub read_count
{
	my $result  = shift;
	my $samples = shift;

	my %meta    = ();
	foreach my $x (@{$samples}) {
		my $file = qq{$result/$x/$x.FusionReport.txt};
		open FUSION, $file or die "Can't open $file!\n";
		while (<FUSION>) {
			s/\n$//;
			my @arr = split /\t/;
			my $key = qq{$arr[9]:$arr[13]};
			$meta{$key}{$x} = [$arr[1], $_];
		}
		close FUSION;
	}
	return %meta;	
}

sub cal_sum
{
	my $cnt = 0;
	foreach my $x (@_) {
		$cnt += $x;
	}
	return $cnt;
}

1;
