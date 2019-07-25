use Excel::Writer::XLSX;
use Encode qw/decode/;
my ($result_dir, $out_dir) = @ARGV;

my %A3SS_only   = read_file(qq{$result_dir/A3SS.MATS.JunctionCountOnly.xls});
my %A3SS_target = read_file(qq{$result_dir/A3SS.MATS.ReadsOnTargetAndJunctionCounts.xls});

my %A5SS_only   = read_file(qq{$result_dir/A5SS.MATS.JunctionCountOnly.xls});
my %A5SS_target = read_file(qq{$result_dir/A5SS.MATS.ReadsOnTargetAndJunctionCounts.xls});

my %MXE_only    = read_file(qq{$result_dir/MXE.MATS.JunctionCountOnly.xls});
my %MXE_target  = read_file(qq{$result_dir/MXE.MATS.ReadsOnTargetAndJunctionCounts.xls});

my %RI_only     = read_file(qq{$result_dir/RI.MATS.JunctionCountOnly.xls});
my %RI_target   = read_file(qq{$result_dir/RI.MATS.ReadsOnTargetAndJunctionCounts.xls});

my %SE_only     = read_file(qq{$result_dir/SE.MATS.JunctionCountOnly.xls});
my %SE_target   = read_file(qq{$result_dir/SE.MATS.ReadsOnTargetAndJunctionCounts.xls});

my %A3SS_both   = cal_overlap(\%A3SS_only, \%A3SS_target);
my %A5SS_both   = cal_overlap(\%A5SS_only, \%A5SS_target);
my %MXE_both    = cal_overlap(\%MXE_only,  \%MXE_target);
my %RI_both     = cal_overlap(\%RI_only,   \%RI_target);
my %SE_both     = cal_overlap(\%SE_only,   \%SE_target);

my @A3SS_tiltle = qw/ID	GeneID	geneSymbol	chr	strand	isoform_exon	longExonStart_0base	longExonEnd	isoform_up	shortES	shortEE	isoform_down	flankingES	flankingEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @A5SS_tiltle = qw/ID	GeneID	geneSymbol	chr	strand	isoform_exon	longExonStart_0base	longExonEnd	isoform_up	shortES	shortEE	isoform_down	flankingES	flankingEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @MXE_tiltle  = qw/ID	GeneID	geneSymbol	chr	strand	isoform_exon1	1stExonStart_0base	1stExonEnd	isoform_exon2	2ndExonStart_0base	2ndExonEnd	isoform_up	upstreamES	upstreamEE	isoform_down	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @RI_tiltle   = qw/ID	GeneID	geneSymbol	chr	strand	isoform_exon	riExonStart_0base	riExonEnd	isoform_up	upstreamES	upstreamEE	isoform_down	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @SE_tiltle   = qw/ID	GeneID	geneSymbol	chr	strand	isoform_exon	exonStart_0base	exonEnd	isoform_up	upstreamES	upstreamEE	isoform_down	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;

my $excel     = qq{$out_dir/alternative_splicing.xlsx};
my $workbook  = Excel::Writer::XLSX->new($excel);
my $format = $workbook->add_format();
$format->set_bold();
$format->set_align('center');
$format->set_bg_color('yellow');
$format->set_border();
$format->set_border_color('black');

my $cell_format = $workbook->add_format();
$cell_format->set_align('center');
$cell_format->set_border();
$cell_format->set_border_color('black');


my $worksheet1  = $workbook->add_worksheet("A3SS");
my $worksheet2  = $workbook->add_worksheet("A5SS");
my $worksheet3  = $workbook->add_worksheet("MXE");
my $worksheet4  = $workbook->add_worksheet("RI");

my $worksheet5  = $workbook->add_worksheet("SE");
my $worksheet6  = $workbook->add_worksheet("README");

my $row = 1;
$worksheet1->write_row( 0, 0, \@A3SS_tiltle, $format);

open SAVE, qq{>$result_dir/A3SS.overlap.txt};
my $tilte = join "\t",  @A3SS_tiltle;
print SAVE qq{$tilte\n};

foreach my $x (keys %A3SS_both) {
	my $line1 = $A3SS_both{$x}->[0];
	my $line2 = $A3SS_both{$x}->[1];

	print SAVE qq{$line1\n};
	my @line1 = split /\t/, $line1;
	my @line2 = split /\t/, $line2;

	$worksheet1->write_row( $row, 0, \@line1, $cell_format);
	$row++;

}
close SAVE;



my $row = 1;
$worksheet2->write_row( 0, 0, \@A5SS_tiltle, $format);
open SAVE, qq{>$result_dir/A5SS.overlap.txt};
my $tilte = join "\t",  @A5SS_tiltle;
print SAVE qq{$tilte\n};
foreach my $x (keys %A5SS_both) {
	my $line1 = $A5SS_both{$x}->[0];
	my $line2 = $A5SS_both{$x}->[1];
	print SAVE qq{$line1\n};
	my @line1 = split /\t/, $line1;
	my @line2 = split /\t/, $line2;

	$worksheet2->write_row( $row, 0, \@line1, $cell_format);
	$row++;

}
close SAVE;


my $row = 1;
$worksheet3->write_row( 0, 0, \@MXE_tiltle, $format);
open SAVE, qq{>$result_dir/MXE.overlap.txt};
my $tilte = join "\t",  @MXE_tiltle;
print SAVE qq{$tilte\n};
foreach my $x (keys %MXE_both) {
	my $line1 = $MXE_both{$x}->[0];
	my $line2 = $MXE_both{$x}->[1];
	print SAVE qq{$line1\n};
	my @line1 = split /\t/, $line1;
	my @line2 = split /\t/, $line2;


	$worksheet3->write_row( $row, 0, \@line1, $cell_format);
	$row++;

}

close SAVE;
my $row = 1;
$worksheet4->write_row( 0, 0, \@RI_tiltle, $format);
open SAVE, qq{>$result_dir/RI.overlap.txt};
my $tilte = join "\t",  @RI_tiltle;
print SAVE qq{$tilte\n};
foreach my $x (keys %RI_both) {
	my $line1 = $RI_both{$x}->[0];
	my $line2 = $RI_both{$x}->[1];
	print SAVE qq{$line1\n};
	my @line1 = split /\t/, $line1;
	my @line2 = split /\t/, $line2;

	$worksheet4->write_row( $row, 0, \@line1, $cell_format);
	$row++;

}
close SAVE;


my $row = 1;
$worksheet5->write_row( 0, 0, \@SE_tiltle, $format);
open SAVE, qq{>$result_dir/SE.overlap.txt};
my $tilte = join "\t",  @SE_tiltle;
print SAVE qq{$tilte\n};
foreach my $x (keys %SE_both) {
	my $line1 = $SE_both{$x}->[0];
	my $line2 = $SE_both{$x}->[1];
	print SAVE qq{$line1\n};
	my @line1 = split /\t/, $line1;
	my @line2 = split /\t/, $line2;

	$worksheet5->write_row( $row, 0, \@line1, $cell_format);
	$row++;

}
close SAVE;

my $A1 = decode("utf8", "标题"); my $A2 = decode("utf8", "说明");
my $B1 = decode("utf8", "ID"); my $B2 = decode("utf8", "可变剪切事件的编号");
my $C1 = decode("utf8", "GeneID"); my $C2 = decode("utf8", "基因ID");
my $D1 = decode("utf8", "geneSymbol"); my $D2 = decode("utf8", "基因名称");
my $E1 = decode("utf8", "chr"); my $E2 = decode("utf8", "染色体");
my $F1 = decode("utf8", "strand"); my $F2 = decode("utf8", "基因所在链的方向");
my $G1 = decode("utf8", "ID"); my $G2 = decode("utf8", "可变剪切事件的编号");
my $H1 = decode("utf8", "IJC_SAMPLE_1"); my $H2 = decode("utf8", "Group1中所有样本的inclusion junction counts,生物学重复用逗号分隔");
my $I1 = decode("utf8", "SJC_SAMPLE_1"); my $I2 = decode("utf8", "Group1中所有样本的skipping junction counts,生物学重复用逗号分隔");
my $J1 = decode("utf8", "IJC_SAMPLE_2"); my $J2 = decode("utf8", "Group2中所有样本的inclusion junction counts,生物学重复用逗号分隔");
my $K1 = decode("utf8", "SJC_SAMPLE_2");my $K2 = decode("utf8", "Group2中所有样本的skipping junction counts,生物学重复用逗号分隔");
my $L1 = decode("utf8", "IncFormLen");my $L2 = decode("utf8", "Inclusion Form的长度，用于归一化表达量");
my $M1 = decode("utf8", "SkipFormLen");my $M2 = decode("utf8", "Skipping Form的长度，用于归一化表达量");
my $N1 = decode("utf8", "PValue");my $N2 = decode("utf8", "差异检验的P值");
my $O1 = decode("utf8", "FDR");my $O2 = decode("utf8", "校正后的P值");
my $P1 = decode("utf8", "IncLevel1");my $P2 = decode("utf8", "Group1中所有样本的inclusion level,生物学重复用逗号分隔，用IncRormLen归一化之后的表达量");
my $Q1 = decode("utf8", "IncLevel2");my $Q2 = decode("utf8", "Group2中所有样本的inclusion level,生物学重复用逗号分隔，用IncRormLen归一化之后的表达量");
my $R1 = decode("utf8", "IncLevelDifference");my $R2 = decode("utf8", "average(IncLevel1) - average(IncLevel2)");

my $S1 = decode("utf8", "A3SS");my $S2 = decode("utf8", "");
my $T1 = decode("utf8", "isoform_exon");my $T2 = decode("utf8", "完整的外显子对应的转录本和外显子编号，比如对应转录本NM_0001的第12号外显子，写作NM_0001(12)");
my $U1 = decode("utf8", "longExonStart_0base");my $U2 = decode("utf8", "完整的外显子的起始位置，坐标从0开始计算");
my $V1 = decode("utf8", "longExonEnd");my $V2 = decode("utf8", "完整的外显子的终止位置");
my $W1 = decode("utf8", "isoform_up");my $W2 = decode("utf8", "缺失掉的外显子对应的转录本和外显子编号");
my $X1 = decode("utf8", "shortES");my $X2 = decode("utf8", "缺失掉的外显子的起始位置");
my $Y1 = decode("utf8", "shortEE");my $Y2 = decode("utf8", "缺失掉的外显子的终止位置");
my $Z1 = decode("utf8", "isoform_down");my $Z2 = decode("utf8", "上下游的外显子对应的转录本和外显子编号");

my $A3 = decode("utf8", "flankingES");my $A4 = decode("utf8", "上游的最邻近的外显子的起始位置");
my $B3 = decode("utf8", "flankingEE");my $B4 = decode("utf8", "下游的最邻近的外显子的终止位置");

my $C3 = decode("utf8", "A5SS");my $C4 = decode("utf8", "");
my $D3 = decode("utf8", "isoform_exon");my $D4 = decode("utf8", "完整的外显子对应的转录本和外显子编号，比如对应转录本NM_0001的第12号外显子，写作NM_0001(12)");
my $E3 = decode("utf8", "longExonStart_0base");my $E4 = decode("utf8", "完整的外显子的起始位置");
my $F3 = decode("utf8", "longExonEnd");my $F4 = decode("utf8", "完整的外显子的终止位置");
my $G3 = decode("utf8", "isoform_up");my $G4 = decode("utf8", "缺失掉的外显子对应的转录本和外显子编号");
my $H3 = decode("utf8", "shortES");my $H4 = decode("utf8", "缺失掉的外显子的起始位置");
my $I3 = decode("utf8", "shortEE");my $I4 = decode("utf8", "缺失掉的外显子的终止位置");
my $J3 = decode("utf8", "isoform_down");my $J4 = decode("utf8", "上下游的外显子对应的转录本和外显子编号");
my $K3 = decode("utf8", "flankingES");my $K4 = decode("utf8", "上游的最邻近的外显子的起始位置");
my $L3 = decode("utf8", "flankingEE");my $L4 = decode("utf8", "下游的最邻近的外显子的终止位置");

my $M3 = decode("utf8", "MXE");my $M4 = decode("utf8", "");
my $N3 = decode("utf8", "isoform_exon1");my $N4 = decode("utf8", "第一个外显子对应的转录本和外显子编号");
my $O3 = decode("utf8", "1stExonStart_0base");my $O4 = decode("utf8", "第一个外显子的起始位置，坐标从0开始计算");
my $P3 = decode("utf8", "1stExonEnd");my $P4 = decode("utf8", "第一个外显子的终止位置");
my $Q3 = decode("utf8", "isoform_exon2");my $Q4 = decode("utf8", "第二个外显子对应的转录本和外显子编号");
my $R3 = decode("utf8", "2ndExonStart_0base");my $R4 = decode("utf8", "第二个外显子的起始位置，坐标从0开始计算");
my $S3 = decode("utf8", "2ndExonEnd");my $S4 = decode("utf8", "第二个外显子的终止位置");
my $T3 = decode("utf8", "isoform_up");my $T4 = decode("utf8", "上游最邻近的外显子对应的转录本和外显子编号");
my $U3 = decode("utf8", "upstreamES");my $U4 = decode("utf8", "上游最邻近的外显子的起始位置，坐标从0开始计算");
my $V3 = decode("utf8", "upstreamEE");my $V4 = decode("utf8", "上游最邻近的外显子的终止位置");
my $W3 = decode("utf8", "isoform_down");my $W4 = decode("utf8", "下游最邻近的外显子对应的转录本和外显子编号");
my $X3 = decode("utf8", "downstreamES");my $X4 = decode("utf8", "下游最邻近的外显子的起始位置，坐标从0开始计算");
my $Y3 = decode("utf8", "downstreamEE");my $Y4 = decode("utf8", "下游最邻近的外显子的终止位置");


my $Z3 = decode("utf8", "RI");my $Z4 = decode("utf8", "");
my $A5 = decode("utf8", "isoform_exon");my $A6 = decode("utf8", "保留的内含子将上下游两个外显子连接起来作为一个外显子进行转录，该外显子对应的转录本和外显子编号");
my $B5 = decode("utf8", "riExonStart_0base");my $B6 = decode("utf8", "保留的内含子将上下游两个外显子连接起来作为一个外显子进行转录，该外显子的起始位置，坐标从0开始计算");
my $C5 = decode("utf8", "riExonEnd");my $C6 = decode("utf8", "保留的内含子将上下游两个外显子连接起来作为一个外显子进行转录，该外显子的终止位置");
my $D5 = decode("utf8", "isoform_up");my $D6 = decode("utf8", "保留下来的内含子区域上游最邻近的外显子对应的转录本和外显子编号");
my $E5 = decode("utf8", "upstreamES");my $E6 = decode("utf8", "保留下来的内含子区域上游最邻近的外显子的起始位置，坐标从0开始计算");
my $F5 = decode("utf8", "upstreamEE");my $F6 = decode("utf8", "保留下来的内含子区域上游最邻近的外显子的终止位置");
my $G5 = decode("utf8", "isoform_down");my $G6 = decode("utf8", "保留下来的内含子区域下游最邻近的外显子对应的转录本和外显子编号");
my $H5 = decode("utf8", "downstreamES");my $H6 = decode("utf8", "保留下来的内含子区域下游最邻近的外显子的起始位置，坐标从0开始计算");
my $I5 = decode("utf8", "downstreamEE");my $I6 = decode("utf8", "保留下来的内含子区域下游最邻近的外显子的终止位置");


my $J5 = decode("utf8", "SE");my $J6 = decode("utf8", "");
my $K5 = decode("utf8", "isoform_exon");my $K6 = decode("utf8", "不参与转录的外显子对应的转录本和外显子编号");
my $L5 = decode("utf8", "exonStart_0base");my $L6 = decode("utf8", "不参与转录的外显子对应的起始位置，坐标从0开始计算");
my $M5 = decode("utf8", "exonEnd");my $M6 = decode("utf8", "不参与转录的外显子对应的终止位置");
my $N5 = decode("utf8", "isoform_up");my $N6 = decode("utf8", "上游最邻近的外显子对应的转录本和外显子编号");
my $O5 = decode("utf8", "upstreamES");my $O6 = decode("utf8", "上游最邻近的外显子的起始位置，坐标从0开始计算");
my $P5 = decode("utf8", "upstreamEE");my $P6 = decode("utf8", "上游最邻近的外显子的终止位置");
my $Q5 = decode("utf8", "isoform_down");my $Q6 = decode("utf8", "下游最邻近的外显子对应的转录本和外显子编号");
my $R5 = decode("utf8", "downstreamES");my $R6 = decode("utf8", "下游最邻近的外显子的起始位置，坐标从0开始计算");
my $S5 = decode("utf8", "downstreamEE");my $S6 = decode("utf8", "下游最邻近的外显子的终止位置");


$worksheet6->write_row( 0, 0, [$A1, $A2], $format);
$worksheet6->write_row( 1, 0, [$B1, $B2], $cell_format);
$worksheet6->write_row( 2, 0, [$C1, $C2], $cell_format);
$worksheet6->write_row( 3, 0, [$D1, $D2], $cell_format);
$worksheet6->write_row( 4, 0, [$E1, $E2], $cell_format);
$worksheet6->write_row( 5, 0, [$F1, $F2], $cell_format);
$worksheet6->write_row( 6, 0, [$G1, $G2], $cell_format);
$worksheet6->write_row( 7, 0, [$H1, $H2], $cell_format);
$worksheet6->write_row( 8, 0, [$I1, $I2], $cell_format);
$worksheet6->write_row( 9, 0, [$J1, $J2], $cell_format);
$worksheet6->write_row( 10, 0, [$K1, $K2], $cell_format);
$worksheet6->write_row( 11, 0, [$L1, $L2], $cell_format);
$worksheet6->write_row( 12, 0, [$M1, $M2], $cell_format);
$worksheet6->write_row( 13, 0, [$N1, $N2], $cell_format);
$worksheet6->write_row( 14, 0, [$O1, $O2], $cell_format);
$worksheet6->write_row( 15, 0, [$P1, $P2], $cell_format);
$worksheet6->write_row( 16, 0, [$Q1, $Q2], $cell_format);
$worksheet6->write_row( 17, 0, [$R1, $R2], $cell_format);
$worksheet6->write_row( 18, 0, [$S1, $S2], $cell_format);
$worksheet6->write_row( 19, 0, [$T1, $T2], $cell_format);
$worksheet6->write_row( 20, 0, [$U1, $U2], $cell_format);
$worksheet6->write_row( 21, 0, [$V1, $V2], $cell_format);
$worksheet6->write_row( 22, 0, [$W1, $W2], $cell_format);
$worksheet6->write_row( 23, 0, [$X1, $X2], $cell_format);
$worksheet6->write_row( 24, 0, [$Y1, $Y2], $cell_format);
$worksheet6->write_row( 25, 0, [$Z1, $Z2], $cell_format);

$worksheet6->write_row( 26, 0, [$A3, $A4], $cell_format);
$worksheet6->write_row( 27, 0, [$B3, $B4], $cell_format);
$worksheet6->write_row( 28, 0, [$C3, $C4], $cell_format);
$worksheet6->write_row( 29, 0, [$D3, $D4], $cell_format);
$worksheet6->write_row( 30, 0, [$E3, $E4], $cell_format);
$worksheet6->write_row( 31, 0, [$F3, $F4], $cell_format);
$worksheet6->write_row( 32, 0, [$G3, $G4], $cell_format);
$worksheet6->write_row( 33, 0, [$H3, $H4], $cell_format);
$worksheet6->write_row( 34, 0, [$I3, $I4], $cell_format);
$worksheet6->write_row( 35, 0, [$J3, $J4], $cell_format);
$worksheet6->write_row( 36, 0, [$K3, $K4], $cell_format);
$worksheet6->write_row( 37, 0, [$L3, $L4], $cell_format);
$worksheet6->write_row( 38, 0, [$M3, $M4], $cell_format);
$worksheet6->write_row( 39, 0, [$N3, $N4], $cell_format);
$worksheet6->write_row( 40, 0, [$O3, $O4], $cell_format);
$worksheet6->write_row( 41, 0, [$P3, $P4], $cell_format);
$worksheet6->write_row( 42, 0, [$Q3, $Q4], $cell_format);
$worksheet6->write_row( 43, 0, [$R3, $R4], $cell_format);
$worksheet6->write_row( 44, 0, [$S3, $S4], $cell_format);
$worksheet6->write_row( 45, 0, [$T3, $T4], $cell_format);
$worksheet6->write_row( 46, 0, [$U3, $U4], $cell_format);
$worksheet6->write_row( 47, 0, [$V3, $V4], $cell_format);
$worksheet6->write_row( 48, 0, [$W3, $W4], $cell_format);
$worksheet6->write_row( 49, 0, [$X3, $X4], $cell_format);
$worksheet6->write_row( 50, 0, [$Y3, $Y4], $cell_format);
$worksheet6->write_row( 51, 0, [$Z3, $Z4], $cell_format);


$worksheet6->write_row( 52, 0, [$A5, $A6], $cell_format);
$worksheet6->write_row( 53, 0, [$B5, $B6], $cell_format);
$worksheet6->write_row( 54, 0, [$C5, $C6], $cell_format);
$worksheet6->write_row( 55, 0, [$D5, $D6], $cell_format);
$worksheet6->write_row( 56, 0, [$E5, $E6], $cell_format);
$worksheet6->write_row( 57, 0, [$F5, $F6], $cell_format);
$worksheet6->write_row( 58, 0, [$G5, $G6], $cell_format);
$worksheet6->write_row( 59, 0, [$H5, $H6], $cell_format);
$worksheet6->write_row( 60, 0, [$I5, $I6], $cell_format);
$worksheet6->write_row( 61, 0, [$J5, $J6], $cell_format);
$worksheet6->write_row( 62, 0, [$K5, $K6], $cell_format);
$worksheet6->write_row( 63, 0, [$L5, $L6], $cell_format);
$worksheet6->write_row( 64, 0, [$M5, $M6], $cell_format);
$worksheet6->write_row( 65, 0, [$N5, $N6], $cell_format);
$worksheet6->write_row( 66, 0, [$O5, $O6], $cell_format);
$worksheet6->write_row( 67, 0, [$P5, $P6], $cell_format);
$worksheet6->write_row( 68, 0, [$Q5, $Q6], $cell_format);
$worksheet6->write_row( 69, 0, [$R5, $R6], $cell_format);
$worksheet6->write_row( 70, 0, [$S5, $S6], $cell_format);




$workbook->close();




sub read_file
{
	my $file = shift;

	my %meta  = ();
	open EXO, $file or die "Can't open $file!\n";
	while (<EXO>) {
		chomp;
		next if /^ID/;
		my @arr = split /\t/;
		s/"//g;
		my $index  = $#arr - 4;
		next if $arr[$index] > 0.05;
		$meta{$arr[0]} = $_;
	}
	close EXO;
	return %meta;
}

sub cal_overlap
{
	my $case1 = shift;
	my $case2 = shift;

	my %meta  = ();
	foreach my $x (keys %{$case1}) {
		next if not exists $case2->{$x};
		$meta{$x} = [$case1->{$x}, $case2->{$x}];
	}

	return %meta;
}
