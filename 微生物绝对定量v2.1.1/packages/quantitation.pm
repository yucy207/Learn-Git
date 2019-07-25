package quantitation;

use Sort::Key::Natural qw(natsort);
use Parallel::ForkManager;

sub run
{
	my $metadata    = shift;
	my $base        = shift;
	# my $quantitation = qq{$metadata->{quantitation}};
	my $level 		 = qq{$metadata->{level}};
	my $IR_copeis    = qq{$metadata->{IR_copeis}};
	my $DNA_quantity = qq{$metadata->{DNA_quantity}};
	my $path         = qq{$metadata->{intm_result}};
	my $result       = qq{$metadata->{result}};

	system "sed -i 's/\r//g' $IR_copeis";
	system "sed -i 's/\r//g' $DNA_quantity";

	##绝对定量计算
	
	Internal_Reference($level, $IR_copeis, $DNA_quantity, $path, $result, $base);

	print qq{绝对拷贝数计算完成\n};

}


#################################################子程序

sub Internal_Reference {

	my $level        = shift;
	my $IR_copeis    = shift;
	my $DNA_quantity = shift;
	my $Path         = shift;
	my $result       = shift;
	my $base         = shift;	
	my $tmpPath      = qq{$Path/tmp};
	my $inter_ref_db = qq{$base->{SPIKE_IN}{fasta}};
	my $scriptdir    = qq{$base->{scriptdir}/sub_script};
	my $Rscript_xlsx = qq{$base->{Rscript_xlsx}};

	# return if -e qq{$tmpPath/Internal_Reference.finish};
	
	my $Back_Before_IR_Detection = "$tmpPath/Back_Before_IR_Detection";   ###备份数据文件夹
	my $IR_Detection = "$tmpPath/IR_Detection";                           ###存储IR检测结果

	system qq{mkdir -p $Back_Before_IR_Detection} if not -d $Back_Before_IR_Detection;
	system qq{mkdir -p $IR_Detection} if not -d $IR_Detection;

	###原始数据备份
	single_repseq_fasta("$tmpPath/otu.repseq.fasta", "$Back_Before_IR_Detection/otu.repseq.fasta") if not -e qq{$Back_Before_IR_Detection/otu.repseq.fasta}; ###代表序列备份并转成单行
	system qq{cp $tmpPath/otu.tax.0.03.xls  $Back_Before_IR_Detection/otu.tax.0.03.xls} if not -e qq{$Back_Before_IR_Detection/otu.tax.0.03.xls};
	system qq{cp $tmpPath/map.ID            $Back_Before_IR_Detection/map.ID} if not -e qq{$Back_Before_IR_Detection/map.ID};
	system qq{cp $tmpPath/map.txt           $Back_Before_IR_Detection/map.txt} if not -e qq{$Back_Before_IR_Detection/map.txt};
	system qq{cp $tmpPath/reads.clean.fasta $Back_Before_IR_Detection/bakIR_reads.clean.fasta} if not -e qq{$Back_Before_IR_Detection/bakIR_reads.clean.fasta};

	###（1）查找内参序列对应OTU
	my $repseq           = "$Back_Before_IR_Detection/otu.repseq.fasta";
	my $inter_otu        = "$IR_Detection/IR_otu_list.xls";  ###储存匹配信息

	system qq{/home/genesky/software/blast+/2.7.1/bin/blastn -query $repseq -db $inter_ref_db -evalue 1e-50 -outfmt 6 -max_target_seqs 100000000 -perc_identity 100 -num_threads 8 > $inter_otu};

	###（2）提取内参序列对应的otuID, 
	my %hash;
	open INTER_OTU, $inter_otu;
	while (<INTER_OTU>) {
		chomp;
		my @arr = split /\t/, $_;
		$hash{$arr[0]} = $arr[1] if ($arr[2] == 100 && $arr[10] == 0);  ##$arr[2]为score值，满分100，$arr[10]为E值，最低0
	}
	close INTER_OTU;

	####（3）删除otu.tax.0.03.xls中内参otu
	####     标注with.IR.otu.tax.0.03.xls中内参otu，用于计算Absolute_abundance
	my $otutable         = "$Back_Before_IR_Detection/otu.tax.0.03.xls";
	my $otu_no_inter_ref = "$IR_Detection/no.IR.otu.tax.0.03.xls";
	my $otu_new          = "$IR_Detection/with.IR.otu.tax.0.03.xls";

	open OTU_O, $otutable;
	open OTU_N,">$otu_new";
	open OTU_NO_IR, ">$otu_no_inter_ref";
	while (<OTU_O>) {
		chomp;
		my @otu = split /\t/, $_;
		if (exists $hash{$otu[0]}) {
			$otu[0] = $hash{$otu[0]};
			print OTU_N (join "\t", @otu)."\n";
		} else {
			print OTU_NO_IR "$_\n";
			print OTU_N "$_\n";
		}
	}
	close OTU_O;
	close OTU_N;
	close OTU_NO_IR;


	###（4）删除otu.repseq.fasta中的内参序列
	open SEQ_O, $repseq;
	open SEQ_N,">$IR_Detection/no.IR.otu.repseq.fasta";
	while (my $line1 = <SEQ_O>) {
		my $line2 = <SEQ_O>;
		my ($otu_id) = $line1 =~  /^>(\w+)/;
		next if exists $hash{$otu_id};
		print SEQ_N "$line1$line2";
	}
	close SEQ_O;
	close SEQ_N;


	###（5）计算IR比例
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.IR.persent.stat.r $otu_new $IR_Detection};

	# ###（6）去除reads.clean.fasta中的内参序列
	##提取内参OTU对应的序列编号
	my $mapID = "$Back_Before_IR_Detection/map.ID";
	my %mapID;
	open MAPID, $mapID;
	while (<MAPID>) {
		chomp;
		my @arr1 = split /\t/, $_;
		if (exists $hash{$arr1[1]}){
			$mapID{$arr1[0]} = $arr1[1];
		}
	}
	close MAPID;

	my %clean_id;
	my $map = "$Back_Before_IR_Detection/map.txt";  ##对应clean.fasta
	open MAP, $map;
	open READ_ID,">$IR_Detection/readID_IR.txt";
	while (<MAP>) {
		chomp;
		my @arr1 = split /\t/, $_;
		if (exists $mapID{$arr1[1]}){
			$clean_id{$arr1[0]} = "-";
			print READ_ID "$_\n";
		}
	}
	close MAP;
	close READ_ID;

	###删除reads.clean.fasta中的内参序列
	my $reads_clean = "$Back_Before_IR_Detection/bakIR_reads.clean.fasta";
	open CLEAN_O, $reads_clean;
	open CLEAN_N, ">$IR_Detection/no.IR.reads.clean.fasta";
	while (my $line3 = <CLEAN_O>) {
		my $line4 = <CLEAN_O>;
		my ($id) = $line3 =~ /^>(\S+)/;
		next if exists $clean_id{$id};
		print CLEAN_N "$line3$line4";
	}
	close CLEAN_O;
	close CLEAN_N;
	

	###（7）计算OTU绝对丰度
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Absolute_Abundance.r $otu_new $IR_copeis $tmpPath >> $tmpPath/IRstat.log};
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.OTU.Copies.unit.DNA.r $tmpPath/AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $DNA_quantity $tmpPath/AA_Calculate >> $tmpPath/IRstat.log};

	###若提供提取DNA量和样本量
	my $ncol = `head -n 2 $DNA_quantity |awk 'NR==2 {print NF}'`;
	$ncol =~ s/[\r\n]//g;
	# print "$ncol\n";
	if ($ncol == 4) {
		system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.OTU.Copies.unit.sample.r $tmpPath/AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $DNA_quantity $tmpPath/AA_Calculate >> $tmpPath/IRstat.log};

	}elsif($ncol > 4){
		die "无法计算样本水平拷贝数，请检查文件 DNA_quantity.txt 的数据是否准确\n";
	}

	# print "无法计算样本水平拷贝数，请检查文件 DNA_quantity.txt 的数据是否准确\n" if $level ne "DNA" and not -d "$tmpPath/AA_Calculate/otu_copies_unit_sample.xls";
	# if ($level eq "DNA") {
	# 	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Absolute_Abundance.r $otu_new $IR_copeis $tmpPath >> $tmpPath/IRstat.log};
	# 	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.OTU.Copies.unit.DNA.r $tmpPath/AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $DNA_quantity $tmpPath/AA_Calculate >> $tmpPath/IRstat.log};
	# }else{
	# 	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Absolute_Abundance.r $otu_new $IR_copeis $tmpPath >> $tmpPath/IRstat.log};
	# 	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.OTU.Copies.unit.DNA.r $tmpPath/AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $DNA_quantity $tmpPath/AA_Calculate >> $tmpPath/IRstat.log};
	# 	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.OTU.Copies.unit.sample.r $tmpPath/AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $DNA_quantity $tmpPath/AA_Calculate >> $tmpPath/IRstat.log};
	# }
	

	###替换原始数据
	system qq{cp $IR_Detection/no.IR.reads.clean.fasta        $result/reads.clean.fasta};   ###spike-in 序列已删除，用于提供给客户



	# =============== 此结果用于community/Gene_Copies_Analysis的基因拷贝数矫正分析，因为耗时长，所以放在这里 =============== 
	system qq{perl $scriptdir/Absolute_Abundance/16S.Community.rrnDB.otu_copy_count.pl $tmpPath/otu.repseq.fasta $tmpPath}; 


}


sub single_repseq_fasta
{
 
	my $infile = shift;
	my $single = shift;
	open INPUT, "$infile";
	open OUTPUT, ">$single";

	my $tmp;
	while (<INPUT>) {
		chomp;
		$tmp = join ";", $tmp, $_;
		
	}
	close INPUT;

	my @arry = split />/, $tmp;
	for my $x (1..@arry-1){
		my @fasta = split /;/, $arry[$x];
		my $otuid = $fasta[0];
		my $seq;
		for my $y (1..@fasta-1){
			$seq   = $seq.$fasta[$y];
		}
		print OUTPUT ">".$otuid."\n".$seq."\n";
		
	}
	close OUTPUT;

}

1

