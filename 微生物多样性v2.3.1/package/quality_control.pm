package quality_control;

use Parallel::ForkManager;

sub run
{
	my $metadata    = shift;
	my $base        = shift;

	my $inputpath   = qq{$metadata->{raw_data}};
	my $F_primer    = qq{$metadata->{F_primer}};
	my $R_primer    = qq{$metadata->{R_primer}};
	my $qcpath      = qq{$metadata->{intm_result}};
	my $projectpath = qq{$metadata->{result}};
	my $clean_reads = qq{$metadata->{Clean_reads}};

	my $seq_crumbs     = qq{$base->{seq_crumbs_bin}};
	my $trim_galore    = qq{$base->{trim_galore_bin}};
	my $flash2         = qq{$base->{flash2_bin}};
	my $fastx_toolkit  = qq{$base->{fastx_toolkit_bin}};
	my $mothur         = qq{$base->{mothur_bin}};
	my $usearch        = qq{$base->{usearch_bin}};
	my $Rscript        = qq{$base->{Rscript_bin}};
	my $scriptdir      = qq{$base->{scriptdir}};
	my $primerdb       = qq{$base->{primerdb}};

	system qq{mkdir -p $qcpath}      if not -d $qcpath;
	system qq{mkdir -p $projectpath} if not -d $projectpath;

	#pre_check
	die "ERROR: There is no $qcpath, Please Check It! $!\n" if not -d $qcpath;
	die "ERROR: There is no $inputpath, Please Check It! $!\n" if not -d $inputpath;
	
	seqcheck("$projectpath/Statistics/seqStat.xls", $clean_reads) if -e qq{$projectpath/Statistics/seqStat.xls};
	
	if (-e qq{$qcpath/qc.finish}){
		print qq{质量控制已经运行完成!\n 是否重新运行? [y/n]};
		my $option = <STDIN>;
		$option =~ s/[\r\n]//g;
		return if $option ne "y";
		system qq{rm $qcpath/qc.finish} if -e qq{$qcpath/qc.finish};
		system qq{rm $qcpath/qc.log} if -e qq{$qcpath/qc.log};
		system qq{rm -rf $qcpath/raw_data} if -d qq{$qcpath/raw_data};
		system qq{rm -rf $qcpath/qc} if -d qq{$qcpath/qc};
	}

	my $raw_datapath = qq{$qcpath/raw_data};
	system qq{rm -rf $qcpath/raw_data} if -d qq{$qcpath/raw_data};
	system qq{mkdir -p $raw_datapath} if not -d $raw_datapath;
	system qq{ln -s $inputpath/*.gz $raw_datapath} if not -e qq{$raw_datapath/*.gz};
	
	my @samples = `ls $raw_datapath | grep _R1.fastq.gz | sed 's/_R1.fastq.gz//g'`;
	
	my $num = scalar @samples;

	#check
	my @cs = check($raw_datapath, $num, $F_primer, $R_primer, $primerdb);

	# Prepare primer
	my $QCpath = qq{$qcpath/qc};
	system qq{mkdir -p $QCpath} if not -d $QCpath;
	system qq{echo "forward\t$F_primer" > $QCpath/primer};
	system qq{echo "reverse\t$R_primer" >> $QCpath/primer};
	system qq{echo "forward\t$R_primer" > $QCpath/rev.primer};
	system qq{echo "reverse\t$F_primer" >> $QCpath/rev.primer};

	# Statistics
	my $Statisticspath = qq{$projectpath/Statistics};
	system qq{mkdir -p $Statisticspath} if not -d $Statisticspath;
	open Sta, qq{>$Statisticspath/seqStat.xls} or die "Can not write into $Statisticspath, please check it!\n";
	# Run
	my $MAX_PROCESSES = 18;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);
	foreach my $sample (@samples){
		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;
		chomp $sample;
		my $Fastq1 = qq{$raw_datapath/$sample\_R1.fastq.gz};
		my $Fastq2 = qq{$raw_datapath/$sample\_R2.fastq.gz};
		my $sampledir = qq{$QCpath/$sample};
		system qq{mkdir -p $sampledir} if not -d $sampledir;

		#################precheck
		if (-e "$sampledir/$sample.raw.stat" and "$sampledir/$sample.primer.fasta" and "$sampledir/$sample.clean.fasta" and "$sampledir/$sample.clean.stat"){

			my ($arr, $totalrawreads) = seq_stat("$sampledir/$sample.raw.stat");
			my @arr = split /\t/, $arr;

			my $prim = `less $sampledir/$sample.primer.fasta | grep '>' | wc -l`;
			chomp $prim;

			my ($arr2, $totalcleanreads) = seq_stat("$sampledir/$sample.clean.stat");
			my @arr2 = split /\t/, $arr2;

			if ($arr[2] > 70 and $prim > 0.9*$clean_reads and $arr2[2] > 80){
				system qq{touch $sampledir/$sample.finish} if not -e "$sampledir/$sample.finish";
			}else{
				system qq{rm -rf $sampledir/$sample.finish} if -e "$sampledir/$sample.finish";
			}
		}
		################precheck

		if (not -e "$sampledir/$sample.finish"){
			### raw seq Statistics (R1)
			system qq{$seq_crumbs calculate_stats $raw_datapath/$sample\_R1.fastq.gz > $sampledir/$sample.raw.stat};

			# 去除adapter序列;去除末端质量低于20的序列;去除长度小于100的序列
			system qq{$trim_galore trim_galore --dont_gzip --suppress_warn --no_report_file --output_dir $sampledir -q 20 --phred33 --length 100 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -a2 CTGTCTCTTATACACATCTGACGCTGCCGACGA --paired $Fastq1 $Fastq2 &>> $qcpath/qc.log};
			
			# 双端合并
			system qq{$flash2 flash2 $sampledir/$sample\_R1_val_1.fq $sampledir/$sample\_R2_val_2.fq -Q 20 -C 90 -d $sampledir -m 10 -M 250 -t 10 2>&1 | tee $sampledir/flash.log &>> $qcpath/qc.log};
			system qq{mv $sampledir/out.extendedFrags.fastq $sampledir/$sample.merge.fastq};
			
			# fastq2fasta
			system qq{$fastx_toolkit fastq_to_fasta -n -v -i $sampledir/$sample.merge.fastq -o $sampledir/$sample.merge.fasta -Q 33 &>> $qcpath/qc.log};

			### 反向序列去引物及质控
			system qq{$mothur mothur "#trim.seqs(fasta=$sampledir/$sample.merge.fasta, oligos=$QCpath/rev.primer, maxambig=1, pdiffs=2, processors=30)" &>> $qcpath/qc.log};
			system qq{mv $sampledir/$sample.merge.trim.fasta $sampledir/$sample.merge.trim1.fasta};
			
			# fasta2fastq
			my $fastq = qq{$sampledir/$sample.merge.fastq};
			my $trim1 = qq{$sampledir/$sample.merge.trim1.fasta};
			fasta2fastq($fastq, $trim1, $sampledir);

			# filter
			system qq{$usearch /usearch -fastq_filter $sampledir/$sample.merge.trim1.fastq -fastaout $sampledir/$sample.rev.clean.fasta -fastqout $sampledir/$sample\.rev.clean.fastq -sample $sample -fastq_maxee 2 -fastq_minlen 100 &>> $qcpath/qc.log};
			
			# reverse_complement for dereplication
			system qq{$fastx_toolkit fasta_formatter -i $sampledir/$sample.rev.clean.fasta -o $sampledir/$sample.rev.oneline.clean.fasta};
			system qq{$fastx_toolkit fastx_reverse_complement -i $sampledir/$sample.rev.oneline.clean.fasta -o $sampledir/$sample.for2.clean.fasta -Q 33};
			system qq{$fastx_toolkit fastx_reverse_complement -i $sampledir/$sample.rev.clean.fastq -o $sampledir/$sample.for2.clean.fastq -Q 33};
			
			### 正向序列去引物及质控
			system qq{$mothur mothur "#trim.seqs(fasta=$sampledir/$sample.merge.fasta, oligos=$QCpath/primer, maxambig=1, pdiffs=2, processors=10)" &>> $qcpath/qc.log};
			system qq{mv $sampledir/$sample.merge.trim.fasta $sampledir/$sample.merge.trim2.fasta};
			
			# fasta2fastq
			my $trim2 = qq{$sampledir/$sample.merge.trim2.fasta};
			fasta2fastq($fastq, $trim2, $sampledir);
			
			# filter
			system qq{$usearch /usearch -fastq_filter $sampledir/$sample.merge.trim2.fastq -fastaout $sampledir/$sample.for.clean.fasta -fastqout $sampledir/$sample.for1.clean.fastq -sample $sample -fastq_maxee 2 -fastq_minlen 100 &>> $qcpath/qc.log};
			system qq{$fastx_toolkit fasta_formatter -i $sampledir/$sample.for.clean.fasta -o $sampledir/$sample.for1.clean.fasta};

			### 合并
			# primer
			system qq{cat $sampledir/$sample.merge.trim1.fasta $sampledir/$sample.merge.trim2.fasta > $sampledir/$sample.primer.fasta};
			system qq{cat $sampledir/$sample.merge.trim1.fastq $sampledir/$sample.merge.trim2.fastq > $sampledir/$sample.primer.fastq};
			# 合并
			system qq{cat $sampledir/$sample.for1.clean.fasta $sampledir/$sample.for2.clean.fasta > $sampledir/$sample.clean.fasta};
			system qq{sed -i 's/sample/barcodelabel/g' $sampledir/$sample.clean.fasta};
			system qq{cat $sampledir/$sample.for1.clean.fastq $sampledir/$sample.for2.clean.fastq > $sampledir/$sample.clean.fastq};
			system qq{sed -i 's/sample/barcodelabel/g' $sampledir/$sample.clean.fastq};
			system qq{rm -f $sampledir/*for* $sampledir/*rev* $sampledir/*trim* $sampledir/*scrap* $sampledir/out*};
			
			### clean seq Statistics
			system qq{$seq_crumbs calculate_stats $sampledir/$sample.clean.fastq > $sampledir/$sample.clean.stat};
		}
		### seq Statistics
		my ($arr, $totalrawreads) = seq_stat("$sampledir/$sample.raw.stat");
		my @arr = split /\t/, $arr;
		my $stat = qq{$sample\t$arr};

		my $merge = `less $sampledir/$sample.merge.fastq | grep '^+\$' | wc -l`;
		chomp $merge;
		my $perc_merge = sprintf "%0.2f",100*$merge/$totalrawreads;

		my $prim = `less $sampledir/$sample.primer.fasta | grep '>' | wc -l`;
		chomp $prim;

		my ($arr2, $totalcleanreads) = seq_stat("$sampledir/$sample.clean.stat");
		my @arr2 = split /\t/, $arr2;
		my $perc_clean = sprintf "%0.2f",100*$totalcleanreads/$totalrawreads;
		my $stat2 = qq{$arr2[0]\t$perc_clean\t$arr2[1]\t$arr2[2]};

		print Sta qq{$stat\t$merge\t$perc_merge\t$prim\t$stat2\n};
		
		if ($arr[2] > 70 and $prim > 0.9*$clean_reads and $arr2[2] > 80){
			system qq{touch $sampledir/$sample.finish} if not -e "$sampledir/$sample.finish";
		}else{
			system qq{rm -rf $sampledir/$sample.finish} if -e "$sampledir/$sample.finish";
		}

		$pm->finish; # Terminates the child process
	}
	$pm->wait_all_children;  # blocks until all forked processes have exited
	system qq{rm mothur.*.logfile} if -e qq{mothur.*.logfile};
	# add title to seqStat
	system qq{sed -i '1iSamples\tRaw_reads\tQ20(%)\tQ30(%)\tMerge\tMerge(%)\tPrimer\tClean_reads\tClean_reads(%)\tQ20(%)\tQ30(%)' $Statisticspath/seqStat.xls};
	close Sta;

	#样本质控结果
	seqcheck("$Statisticspath/seqStat.xls", $clean_reads);

	# generate reads.clean.fasta
	my $str = (join "\\|", @cs)."\\|"." ";
	system qq{find $QCpath -name '*clean.fasta' | grep -v '$str' | xargs cat > $projectpath/reads.clean.fasta};
	
	# seq length distribution
	my $fasta = qq{$projectpath/reads.clean.fasta};
	my ($name) = $fasta =~ /([^\/]*)\.f/;
	system qq{$seq_crumbs calculate_stats $fasta > $Statisticspath/$name.stat};

	seq_distribution($name, $Statisticspath);

	#seq-distribut.pdf
	system qq{$Rscript $scriptdir/seq-distribut.r $Statisticspath/reads.clean.length.distribution.xls $Statisticspath &>> $qcpath/qc.log};
	
	system qq{touch $qcpath/qc.finish};
	print qq{质量控制已经运行完成!\n};
}


sub check
{	
	my $raw_datapath = shift;
	my $num          = shift;
	my $F_primer     = shift;
	my $R_primer     = shift;
	my $primerdb     = shift;

	# control sample name
	print "======= 输入对照样本的名称，如果没有则输入 NONE, 多个对照样本用空格区分 =======\n";
	my $opt = <STDIN>;
	$opt =~ s/[\r\n]//g;
	my @cs = "NONE";
	if (uc($opt) ne "NONE"){
		@cs = split /\s+/, $opt;
		foreach my $c (@cs){
			die "Can not find $raw_datapath/$c\_R*.fastq.gz, please check it!\n" if not -e qq{$raw_datapath/$c\_R1.fastq.gz};
		}
		$num = $num - (scalar @cs);
	}

	# check the number
	print "======= 本次项目共有 $num 个样本(不包括对照) =======\n  确定? [y/n]";
	my $option = <STDIN>;
	$option =~ s/[\r\n]//g;
	exit if $option ne "y";

	# check the primer
	open PRIMER, $primerdb or die "Can not open $primerdb, please check it!\n";
	my $check = "no";
	while(<PRIMER>){
		chomp;
		next if /^\s*$/;
		my $str1 = $_;
		my $str2 = <PRIMER>;
		my $F;
		my $R;
		my @data1 = split /\s+/, $str1;
		my @data2 = split /\s+/, $str2;
		if ($data1[1] eq "F" and $data2[1] eq "R" and $data1[0] eq $data2[0]){
			$F = $data1[2];
			$R = $data2[2];
			if ($F eq $F_primer and $R eq $R_primer){

				print "======== 您输入的序列为 ".$data1[0]." 的引物序列 ========\n  确定? [y/n]";
				my $option = <STDIN>;
				$option =~ s/[\r\n]//g;
				exit if $option ne "y";
				print "================= Start QC =================\n";

				$check = "yes";
			}
		}elsif ($data2[1] eq "F" and $data1[1] eq "R" and $data1[0] eq $data2[0]){
			$F = $data2[2];
			$R = $data1[2];
			
			if ($F eq $F_primer and $R eq $R_primer){
				print "======== 您输入的序列为 ".$data1[0]." 的引物序列 ========\n  确定? [y/n]";
				my $option = <STDIN>;
				$option =~ s/[\r\n]//g;
				exit if $option ne "y";
				print "================= Start QC =================\n";

				$check = "yes";
			}
		}else{
			die "======== $primerdb 中 $data1[0] 的引物出错 ========\n";
		}
		
		last if $check eq "yes";
	}
	close PRIMER;

	if($check eq "no"){
		print "======= 无法在 $primerdb 中找到输入的引物序列,确定引物序列的准确性并继续吗？ =======\n  继续? [y/n]";
		my $option = <STDIN>;
		$option =~ s/[\r\n]//g;
		exit if $option ne "y";
		print "================= Start QC =================\n";
	}

	return @cs;
}


sub seq_stat
{
	my $seq = shift;

	open FILE, $seq or die "Can not open $seq, please check it!\n";

		my $totbases;
		my $totreads;
		my $average;
		my $Q20;
		my $Q30;
		my $a = 1;
		while(<FILE>){
			
			chomp;
			next if /^\s*$/;
			
			if (/tot. residues/){
				($totbases) = $_ =~ /^tot. residues: (.*)/;
			}
			if (/num. seqs/){
				($totreads) = $_ =~ /^num. seqs.: (.*)/;
			}
			if (/average/ and $a == 1){
				($average) = $_ =~ /^average: (.*)/;
				$a = $a + 1;
			}
			if (/Q20/){
				($Q20) = $_ =~ /^Q20: (.*)/;
			}
			if (/Q30/){
				($Q30) = $_ =~ /^Q30: (.*)/;
			}	
		}
		
		my $stat = qq{$totreads\t$Q20\t$Q30};
	close FILE;

	return ($stat, $totreads);
}

sub seqcheck
{
	my $seq = shift;
	my $clean_reads = shift;

	open FILE, $seq;
	my $i = 0;
	my $j = 0;
	my $sum = 0;
	my @arr;

	while(<FILE>){
		chomp;
		next if /Samples/;
		
		my @str = split /\t/, $_;
		
		if ($str[3] > 70 and $str[6] > 0.9*$clean_reads and $str[10] > 80){
			$i++;
			
		}else{
			$j++;
			push @arr,$str[0];
			#print qq{样本$str[0]\tRaw_data_Q30(%):$str[3]\tClean_reads($str[7])\tClean_data_Q30(%):$str[10]\n};
			print qq{$_\n};

		}
	}
	$sum = $i + $j;
	if ($j >0){
		print qq{样本质控完成,共计$sum个样本,$j个样本质控不合格\n 是否继续? [y/n]};
		my $option = <STDIN>;
		$option =~ s/[\r\n]//g;
		exit if $option ne "y";
	}else{
		print qq{prefect！共计$i个样本,质控全部合格\n};
	}
	return @arr;
	close FILE;
}

sub seq_distribution
{
	my $name            = shift;
	my $Statisticspath  = shift;

	open FILE, qq{$Statisticspath/$name.stat} or die "Can not open $outputpath/$name.stat, please check it!\n";
	open OUT, qq{>$Statisticspath/$name.length.distribution.xls} or die "Can not write into $outputpath, please check it!\n";

	while(<FILE>){

		chomp;
		last if /Quality stats and distribution/;
		
		next if not /^\[/;

		my ($start) = $_ =~ /\[(\d+)\s+,/;
		my ($end) = $_ =~ /,\s+(\d+)\[/;
		my ($len) = $_ =~ /(\d+)\):/;
		
		if ($start + 1 == $end){
			
			print OUT $start."\t".$len."\n";
		
		}else{
		
			print OUT $start."-".$end."\t".$len."\n";

		}
	}
}

sub fasta2fastq
{
	
	my $fastq     = shift;
	my $fasta     = shift;
	my $sampledir = shift;
	my ($name) = $fasta =~ /([^\/]*)\.fa/;
	# read in fasta
	$/ = ">"; 
	open FA, $fasta, or die "Can not open $fasta, please check it!\n";
	my %FASTA;
	while(<FA>){
		
		chomp;
		next if /^\s*$/;
		
		my @data = split /\n/;
		my $id = (split /\s+/, $data[0])[0];
		$id =~ s/:/_/g;

		my $seq = $data[1];
		$FASTA{$id} = $seq;
		
	}
	$/ = "\n";
	close FA;
	#print "[WARNING]: your sequence names contained ':'.  I changed them to '_' to avoid problems in your downstream analysis.\n";

	# read in fastq
	open FQ, $fastq or die "Can not open $fastq, please check it!\n";
	open OTU, qq{>$sampledir/$name\.fastq} or die "Can not write into $sampledir/$name\.fastq, please check it !\n";
	while(<FQ>){

		chomp;
		next if /^\s*$/;
		my $ids = $_;
		my $seq = <FQ>;
		my $cha = <FQ>;
		my $qua = <FQ>;

		my $id = (split /\s+/, $ids)[0];
		$id =~ s/:/_/g;
		$id =~ s/^\@//g;

		if (exists $FASTA{$id}){

			# 子序列（已去除引物）在原始序列中的位置, 提取其对应的质量值
			my $index = index($seq, $FASTA{$id});

			if ($index == -1){		# we need reverse complement this sequence

				#my ($rev) = `perl /home/panrf/CommonScripts/reverse_complement_seq.pl $FASTA{$id}`;
				
				my $dna = $FASTA{$id};
				# reverse the DNA sequence
				my $rev = reverse($dna);
				# complement the reversed DNA sequence
				$rev =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

				chomp $rev;
				$index = index($seq, $rev);

			}

			my $quality = substr($qua, $index, length($FASTA{$id}));

			print OTU qq{\@$id\n$FASTA{$id}\n+\n$quality\n};
			
		}

	}

	close FQ;
	close OTU;
}

1;
