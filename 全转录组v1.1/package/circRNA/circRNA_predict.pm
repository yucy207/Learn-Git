package package::circRNA::circRNA_predict;
use List::Util qw/max min/;
use Parallel::ForkManager;

sub run
{
	my $metadata  = shift;
	my $base      = shift;
	my @samples   = split /,/, $metadata->{'samples'};
	my $organ     = qq{$metadata->{organ}};

	my $trim      = qq{$metadata->{project}/trim/trim_galore};
	my $circrna   = qq{$metadata->{project}/circRNA/circRNA_predict};
	my $report    = qq{$metadata->{report}/05_circRNA_Analysis/01_circRNA_Prediction};

	my $util      = qq{$base->{util}};
	my $rscript   = qq{$base->{rscript_bin}};
	my $bwa       = qq{$base->{bwa_bin}};
	my $ciri_bin  = qq{$base->{ciri_bin}};

	my $ref_fasta = qq{$base->{$organ}{genome_fasta}};
	my $ref_gtf   = qq{$base->{$organ}{genome_mRNA_gtf}};
	my $ref_gtf_1 = qq{$base->{$organ}{circ_ref_gtf}};
	my $circBase  = qq{$base->{$organ}{circbase_db}};
	my $bwa_index = qq{$base->{$organ}{bwa_index}};

	system qq{mkdir -p $circrna}        if not -d $circrna;
	system qq{mkdir -p $circrna/log}    if not -d qq{$circrna/log};
	system qq{mkdir -p $circrna/run}    if not -d qq{$circrna/run};
	system qq{mkdir -p $circrna/result} if not -d qq{$circrna/result};

	system qq{mkdir -p $report}         if not -d $report;
	system qq{mkdir -p $report/Predict_circRNA_length_hist_Figure} if not -d qq{$report/Predict_circRNA_length_hist_Figure};
	system qq{mkdir -p $report/Predict_circRNA_Raw_Result}         if not -d qq{$report/Predict_circRNA_Raw_Result};


	my @res_samples = res_check(qq{$circrna/result}, \@samples);

	if (scalar @res_samples == 0) {
	 	print qq{circRNA 预测分析已经完成!\n};
	 	#return;
	}

	pre_check($metadata, $base);

	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@res_samples) {

		my $pid = $pm->start and next;	
		# CIRI
		my $bwa_map = qq{$bwa mem -T 19 -t 5 $bwa_index $trim/$x/$x\_R1_val_1.fq.gz $trim/$x/$x\_R2_val_2.fq.gz > $circrna/result/$x.sam\n};
		my $ciri    = qq{$ciri_bin CIRI2.pl -T 20 -I $circrna/result/$x.sam -O $circrna/result/$x.predict_circRNA.xls -F $ref_fasta -A $ref_gtf\n};
		my $cp      = qq{perl $util/convert_txt2excel.pl $circrna/result/$x.predict_circRNA.xls $report/Predict_circRNA_Raw_Result/$x.predict_circRNA.xlsx};
		my $rm      = qq{rm $circrna/result/$x.sam\n};
		my $finish  = qq{touch $circrna/result/$x.finish};

		open SAVE, qq{>$circrna/run/$x.sh} or die "Can't open $circrna/run/$x.sh\n";
		print SAVE qq{$bwa_map\n$ciri\n};
		print SAVE qq{$cp\n$rm\n$finish\n};
		close SAVE;

		system qq{bash $circrna/run/$x.sh &> $circrna/log/$x.log};

		$pm->finish;
	}

	$pm->wait_all_children;

	################################### report ###################################
	my $excel     = qq{$report/circRNA.Summary.xlsx};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);

	# add_worksheet:Summary
	my $worksheet = $workbook->add_worksheet("Summary");
	my $title     = qq{sample\tcircRNA Count\tMax Length\tMin Length\tAverage Length};
	my @head      = split /\t/, $title;
	$worksheet->write_row( 0, 0, \@head, $format{'title'});

	my %cnt         = ();
	my %length      = ();
	my %all_circRNA = ();
	foreach my $x (@samples) {
		open EXO, qq{$circrna/result/$x.predict_circRNA.xls} or die "Can't open $circrna/result/$x.predict_circRNA.xls\n";
		open SAVE, qq{>$circrna/result/$x.length.xls} or die "Can't open $circrna/result/$x.length.xls\n"; 
		while (<EXO>) {
			chomp;
			my @arr = split /\t/;
			next if /circRNA_start/;
			my $id  = qq{$arr[0]:$arr[10]};
			$all_circRNA{$id}{$x} = $_;
			$cnt{$x}++;
			my $length = $arr[3] - $arr[2];
			push @{$length{$x}}, $length;
			print SAVE qq{$length\n};
		}
		close EXO;
		close SAVE;
	}

	my $row = 1;
	foreach my $x (@samples) {

		my $max_length = max(@{$length{$x}});
		my $min_length = min(@{$length{$x}});
		my $mean_length = cal_mean(@{$length{$x}});
		my $line =  qq{$x\t$cnt{$x}\t$max_length\t$min_length\t$mean_length};
		my @result = split /\t/, $line;
		$worksheet->write_row($row, 0, \@result, $format{'normal'});
		$row++;
		# length_hist
		system qq{$rscript Rscript $util/ggplot2_hist.R $circrna/result/$x.length.xls $report/Predict_circRNA_length_hist_Figure/$x.circRNA.lenght.hist.pdf &> /home/tmp/$x.log\n};
	
	}

	# get all circRNA bed 
	open BED, qq{>$circrna/result/circRNA.bed} or die "Can't open $circrna/result/circRNA.bed!\n";
	foreach my $x (keys %all_circRNA) {
		my ($chr, $start, $end, $strand) = split /[:\|]+/, $x;
		print BED qq{$chr\t$start\t$end\t$x\t.\t$strand\n};

	}
	close BED;

	system qq{perl $util/circRNA_splice_seq_annotation.pl -gtf $ref_gtf_1 -bed $circrna/result/circRNA.bed -fasta $ref_fasta -out_dir $circrna/result};
	system qq{cp $circrna/result/circrnas_splice_transcript.fasta $report};

	my %hash       = parse_fasta(qq{$circrna/result/circrnas_splice_transcript.fasta});
	my %db         = parse_db($circBase) if -e $circBase;
	my %trans_type = parse_type(qq{$circrna/result/fmt.overlap.xls});

	# add_worksheet:circRNA_predict
	my $row1 = 1;
	my $worksheet1 = $workbook->add_worksheet("circRNA_predict");
	my $title1     = qq{circRNA_ID\tcircBase_ID\tchr\tstart\tend\tstrand\tgenomic length\tspliced length\tsamples\tcounts\tcircRNA_type\tbest_transcript\tgene symbol\n};
	
	open SUM, qq{>$circrna/result/Summary.xls} or die "Can't open $circrna/result/Summary.xls!\n";	
	my $title2 = qq{circRNA_ID\tcircBase_ID\tchr\tstart\tend\tstrand\tgenomic length\tspliced length\tcircRNA_type\tbest_transcript\tgene symbol\n};
	print SUM "$title2";
	my @head1  = split /\t/, $title1;
	$worksheet1->write_row( 0, 0, \@head1, $format{'title'});

	foreach my $x (keys %hash) {
		my $anno = exists $db{$x} ? $db{$x} : 'NA';
		my ($genome_len, $splice_len, $trans, $gene) = @{$hash{$x}};
		my ($chr, $start, $end, $strand, $type, $gene_id);

		if(exists $trans_type{$trans}){
			$type = $trans_type{$trans};
		}else{
			$type = "intergenic_region";
		}

		my @samples = ();
		my @counts  = ();
		foreach my $y (keys %{$all_circRNA{$x}}) {
			push @samples, $y;
			my @arr = split /\t/, $all_circRNA{$x}{$y};
			$chr    = $arr[1];
			$start  = $arr[2];
			$end    = $arr[3];
			$strand = $arr[10];
			push @counts, $arr[4];
		}
		my $sam = join ",", @samples;
		my $cnt = join ",", @counts;

		my @result = ($x, $anno, $chr, $start, $end, $strand, $genome_len, $splice_len, $sam, $cnt, $type, $trans, $gene);
		my $line   = join "\t", @result;
		my @res    = split /\t/, $line;
		$worksheet1->write_row($row1, 0, \@res, $format{'normal'});
		$row1++;

		my @result_1 = ($x, $anno, $chr, $start, $end, $strand, $genome_len, $splice_len, $type, $trans, $gene);
		my $line_1   = join "\t", @result_1;
		print SUM "$line_1\n";

	}

	$workbook->close();
	close SUM;

	print qq{circRNA 预测分析运行完成!\n};

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	# requir sample clean reads
	my $trim     = qq{$metadata->{project}/trim/trim_galore};
	foreach my $x (@samples) {
		die "$x R1 clean reads isn't exist!\n" if not -e qq{$trim/$x/$x\_R1_val_1.fq.gz};
		die "$x R2 clean reads isn't exist!\n" if not -e qq{$trim/$x/$x\_R2_val_2.fq.gz};
	}
	# reference genome bowtie1 index
	my $organ     = qq{$metadata->{organ}};
	my $ref_fasta = qq{$base->{$organ}{genome_fasta}};
	my $ref_gtf   = qq{$base->{$organ}{genome_mRNA_gtf}};
	my $bwa_index = qq{$base->{$organ}{bwa_index}};

	die "reference genome fasta file isn't exist!\n" if not -e qq{$ref_fasta};
	die "reference genome mRNA gtf file isn't exist!\n" if not -e qq{$ref_gtf};
	die "reference genome bwa index file isn't exist!\n" if not -e qq{$bwa_index.bwt};

}

sub res_check
{
	my $output = shift;
	my $sample = shift;
	my @temp_sample = ();
	foreach my $x (@{$sample}) {
		next if -e qq{$output/$x.finish};
		push @temp_sample, $x;
	}
	return @temp_sample;
}


sub cal_mean
{
        my $sum = 0;
        foreach my $x (@_) {
                $sum += $x;
        }
        my $length = scalar @_;
        return sprintf "%d", $sum / $length ;
}

sub parse_fasta
{
	my $fasta = shift;
	my %hash  = ();
	local $/  = ">";
	open FASTA, $fasta or die "Can't open $fasta!\n";
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($symbol, $seq) = split /\n/, $_, 2;
		my @arr    = split /\|/, $symbol;
		my $id     = qq{$arr[0]|$arr[1]};
		my $strand = $arr[2];
		my $trans  = $arr[3];
		my $gene   = $arr[4];
		$seq       =~ s/\s+//g;
		my $splice_len    = length($seq);
		my ($start, $end) = $id =~ /(\d+)\|(\d+)/;
		my $genome_len    = $end - $start + 1;
		$hash{$id} = [$genome_len, $splice_len, $trans, $gene]; 
		
	}
	close FASTA;
	return %hash;
}

sub parse_db
{
	my $txt = shift;

	my %hash = ();
	open TXT, $txt or die "Can't open $txt!\n";
	while (<TXT>) {
		chomp;
		my @arr = split /\t/;
		my $chr = $arr[0];
		my $start = $arr[1] + 1;
		my $end   = $arr[2];
		my $strand = $arr[3];
		my $val    = $arr[4];
		my $key = qq{$chr:$start|$end:$strand};

		$hash{$key} = $val;
	}
	close TXT;

	return %hash;
}

sub parse_type
{
	my $overlap = shift;
	my %trans   = ();
	open IN, $overlap or die "Can't open $overlap!\n";
	while(<IN>){
		my @arr = split/\t/,$_;
		if(exists $trans{$arr[9]}){
			if ($arr[10] eq "exon"){$trans{$arr[9]}++;}
		}else{
			if ($arr[10] eq "exon"){$trans{$arr[9]} = 1;}
			if ($arr[10] eq "intron"){$trans{$arr[9]} = 0;}
		}
	}
	close IN;
	foreach my $k(keys %trans){
		if($trans{$k} < 2){
			$trans{$k} = "intron";
		}else{
			$trans{$k} = "exon";
		}
	}
	return %trans;
}

1;
