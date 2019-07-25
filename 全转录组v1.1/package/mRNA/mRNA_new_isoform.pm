package package::mRNA::mRNA_new_isoform;
use Parallel::ForkManager;
use package::format;
use Excel::Writer::XLSX;
use Encode qw/decode/;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my $organ    = qq{$metadata->{organ}};

	my $merge    = qq{$metadata->{project}/merge_transcript/result};
	my $result   = qq{$metadata->{project}/mRNA/novel_mRNA};
	my $report   = qq{$metadata->{report}/03_mRNA_Analysis/05_Novel_mRNA};
	my $map      = qq{$metadata->{project}/mapping/result};

	my $util         = qq{$base->{util}};
	my $transcoder   = qq{$base->{transdedcoder_bin}};
	my $gtf_to_fasta = qq{$base->{gtf_to_fasta_bin}};
	my $diamond      = qq{$base->{diamond_bin}};	
	my $cuffquant    = qq{$base->{cufflinks_bin}};
	my $cuffnorm     = qq{$base->{cufflinks_bin}};

	my $ref_fasta     = qq{$base->{$organ}{genome_fasta}};
	my $nr_db         = qq{$base->{nr_db}};
	my $gi_to_protein = qq{$base->{gi_to_protein}};

	if (not -d $merge) {

		print "[Error] : mRNA_new_isoform analysis is failed, you must assemble and merge transcript firstly!\n";
		exit;

	}

	if (-e qq{$report/New_Isoform_Expression_Summary.xlsx} ) {

		print qq{mRNA 新转录本分析已经完成!\n};
		#return;
	}

	pre_check($metadata, $base);

	system qq{mkdir -p $result/run}       if not -d qq{$result/run};
	system qq{mkdir -p $result/result}    if not -d qq{$result/result};
	system qq{mkdir -p $result/log}       if not -d qq{$result/log};
	system qq{mkdir -p $result/result/cuffquant} if not -d qq{$result/result/cuffquant};
	system qq{mkdir -p $result/result/cuffnorm}  if not -d qq{$result/result/cuffnorm};
	system qq{mkdir -p $report}                  if not -d $report;
	system qq{mkdir -p $report/New_Isoforms_Raw_Results} if not -d qq{$report/New_Isoforms_Raw_Results};

	########### filter exon num big than 1 transcript ###############
	predict_mRNA(qq{$merge/merged.gtf}, qq{$result/result/filter.gtf}) if not -e qq{$result/result/filter.gtf};
	die "all new isoform exon numbers small than 1!\n " if -s qq{$result/result/filter.gtf} == 0;

	open SAVE, qq{>$result/run/run.sh} or die "can't open $result/run/run.sh\n";
	my $fasta = qq{$gtf_to_fasta gtf_to_fasta $result/result/filter.gtf $ref_fasta $result/result/filter.fasta.tmp};
	my $cd    = qq{cd $result/result};
	my $fmt   = qq{perl $util/transcript_fasta_fmt.pl $result/result/filter.fasta.tmp $result/result/filter.fasta};
	my $trans_cmd = qq{$transcoder -t $result/result/filter.fasta};
	print SAVE qq{$fasta\n$cd\n$fmt\n$trans_cmd\n};
	close SAVE;

	system qq{bash $result/run/run.sh &> $result/log/run.log\n} if not -e qq{$result/result/filter.fasta.transdecoder_dir/longest_orfs.gff3};

	########### filter orf length big than 300 bp transcript #########
	filter_orf(qq{$result/result/filter.fasta}, qq{$result/result/filter.fasta.transdecoder_dir/longest_orfs.gff3},  qq{$result/result/filter_orf.fasta}) if not -e qq{$result/result/filter_orf.fasta};
	die "all new isoform orf length small than 300!\n " if -s qq{$result/result/filter_orf.fasta} == 0;

	########### filter reads can blast nr transcript #################
	system qq{rm -f /home/tmp/log} if -e qq{/home/tmp/log};
	my $blastp = qq{$diamond blastx -p 20 -d $nr_db -t /dev/shm -q $result/result/filter_orf.fasta -o $result/result/blast_nr.xls -e 10 &> /home/tmp/log};
	system qq{$blastp\n} if not -e qq{$result/result/blast_nr.xls};

	filter_blast(qq{$result/result/blast_nr.xls}, qq{$result/result/filter_orf.fasta}, qq{$result/result/filter.gtf}, qq{$result/result}) if not -e qq{$result/result/novel.mRNA.gtf};
	die "all new isoform blast nr identity small than  95%!\n " if -s qq{$result/result/novel.mRNA.gtf} == 0;

	######################### expression ############################
	my $max_threads = $base->{'max_threads'};
	my $pm = Parallel::ForkManager->new($max_threads);

	foreach my $x (@samples){

     	my $pid = $pm->start and next;

     	open SAVE, qq{>$result/run/$x.profile.sh} or die "Can't open $result/run/$x.profile.sh!\n";
     	my $cuffquant_cmd = qq{$cuffquant cuffquant -o $result/result/cuffquant/$x $result/result/novel.mRNA.gtf $map/$x/accepted_hits.bam\n};
	 	my $cuffnorm_cmd  = qq{$cuffnorm cuffnorm -o $result/result/cuffnorm/$x -L $x,$x $result/result/novel.mRNA.gtf $result/result/cuffquant/$x/abundances.cxb $result/result/cuffquant/$x/abundances.cxb\n};
	 	print SAVE qq{$cuffquant_cmd\n$cuffnorm_cmd\n};
	 	close SAVE;
	 	system qq{bash $result/run/$x.profile.sh &> $result/log/$x.profile.log} if not -e qq{$result/result/cuffquant/$x/abundances.cxb};
	 	
	 	$pm->finish;

	 }

	$pm->wait_all_children;

	merge_profile(\@samples, qq{$result/result/cuffnorm}, qq{$result/result/gene_expression.xls}) if not -e qq{$result/result/gene_expression.xls};
	filter_trans(qq{$result/result/gene_expression.xls}, qq{$result/result/high_expression.xls} ) if not -e qq{$result/result/high_expression.xls};
	select_fasta(qq{$result/result/high_expression.xls}, qq{$result/result/novel.mRNA.fasta}, qq{$report/New_Isoforms_Raw_Results/new.isoform.fasta}) if not -e qq{$report/New_Isoforms_Raw_Results/new.isoform.fasta};
	select_gtf(qq{$result/result/high_expression.xls}, qq{$result//result/novel.mRNA.gtf}, qq{$report/New_Isoforms_Raw_Results/new.isoform.gtf}) if not -e qq{$report/New_Isoforms_Raw_Results/new.isoform.gtf};
	assign_gi(qq{$result/result/novel.mRNA.fasta}, qq{$result/result/blast_nr.xls}, qq{$result/result/high_expression.xls}, qq{$gi_to_protein}, qq{$report/New_Isoform_Expression_Summary.xlsx});
	
	print qq{mRNA 新转录本分析已经完成!\n};

}


sub pre_check
{
	my $metadata  = shift;
	my $base      = shift;

	my $merge     = qq{$metadata->{project}/merge_transcript/result};
	my $organ     = qq{$metadata->{organ}};
	my $ref_fasta = qq{$base->{$organ}{genome_fasta}};
	my $nr_db     = qq{$base->{nr_db}};

	die "merged gtf file isn't exist!\n"   if not -e $merge;
	die "genome fasta file isn't exist!\n" if not -e $ref_fasta;
	die "genome diamond database isn't exist!\n" if not -e qq{$nr_db.dmnd};

}

sub predict_mRNA
{
	my $gtf  = shift;
	my $out  = shift;
	my %meta = ();
	open GTF, $gtf or die "Can't open $gtf\n";
	while (<GTF>) {
		my @arr = split /\t/;
		next if not /class_code\s+\"u\";/;
		my ($trans) = $_ =~ /transcript_id\s+\"(.+?)\";/;
		if (exists $meta{$trans}) {
			$meta{$trans}->[0] .= $_;
			$meta{$trans}->[2] = $arr[4];
		} else {
			$meta{$trans} = [$_, $arr[3], $arr[4]];
		} 
		
	}
	close GTF;

	open SAVE, qq{>$out} or die "Can't open $out\n";
	foreach my $x (keys %meta) {
		my $len = $meta{$x}->[2] - $meta{$x}->[1];
		my @lines = split /\n/, $meta{$x}->[0];
		my $exon_num = scalar @lines;
		next if $exon_num <= 1;
		print SAVE qq{$meta{$x}->[0]};
	}
	close SAVE;
}

sub filter_orf
{
	my $fasta = shift;
	my $orf   = shift;
	my $out   = shift;
 	my %trans = ();
	open ORF, $orf or die "Can't open $orf\n";
	while (<ORF>) {
		chomp;
		next if not /CDS/;
		my @arr  = split /\t/;
		my $len  = $arr[4] - $arr[3];
		$trans{$arr[0]} = $len if $len > 300;
		#$trans{$arr[0]} = $len;
	}
	close ORF;

	local $/ = ">";
	my %meta = ();
	open FASTA, $fasta or die "Can't open $fasta\n";
	
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($id, $seq) = split /\n/, $_, 2;

		next if not exists $trans{$id};
		$meta{$id} = $seq;

	}
	close FASTA;
	local $/ = "\n";
	open SAVE, qq{>$out} or die "Can't open $out\n";
	foreach my $x (keys %meta) {
		print SAVE qq{>$x\n$meta{$x}};
	}
	close SAVE;

}

sub filter_blast
{
	my $blast   = shift;
	my $fasta   = shift;
	my $gtf     = shift;
	my $out_dir = shift;

	my %trans  = ();
	open BLAST, $blast or die "Can't open $blast\n";
	while (<BLAST>) {
		chomp;
		my @arr = split /\t/;
		next if $arr[2] <= 95;
		$trans{$arr[0]} = "-";
	}
	close BLAST;

	local $/ = ">";
	my %meta = ();
	my %trans_name = ();
	open FASTA, $fasta or die "Can't open $fasta\n";
	
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($id, $seq) = split /\n/, $_, 2;


		next if not exists $trans{$id};
		$trans_name{$id} = "-";
		$meta{$id} = $seq;

	}
	close FASTA;

	local $/ = "\n";
	open SAVE, qq{>$out_dir/novel.mRNA.fasta} or die "Can't open $out_dir/novel.mRNA.fasta!\n";
	foreach my $x (keys %meta) {
		print SAVE qq{>$x\n$meta{$x}};
	}
	close SAVE;


	open GTF, $gtf or die "Can't open $gtf!\n";
	open SAVE, qq{>$out_dir/novel.mRNA.gtf} or die "Can't open $out_dir/novel.mRNA.gtf!\n";
	while (<GTF>) {
		chomp;
		my ($trans) = $_ =~ /transcript_id\s+\"(.+?)\";/;
		next if not exists $trans_name{$trans};
		print SAVE qq{$_\n};
	}
	close SAVE;
	close GTF;

}

sub merge_profile
{
	my $samples = shift;
	my $dir     = shift;
	my $out     = shift;

	my %meta    = ();
	foreach my $x (@{$samples}) {
		my $count = qq{$dir/$x/isoforms.fpkm_table\n};
		open EXO, $count or die "Can't open $count\n";
		while (<EXO>) {
			chomp;
			next if /^tracking_id/;
			my @arr = split /\s+/;
			$meta{$arr[0]}{$x} = $arr[1];
		}
		close EXO;
	}

	open SAVE, qq{>$out} or die "Can't open $out\n";
	my $title = join "\t", @{$samples};
	print SAVE qq{gene\t$title\n};
	foreach my $key (sort keys %meta) {
		my $res = join "\t", map { $meta{$key}{$_} } @{$samples};
		print SAVE qq{$key\t$res\n};
	}

	close SAVE;
}

sub filter_trans
{
	my $count = shift;
	my $out   = shift;
	open EXO, $count or die "Can't open $count!\n";
	open SAVE, qq{>$out} or die "Can't open $out!\n";
	while (<EXO>) {
		chomp;
		if (/^gene/) {
			print SAVE qq{$_\n};
		} else {
			my @arr = split /\t/;
			my $sum = 0;
			map {$sum += $arr[$_]} 1..($#arr);
			my $mean = $sum / $#arr;
			next if $mean < 10;
			print SAVE qq{$_\n};
		}
	}
	close EXO;
	close SAVE;
}

sub select_gtf
{
	my $count = shift;
	my $gtf   = shift;
	my $out   = shift;

	my %ids   = ();
	open EXO, $count or die "Can't open $count!\n";
	while (<EXO>) {
		chomp;
		next if /^gene/;
		my @arr = split /\t/;
		$ids{$arr[0]} = "-";
	}
	close EXO;


	open GTF, $gtf or die "Can't open $gtf!\n";
	open SAVE, qq{>$out} or die "Can't open $out!\n";
	while (<GTF>) {
		chomp;
		my ($symbol) = $_ =~ /(TCONS\_\d+)/;
		next if not exists $ids{$symbol};
		print SAVE qq{$_\n};

	}
	close GTF;
	close SAVE;
}

sub select_fasta
{
	my $count = shift;
	my $fasta = shift;
	my $out   = shift;

	my %ids   = ();
	open EXO, $count or die "Can't open $count!\n";
	while (<EXO>) {
		chomp;
		next if /^gene/;
		my @arr = split /\t/;
		$ids{$arr[0]} = "-";
	}
	close EXO;

	local $/ = ">";
	open FASTA, $fasta or die "Can't open $fasta!\n";
	open SAVE, qq{>$out} or die "Can't open $out!\n";
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($id, $seq) = split /\n/, $_, 2;
		my ($symbol) = $id =~ /(TCONS\_\d+)/;
		next if not exists $ids{$symbol};
		print SAVE qq{>$symbol\n$seq};

	}
	close FASTA;
	close SAVE;
}

sub assign_gi
{
	my $fasta = shift;
	my $blast = shift;
	my $count = shift;
	my $gi    = shift;
	my $out   = shift;


	my %trans2id = ();
	local $/ = ">";
	open FASTA, $fasta or die "Can't open $fasta!\n";
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($id, $seq) = split /\n/, $_, 2;
		my ($acc) = $id =~ /(\d+)/;
		my ($trans) = $id =~ /(TCONS_\d+)/;
		$trans2id{$trans} = $trans;
	}
	close FASTA;
	$/ = "\n";


	my %trans2gi = ();
	open EXO, $blast or die "Can't open $blast!\n";
	while (<EXO>) {
		chomp;
		my @arr = split /\t/;
		next if $arr[2] < 80;
		next if $arr[10] > 1e-5;
		my $gi = $arr[1];
		$trans2gi{$arr[0]} = $gi if not exists $trans2gi{$arr[0]};
	}
	close EXO;

	my %meta  = ();
	open EXO, $gi or die "Can't open $gi!\n";
	while (<EXO>) {
		chomp;
		my ($gi, $id) = split /\t/;
		next if  $gi =~ /\s+/;
		$meta{$gi} = $id;
	}
	close EXO;

	open COUNT, $count or die "Can't open $count!\n";
	my $excel     = qq{$out};
	my $workbook  = Excel::Writer::XLSX->new($excel);
	my %format    = package::format::run($workbook);

	my $worksheet = $workbook->add_worksheet("New Isoform");
	my $row = 0;
	while (<COUNT>) {
		chomp;
		my @arr      = split /\t/;
		my @counts   = @arr[1..$#arr];
		my $title    = join "\t", @counts;
		if (/^gene/) {
			my $line =  qq{gene\tprotein\t$title};
			my @head = split /\t/, $line;
			$worksheet->write_row($row, 0, \@head, $format{'title'});
		} else { 
			my $id   = qq{$trans2id{$arr[0]}};	
			if (exists $trans2gi{$id}) {
				my $gi      = qq{$trans2gi{$id}};
				my $anno    = (split /\|/, $gi)[0];
				my $protein = $meta{$anno};
				my $line = qq{$arr[0]\t$anno|$protein\t$title};
				my @res = split /\t/, $line;
				$worksheet->write_row($row, 0, \@res, $format{'normal'});

			} else {
				my @res = @counts;
				unshift @res, "";
				unshift @res, $arr[0];
				$worksheet->write_row($row, 0, \@res, $format{'normal'});
			}
		}
		$row++;
	}
	close COUNT;
	$workbook->close();

}

1;
