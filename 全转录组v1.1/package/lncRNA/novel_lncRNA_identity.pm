package package::lncRNA::novel_lncRNA_identity;
use Parallel::ForkManager;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my $organ    = qq{$metadata->{organ}};

	my $merge    = qq{$metadata->{project}/merge_transcript/result};
	my $report   = qq{$metadata->{report}/04_lncRNA_Analysis/05_Novel_lncRNA};
	my $result   = qq{$metadata->{project}/lncRNA/novel_lncRNA};

	my $util           = qq{$base->{util}};
	my $cuffcompare    = qq{$base->{cufflinks_bin}};
	my $transcoder     = qq{$base->{transdedcoder_bin}};
	my $gtf_to_fasta   = qq{$base->{gtf_to_fasta_bin}};
	my $diamond        = qq{$base->{diamond_bin}};

	my $ref_lncRNA_gtf = qq{$base->{$organ}{genome_lncRNA_gtf}};
	my $ref_fasta      = qq{$base->{$organ}{genome_fasta}};
	my $nr_db          = qq{$base->{nr_db}};

	if (-e qq{$report/New_Isoforms_Raw_Results/new.isoform.fasta} ) {
		print qq{lncRNA 新转录本分析已经完成!\n};
		return;
	}

	pre_check($metadata, $base);

	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $report/New_Isoforms_Raw_Results}  if not -d qq{$report/New_Isoforms_Raw_Results};

	# filter exon 1 and length samll than 200bp transcript
	predict_lncRNA(qq{$merge/merged.gtf}, qq{$result/result/filter.gtf}) if not -e qq{$result/result/filter.gtf};
	die "all new isoform exon numbers small than 1!\n " if -s qq{$result/result/filter.gtf} == 0;

	open SAVE, qq{>$result/run/run.sh} or die "can't open $result/run/run.sh\n";
	my $cmd       = qq{$cuffcompare cuffcompare -r $ref_lncRNA_gtf $merge/merged.gtf -o $result/result/compare};
	my $grep      = qq{grep "class_code \\"u\\"" $result/result/compare.combined.gtf > $result/result/compared.gtf\n};
	my $fasta     = qq{$gtf_to_fasta gtf_to_fasta $result/result/compared.gtf $ref_fasta $result/result/compared.fasta.tmp};
	my $cd        = qq{cd $result/result};
	my $fmt       = qq{perl $util/transcript_fasta_fmt.pl $result/result/compared.fasta.tmp $result/result/compared.fasta};
	my $trans_cmd = qq{$transcoder -t $result/result/compared.fasta};
	print SAVE qq{$cmd\n$grep\n$fasta\n$cd\n$fmt\n$trans_cmd\n};
	close SAVE;
	system qq{bash $result/run/run.sh &> $result/log/run.log\n} if not -e qq{$result/result/filter.fasta};

	filter_orf(qq{$result/result/compared.fasta}, qq{$result/result/compared.fasta.transdecoder_dir/longest_orfs.gff3},  qq{$result/result/filter_orf.fasta}) if not -e qq{$result/result/filter_orf.fasta};
	die "all new isoform orf length big than 300!\n " if -s qq{$result/result/filter_orf.fasta} == 0;

	my $blastp = qq{$diamond blastx -p 20 -d $nr_db -t /dev/shm -q $result/result/filter_orf.fasta -o $result/result/blast_nr.xls -e 1e-5};
	system qq{$blastp &> $result/result/blastp.log\n} if not -e qq{$result/result/blast_nr.xls};

	filter_blast(qq{$result/result/blast_nr.xls}, qq{$result/result/filter_orf.fasta}, qq{$merge/merged.gtf}, qq{$report/New_Isoforms_Raw_Results});

}

sub pre_check
{
	my $metadata = shift;
	my $base     = shift;
	my $organ    =  qq{$metadata->{organ}};

	my $merge    = qq{$metadata->{project}/merge_transcript/result};
	die "merged gtf file isn't exist!\n" if not -e $merge;

	my $ref_lncRNA_gtf =  qq{$base->{$organ}{genome_lncRNA_gtf}};
	die "genome lncRNA gtf file isn't exist!\n" if not -e $ref_lncRNA_gtf;

}

sub predict_lncRNA
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
		next if $len < 100;
		my @lines = split /\n/, $meta{$x}->[0];
		my $exon_num = scalar @lines;
		#next if $exon_num <= 1;
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
		#$trans{$arr[0]} = $len if $len < 300;
		$trans{$arr[0]} = $len;
	}
	close ORF;

	local $/ = ">";
	my %meta = ();
	open FASTA, $fasta or die "Can't open $fasta\n";
	
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($id, $seq) = split /\n/, $_, 2;
		next if exists $trans{$id};
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
		#next if $arr[2] >= 95;
		$trans{$arr[0]} = "-";
	}
	close BLAST;

	local $/ = ">";
	my %meta = ();
	open FASTA, $fasta or die "Can't open $fasta\n";
	
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($id, $seq) = split /\n/, $_, 2;
		next if exists $trans{$id};
		$meta{$id} = $seq;

	}
	close FASTA;

	local $/ = "\n";
	open SAVE, qq{>$out_dir/new.isoform.fasta} or die "Can't open $out_dir/new.isoform.fasta\n";
	foreach my $x (keys %meta) {
		print SAVE qq{>$x\n$meta{$x}};
	}
	close SAVE;

	open GTF, $gtf or die "Can't open $gtf!\n";
	open SAVE, qq{>$out_dir/new.isoform.gtf} or die "Can't open $out_dir/novel.mRNA.gtf!\n";
	while (<GTF>) {
		chomp;
		my ($trans) = $_ =~ /transcript_id\s+\"(.+?)\";/;
		next if not exists $meta{$trans};
		print SAVE qq{$_\n};
	}
	close SAVE;
	close GTF;

}


1;

