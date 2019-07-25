#!/usr/bin/env perl
use Parallel::ForkManager;
use List::Util qw /min max/;

my ($diff, $gtf_to_fasta, $lncRNA_gtf, $genome_fasta, $mRNA_gtf, $out_dir) = @ARGV;

my %lncRNA_gene  = ();
my %lncRNA_trans = ();

my %mRNA_locus   = ();
my %target_mRNA  = ();

my %lncRNA_locus = ();


system qq{mkdir -p $out_dir} if not -d $out_dir;

extract_diff_lncRNA_gene();
extract_diff_lncRNA_transcript();
extract_lncRNA_location();
extract_neighbor_mRNA_gene_id();
run();

sub extract_diff_lncRNA_gene
{

	open EXO, $diff or die "Can't open $diff!\n";
	while (<EXO>) {
		chomp;
		$lncRNA_gene{$_} = "-";
	}
	close EXO;

}

sub extract_diff_lncRNA_transcript
{
	open GTF, $lncRNA_gtf or die "Can't open $lncRNA_gtf!\n";
	while (<GTF>) {
		chomp;
		my ($trans) = $_ =~ /transcript_id\s+\"(.+?)\"/;
		my ($gene) = $_ =~ /gene_id\s+\"(.+?)\"/;
		next if not exists $lncRNA_gene{$gene};
		next if exists $lncRNA_trans{$gene};
		#print qq{$trans\t$gene\n};
		push @{$lncRNA_trans{$gene}}, $trans if not $trans ~~ @{$lncRNA_trans{$gene}} ; 
	}
	close GTF;
}

sub extract_lncRNA_location
{
	my %hash = ();
	my %meta = ();
	open GTF, $lncRNA_gtf;
	while (<GTF>) {

        my ($trans) = $_ =~ /transcript_id\s+\"(.+?)\"/;
        my ($gene) = $_ =~ /gene_id\s+\"(.+?)\"/;
        my @arr = split /\t/;
        $meta{$gene} = $arr[0];

        next if $arr[2] ne 'exon';
        push @{$hash{$gene}}, $arr[3];
        push @{$hash{$gene}}, $arr[4];

	}

	close GTF;

	foreach my $x (keys %hash) {
		my $chr   = $meta{$x};
		my $start = min @{$hash{$x}};
		my $end   = max @{$hash{$x}};
		#print qq{$x\t$chr:$start-$end\n};
		$lncRNA_locus{$x} = [$chr, $start, $end];
	}

}

sub extract_neighbor_mRNA_gene_id
{
	open GTF, $mRNA_gtf or die "Can't open $mRNA_gtf!\n";
	while (<GTF>) {
		chomp;
		my @arr    = split /\t/;
		my ($gene) = $_ =~ /gene_id\s+\"(.+?)\"/;
		next if not defined $gene;
		foreach my $x (keys %lncRNA_gene) {
			my $chr   = $lncRNA_locus{$x}->[0];
			my $start = $lncRNA_locus{$x}->[1];
			my $end   = $lncRNA_locus{$x}->[2];

			my $step  = 10000;
			next if $arr[0] ne $chr;
			next if $arr[3] > $end + $step ;
			next if $arr[4] < $start - $step;

			map { push@{$target_mRNA{$_}}, $gene if not $gene ~~ @{$target_mRNA{$_}}  }  @{$lncRNA_trans{$x}};

		}

	}
	close GTF;

}

sub run
{
	my $max_threads = 30;
	my $pm = Parallel::ForkManager->new($max_threads);

	my @val = keys %target_mRNA;
	my $cnt = @val;
	print qq{$cnt\n};
	foreach my $trans (keys %target_mRNA) {
		next if -e qq{$out_dir/$trans/$trans.finish};
		my $pid = $pm->start and next;
		my $temp_out = qq{$out_dir/$trans};
		my $mRNA_ids = $target_mRNA{$trans};
		# 获取随机长度字符
		my $random_word = get_random_word(10);
		# extract lncRNA fasta
		extract_diff_lncRNA_transcript_fasta($lncRNA_gtf, $genome_fasta, $trans, $random_word, $temp_out);
		# extract mRNA fasta
		extract_neighbor_mRNA_gene_fasta($mRNA_gtf, $genome_fasta, $mRNA_ids, $trans, $random_word, $temp_out);
		# run lncTar
		system qq{perl /home/xudl/soft/LncTar/LncTar.pl -p 1 -l $temp_out/$trans.$random_word.input -m $temp_out/$trans.$random_word.mRNA.input  -d -0.1 -s F -o $temp_out/$trans.target.txt &> $temp_out/$trans.log\n};
		system qq{rm -f /home/tmp/$trans.$random_word.mRNA.input_mRNA_temp.txt};
		system qq{rm -f /home/tmp/$trans.$random_word.input_lncRNA_temp.txt};
		system qq{touch $out_dir/$trans/$trans.finish};

		$pm->finish;
	}

	$pm->wait_all_children;

}

sub get_random_word
{
	my $length   = shift @_;
	my @aword    = (0..9,'a'..'z','A'..'Z');
	my $password = join '', map { $aword[int rand @aword] } 0..($length-1);
   	return $password;
}


sub extract_diff_lncRNA_transcript_fasta
{
	my $lncRNA_gtf   = shift;
	my $genome_fasta = shift;
	my $trans_ids    = shift;
	my $random_word  = shift;
	my $out_dir      = shift;

	my %meta         = ();

	system qq{mkdir -p $out_dir} if not -d $out_dir;
	open GTF, $lncRNA_gtf or die "Can't open $lncRNA_gtf!\n";
	open SAVE, qq{>$out_dir/lncRNA.gtf} or die "Can't open $out_dir/lncRNA.gtf!\n";
	while (<GTF>) {
		chomp;
		my ($trans) = $_ =~ /transcript_id\s+\"(.+?)\"/;
		my ($gene)  = $_ =~ /gene_id\s+\"(.+?)\"/;
		next if $trans ne $trans_ids;
		print SAVE qq{$_\n};
		$meta{$trans} = $gene;
	}
	close GTF;
	close SAVE;
	system qq{$gtf_to_fasta $out_dir/lncRNA.gtf $genome_fasta $out_dir/lncRNA.fasta};
	parse_fasta(qq{$out_dir/lncRNA.fasta}, qq{$out_dir/$trans_ids.$random_word.input}, \%meta);

}


sub extract_neighbor_mRNA_gene_fasta
{
	my $mRNA_gtf     = shift;
	my $genome_fasta = shift;
	my $gene_ids     = shift;
	my $trans_ids    = shift;
	my $random_word  = shift;
	my $out_dir      = shift;

	my %meta         = ();

	system qq{mkdir -p $out_dir} if not -d $out_dir;
	open GTF, $mRNA_gtf or die "Can't open $mRNA_gtf!\n";
	open SAVE, qq{>$out_dir/$trans_ids.mRNA.gtf} or die "Can't open $out_dir/$trans_ids.mRNA.gtf!\n";
	while (<GTF>) {
		chomp;
		my @arr     = split /\t/;
		my ($trans) = $_ =~ /transcript_id\s+\"(.+?)\"/;
		my ($gene)  = $_ =~ /gene_id\s+\"(.+?)\";/;
		next if not $gene ~~ @{$gene_ids};
		print SAVE qq{$_\n};
		$meta{$trans} = $gene;
		
	}
	close GTF;
	close SAVE;
	system qq{$gtf_to_fasta $out_dir/$trans_ids.mRNA.gtf $genome_fasta $out_dir/$trans_ids.mRNA.fasta};

	parse_fasta(qq{$out_dir/$trans_ids.mRNA.fasta}, qq{$out_dir/$trans_ids.$random_word.mRNA.input}, \%meta);
}



sub parse_fasta
{
	my $fasta  = shift;
	my $out    = shift;
	my $hash   = shift;

	local $/ = ">";
	open FASTA, $fasta or die "Can't open $fasta!\n";
	open SAVE, qq{>$out} or die "Can't open $out!\n";
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($id, $seq) = split /\n/, $_ , 2;
		my ($symbol) = $id =~ /\d+\s+(.+?)\s+/;
		$seq =~ s/\s+//g;
		my $gene = $hash->{$symbol};
		print SAVE qq{>$symbol($gene)\n$seq\n};
	}
	close FASTA;
	close SAVE;
}