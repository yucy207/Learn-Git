#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my ($gtf, $bed, $out_dir, $fasta,  $help);

GetOptions(
	'gtf|g=s'     => \$gtf,
	'bed|b=s'     => \$bed,
	'fasta|f=s'   => \$fasta,
	'out_dir|o=s' => \$out_dir,
	'help|h!'     => \$help
);

if ($help or not $gtf or not $bed or not $out_dir or not $fasta) {
	usage();
	exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl  -gtf refGene_gtf -bed circRNA_bed -fasta genome.fasta -out_dir out_dir
	-gtf     -g  refgene gtf file [required]
	-bed     -b  circRNA bed fiile [required]
	-fasta   -f  genome fasta file [required]
	-out_dir -o  output directory [required]
	-help    -h  print help message
EOF

print $help;

}

make_dir();
parse_gtf();
cal_overlap();
fmt_overlap();


my %circrnas = ();
my %intragenic_circrnas = ();
my %trans2gene  = ();

parse_circrna_bed();
parse_overlap();
intergenic_circrna();
bed2fasta();

sub make_dir
{
	system qq{mkdir -p $out_dir} if not -d $out_dir;

}

sub parse_gtf
{

	my $trans_id = '';
	my $end_pos  = 0;
	open BED, qq{>$out_dir/refgene.bed} or die "Can't open $out_dir/refgene.bed!\n";
	open GTF, $gtf or die "Can't open $gtf!\n";
	while (<GTF>) {
		chomp;
		my @arr = split /\t/;
		next if $arr[2] ne 'exon';
		my $chr     = $arr[0];
		my $start   = $arr[3];
		my $end     = $arr[4];
		my $strand  = $arr[6];
		my ($trans) = $_ =~ /transcript_id\s+\"(\S+)\"/;
		my ($gene)  = $_ =~ /gene_id\s+\"(\S+)\"/;
		my ($exon_num)  = $_ =~ /exon_number\s+\"(\d+)\"/;
		my $name    = qq{$trans|exon_$exon_num|$gene};
		if ($trans ne $trans_id) {
			print BED qq{$chr\t$start\t$end\t$name\t.\t$strand\n};
			$trans_id = $trans;
			$end_pos  = $end;
		} else {
			if ($end_pos < $start){
				my $intron_num = $exon_num  - 1;
				print BED qq{$chr\t$end_pos\t$start\t$trans|intron_$intron_num|$gene\t.\t$strand\n};
			}
			print BED qq{$chr\t$start\t$end\t$name\t.\t$strand\n};
			$trans_id = $trans;
			$end_pos  = $end;
		}

	}
	close GTF;
	close BED;
}


sub cal_overlap
{
	system qq{/home/xudl/soft/bedtools2-2.25.0/bin/bedtools intersect  -a $bed -b $out_dir/refgene.bed  -wb -s > $out_dir/overlap.xls};
}

sub parse_circrna_bed
{
	open BED, $bed or die "Can't open $bed!\n";
	while (<BED>) {
		chomp;
		my @arr  = split /\t/;
		$circrnas{$arr[3]} = $_;
	}
	close BED;
}

sub fmt_overlap
{
	open SAVE, qq{>$out_dir/fmt.overlap.xls} or die "can't open $out_dir/fmt.overlap.xls!\n";
	open TXT,  qq{$out_dir/overlap.xls} or die "Can't open $out_dir/overlap.xls!\n";
	while (<TXT>) {
		chomp;
		my @arr = split /\t/;
		my ($trans, $pos, $gene)  = split /\|/, $arr[9];
		my ($type, $num) = split /\_/, $pos;
		my $line  = join "\t", @arr[0..8];
		print SAVE qq{$line\t$trans\t$type\t$num\t$gene\t$arr[10]\t$arr[11]\n};
	}
	close TXT;
	close SAVE;

	system qq{sort -k4,4 -k10,10 $out_dir/fmt.overlap.xls > $out_dir/overlap.sort.xls};
}



sub parse_overlap
{

	open TXT, qq{$out_dir/overlap.sort.xls} or die "Can't open $out_dir/overlap.sort.xls!\n";
	while (<TXT>) {
		chomp;
		my @arr = split /\t/;
		my $circRNA_id = $arr[3];
		my $trans      = $arr[9];
		my $type       = $arr[10];
		my $gene       = $arr[12];
		my $start      = $arr[1];
		my $end        = $arr[2];

		$trans2gene{$trans} = $gene;
		push @{$intragenic_circrnas{$circRNA_id}{$trans}}, [$type, $start, $end];

	}
	close TXT;

	open SAVE, qq{>$out_dir/intragenic_circrnas.bed} or die "Can't open $out_dir/intragenic_circrnas.bed!\n";
	open BEST, qq{>$out_dir/best_intragenic_circrnas.xls} or die "Can't open $out_dir/best_intragenic_circrnas.xls!\n";
	foreach my $x (keys %intragenic_circrnas) {
		my ($chr, $circRNA_start, $circRNA_end, $circRNA_id, $na, $strand) = split /\t/, $circrnas{$x};


		my %hash = ();
		my %target_bed = ();
		foreach my $y (sort {$a cmp $b} keys %{$intragenic_circrnas{$x}}) {
			my $splice_len = 0;
			my $cnt  = scalar @{$intragenic_circrnas{$x}{$y}};

			# determine a transcript have exons, if exon number > 0, flag_type = 1;
			my $flag_type = 0;
			# determine a transcript exon bundary is exact mathch circrna position


			my $start_exon_type = 0;
			my $end_exon_type   = 0;

			my ($overlap_start, $overlap_end) = (0, 0);
			foreach my $z (1..$cnt) {
				my $index = $z - 1;
				my ($type, $start, $end) = @{$intragenic_circrnas{$x}{$y}[$index]};
				$flag_type = 1 if $type eq 'exon';
				$overlap_start = $start if $z == 1;
				$overlap_end   = $end   if $z == $cnt;
				$start_exon_type  = 1 if $type eq "exon" and  $z == 1    and $start == $circRNA_start;
				$end_exon_type    = 1 if $type eq "exon" and  $z == $cnt and $end   == $circRNA_end;

				next if $type eq 'intron' and $z > 1 and $z < $cnt;
				my $len = ($end - $start + 1);
				$splice_len += $len;
				push @{$target_bed{$y}}, [$start, $end];
			}

			$flag_type = 2 if $start_exon_type == 1 and $end_exon_type == 1;

			if ( $overlap_start != $circRNA_start ) {
				$splice_len += ($overlap_start - $circRNA_start + 1);
				push @{$target_bed{$y}}, [$circRNA_start, $overlap_start];
			}

			if ( $overlap_end   != $circRNA_end ) {
				$splice_len += ($circRNA_end   - $overlap_end   + 1);
				push @{$target_bed{$y}}, [$overlap_end, $circRNA_end];
			}


			print  BEST qq{$x\t$y\t$splice_len\n};
			$hash{$flag_type}{$y} = $splice_len;
		}
		# best transcript id
		my  $opt_trans;
		my  $splice_len;

		foreach my $x (sort keys %hash) {
			next if not exists $hash{$x};
			my  $max_len = 0;
			foreach my $y (sort {$a cmp $b} keys %{$hash{$x}}) {

				if ($hash{$x}{$y} >= $max_len) {
					$opt_trans = $y;
					$max_len   = $hash{$x}{$y};
					$splice_len = $max_len;
				}
			}	
		}



		print BEST qq{BEST:$x\t$opt_trans\t$splice_len\n};

		foreach my $x (@{$target_bed{$opt_trans}}) {
			my ($start, $end) = @{$x};
			my $gene  = $trans2gene{$opt_trans};
			print SAVE qq{$chr\t$start\t$end\t$circRNA_id\t.\t$strand\t$opt_trans\t$gene\n};
		}


	}
	close SAVE;
	close BEST;
}



sub intergenic_circrna
{
	open SAVE, qq{>$out_dir/intergenic_circrnas.bed} or die "Can't open $out_dir/intergenic_circrna.bed!\n";
	foreach my $x (keys %circrnas) {
		next if exists $intragenic_circrnas{$x};
		print SAVE qq{$circrnas{$x}\tNONE\tNONE\n};

	}
	close SAVE;
}

sub reverse_complement
{
	my $seq = shift;

	$seq =~ tr/ATGCN/TACGN/;
	return reverse $seq;

}


sub extract_fasta
{
	my $bed   = shift;
	my $fasta = shift;
	my $out   = shift;



	my %hash  = ();
	local $/  = ">";
	open FASTA, $fasta or die "Can't open $fasta!\n";
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($symbol, $seq) = split /\n/, $_, 2;
		$seq =~ s/\s+//g;
		$hash{$symbol} = $seq;
		next if $symbol !~ /\s+/;
		my ($id, $summary) = split /\s+/, $symbol, 2;
		$hash{$id} = $seq;
	}
	close FASTA;
	$/ = "\n";

	my %circrna_info  = ();
	my %circrna_region = ();
	open BED, $bed or die "Can't open $bed!\n";
	while (<BED>) {
		chomp;
		my @arr = split /\t/;
		$circrna_info{$arr[3]} = qq{$arr[3]|$arr[5]|$arr[6]|$arr[7]};
		push @{$circrna_region{$arr[3]}}, [$arr[0], $arr[1], $arr[2]];
	}
	close BED;

	open SAVE, qq{>$out} or die "Can't open $out!\n";
	foreach my $x (keys %circrna_region) {
		my $id = $circrna_info{$x};
		my @arr = split /\|/, $id;
		my $strand = $arr[$#arr - 2];
		my $seq;
		foreach my $y (@{$circrna_region{$x}}) {
			my ($chr, $start, $end) = @{$y};
			my $part = substr($hash{$chr}, $start - 1, $end - $start +1);
			$seq .= $part;
		}
		$seq = uc($seq);
		$seq = reverse_complement($seq) if $strand eq '-';
		print SAVE qq{>$id\n$seq\n};
	}
	close SAVE;

}

sub bed2fasta
{
	system qq{cat $out_dir/intragenic_circrnas.bed $out_dir/intergenic_circrnas.bed > $out_dir/circRNA.trans.bed\n};
	extract_fasta(qq{$out_dir/circRNA.trans.bed},  $fasta, qq{$out_dir/circrnas_splice_transcript.fasta});

}

