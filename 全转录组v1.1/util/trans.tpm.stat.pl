#!usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

my ($gtflist, $tpmfile) = @ARGV;
$gtflist = abs_path($gtflist);
$tpmfile = abs_path($tpmfile);

my %sample_hash = ();
my %gene_hash   = ();
my @samples;

open(GTFLIST, $gtflist) or "die cannot open gtflist!\n";
while(<GTFLIST>){
	chomp;
	my ($sample, $gtf_dir) = split /\t/, $_;
	open(GTF, $gtf_dir) or "die cannot open gtf_dir!\n";
	while(<GTF>){
		chomp;
		next if /exon/ or /^#/;
		my ($trans_id) = $_ =~/transcript_id\s+\"(.+?)\"\;/;
		my ($gene_id)  = $_ =~/gene_id\s+\"(.+?)\"\;/;
		my ($tpm)      = $_ =~/TPM\s+\"(.+?)\"\;/;
		if(exists $sample_hash{$trans_id}{$sample}){
			if($tpm >= $sample_hash{$trans_id}{$sample}){

				$sample_hash{$trans_id}{$sample} = $tpm;
				$gene_hash{$trans_id} = $gene_id;
			}
		}else{

			$sample_hash{$trans_id}{$sample} = $tpm;
			$gene_hash{$trans_id} = $gene_id;

		}		
	}
	close GTF;
	push @samples, $sample;
}
close GTFLIST;

open(OUT, ">$tpmfile") or "die cannot make tpmfile!\n";
my $head = join("\t", @samples);
print OUT "transcript_ID\tgene_ID\t$head\n";
foreach my $x(keys %sample_hash){
	my @res;
	my $result;
	foreach my $y(@samples){
		push @res, $sample_hash{$x}{$y};
		$result = join("\t", @res);		
	}
	print OUT "$x\t$gene_hash{$x}\t$result\n";
}
