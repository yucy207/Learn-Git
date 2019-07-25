my ($fasta, $out) = @ARGV;
fmt_fasta($fasta, $out);

sub fmt_fasta
{
	my $fasta = shift;
	my $out   = shift;
	my %hash  = ();
	local $/  = ">";
	open FASTA, $fasta or die "Can't open $fasta!\n";
	open SAVE, qq{>$out} or die "can't open $out!\n";
	while (<FASTA>) {
		chomp;
		next if /^\s*$/;
		my ($symbol, $seq) = split /\n/, $_, 2;
		my ($id) = $symbol =~ /^\d+\s+(\S+)/;
		print SAVE qq{>$id\n$seq};
	}
	close FASTA;
	close SAVE;
}