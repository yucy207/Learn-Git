sub read_bed_chr{
    my $bed_file = shift @_;
    my %hashBED;
    open BED, $bed_file;
    my $count = 0;
    while(<BED>)
    {   
        $count++;
        my ($chr, $tmp) = split /\t/, $_;
        $hashBED{$chr}{'Chr_Sort_Order'}         = $count;  ##实现按chr分组
        $hashBED{$chr}{'Chr_Bed_Region'}{$count} = $_;
    }
    close BED;
    return %hashBED;
}