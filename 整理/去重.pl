##每一行的标识为第一列

while (<F>){
    my @arr=split "\t",$_;
    my $s=$arr[0];
    $hashtimes{$s}++;
    next if ($hashtimes{$s}>=2);
}
