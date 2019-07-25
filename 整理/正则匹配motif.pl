sub str_annotation_indel{
    my ($seq, $num) = @_;
    my $str_find = 0;
    my @specs = ([2,$num],
        [3,$num],
        [4,$num],
        [5,$num],
        [6,$num],
    ); 
    my $i = 0;
    my $str_seq_all = "";
    for($i=0; $i<scalar(@specs); $i++){
        my $motiflength = $specs[$i]->[0];
        my $minreps = $specs[$i]->[1] - 1;
        while($seq =~/(([gatc]{$motiflength})\2{$minreps,})/ig){  ##\2第2个捕获的()部分,不使用/2后面每次重复的$motiflength长序列内容可以不一致
            my $str_seq = $1;
            my $positon_end = pos($seq);
            my $positon_start = $positon_end - length($str_seq);
            if($positon_start<=52 && $positon_end>=48 ){
                $str_find++;
                $str_seq_all .= "$str_seq;";
            }
        }
    }
    my $result = "$str_find|$str_seq_all";
    return  $result;       
}