package package::Math::basic;
use strict;
use warnings;

sub Pow{
	my ($a,$b)=@_;
	my $v=1;
	for my $i(1..$b){
		$v*=$a;
	}
	return $v;
}

sub numberFormat{
	my $num=shift @_;
	my $count=shift @_;
	my @split_num_re=split //,reverse($num);
	my $newNum="";
	for my $i(0..@split_num_re-1){
		$newNum.="," if($i!=0 and $i%$count==0);
		$newNum.=$split_num_re[$i];
	}
	return reverse($newNum);
}

sub min{
	my $a=shift @_;
	my $b=shift @_;
	if($a<$b){
		return $a;
	}else{
		return $b;
	}
}

sub max{
	my $a=shift @_;
	my $b=shift @_;
	if($a>$b){
		return $a;
	}else{
		return $b;
	}
}


1