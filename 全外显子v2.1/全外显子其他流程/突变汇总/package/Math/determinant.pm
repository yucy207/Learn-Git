package package::Math::determinant;
use strict;
use warnings;

###
# 求行列式的值
###
sub value{
	my $data=shift @_;
	my $del_I=shift @_;
	my $del_J=shift @_;
	my $row=1;
	my %newData=();
	foreach my $i(sort{$a <=> $b}keys %{$data}){
		next if($i==$del_I);
		my $col=1;
		foreach my $j(sort{$a <=> $b}keys %{$data->{$i}}){
			next if($j==$del_J);
			$newData{$row}{$col}=$data->{$i}{$j};$col++;
		}
		$row++;
	}
	my $n=keys %newData;
	if($n==1){
		return package::Math::basic::Pow(-1,$del_I+$del_J)*$newData{1}{1};
	}
	my $value=0;
	foreach my $j(sort{$a <=>$b}keys %{$newData{1}}){
		$value+=$newData{1}{$j}*value(\%newData,1,$j);
	}
	return package::Math::basic::Pow(-1,$del_I+$del_J)*$value;
}

###
# 求行列式的乘
###
sub multi{
	my $data1=shift @_;
	my $data2=shift @_;
	my $maxRow=0;
	foreach my $row(keys %{$data1}){
		$maxRow++;
	}
	my $maxCol=0;
	foreach my $row(keys %{$data2}){
		foreach my $col(keys %{$data2->{$row}}){
			$maxCol++;
		}
		last;
	}
	my %multi=();
	for my $i(1..$maxRow){
		for my $j(1..$maxCol){
			foreach my $col(keys %{$data1->{$i}}){
				$multi{$i}{$j}+=$data1->{$i}{$col}*$data2->{$col}{$j};
			}
		}
	}
	return %multi;
}



1