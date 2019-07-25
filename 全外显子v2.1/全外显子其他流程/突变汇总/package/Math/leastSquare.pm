package package::Math::leastSquare;
use strict;
use warnings;

sub run{
	my $logisticKey=shift @_;
	my $data=shift @_;
	# 二阶矩阵生成
	my %Array=();
	my $row=1;
	foreach my $id(sort{$a cmp $b}keys %{$data}){
		my $col=1;
		foreach my $key1(sort{$a cmp $b}keys %{$logisticKey}){
			my $value=1;$value=$data->{$id}{$key1} if($key1 ne "CV");
			$Array{"X"}{$row}{$col}=$value;
			$Array{"Xtrans"}{$col}{$row}=$value;
			$col++;
		}
		$Array{"Y"}{$row}{1}=$data->{$id}{"DV"};
		$row++;
	}
	# Xtrans*X
	%{$Array{"Xtrans*X"}}=package::Math::determinant::multi($Array{"Xtrans"},$Array{"X"});
	# 矩阵求逆
	my $ArrayValue=package::Math::determinant::value($Array{"Xtrans*X"},0,0);
	return if($ArrayValue==0);
	foreach my $row(keys %{$Array{"Xtrans*X"}}){
		foreach my $col(keys %{$Array{"Xtrans*X"}{$row}}){
			$Array{"Reverse"}{$row}{$col}=package::Math::determinant::value($Array{"Xtrans*X"},$row,$col)/$ArrayValue;
		}
	}
	# Beita=(Xtrans*X)-1*Xtrans*Y
	%{$Array{"(Xtrans*X)-1*Xtrans"}}=package::Math::determinant::multi($Array{"Reverse"},$Array{"Xtrans"});
	my %tmp=package::Math::determinant::multi($Array{"(Xtrans*X)-1*Xtrans"},$Array{"Y"});
	$row=1;
	my %Beita=();
	foreach my $key1(sort{$a cmp $b}keys %{$logisticKey}){
		$Beita{$key1}=$tmp{$row}{1};$row++;
	}
	return %Beita;
}



1