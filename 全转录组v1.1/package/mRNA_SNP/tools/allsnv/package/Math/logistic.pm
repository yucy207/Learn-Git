package package::Math::logistic;
use strict;
use warnings;

###
# Logistic回归
# my %data=();
# %{$data{"id1"}}=("factor"=>0,"age"=>0,"DV"=>0,"REP"=>15);
# %{$data{"id2"}}=("factor"=>0,"age"=>1,"DV"=>0,"REP"=>10);
###
sub run{
	my $data=shift @_;
	my $ep=0.000000000001;
	my $maxLoop=20;
	my %logisticKey=Logistic_Key($data);
	return if(keys %logisticKey<=1);
	# 初始化 Beita
	my %Beita=package::Math::leastSquare::run(\%logisticKey,$data);
	return if(keys %Beita==0);
	# 最大似然+牛顿迭代
	my %SE=();
	my $Pre="";
	while(($Pre eq "" or abs($Pre)>$ep) and $maxLoop>0){
		# (|v1 v2|+|Lamda*v1 0|)-1 |a b| |v5|
		# (|v3 v4|+|0 Lamda*v4|)  =|c d|*|v6|
		my %Derivative=Logistic_Likelihood(\%Beita,\%logisticKey,$data);
		return if(keys %Derivative==0);
		foreach my $key(sort{$a cmp $b}keys %logisticKey){
			$Pre=$Derivative{"Result"}{$key};
			$Beita{$key}=$Beita{$key}-$Pre;
			$SE{$key}=sqrt(abs($Derivative{"SE"}{$key}{$key}));
		}
		$maxLoop--;
	}
	return if($maxLoop==0);
	my %OR=();
	foreach my $key(keys %logisticKey){
		my $OR=exp($Beita{$key});
		my $SE=$SE{$key};
		return if($SE==0);
		my $OR95CI_Low=$OR/exp(1.96*$SE);
		my $OR95CI_Up=$OR*exp(1.96*$SE);
		my $Z=abs(log($OR)/$SE);
		my $Pvalue=2-2*package::Math::normDist::pvalue($Z);
		$OR{$key}{"OR"}=$OR;
		$OR{$key}{"OR95CI_Low"}=$OR95CI_Low;
		$OR{$key}{"OR95CI_Up"}=$OR95CI_Up;
		$OR{$key}{"Pvalue"}=$Pvalue;
	}
	return %OR;
}
###
# Logistic参数，常数CV，因变量DV，重复次数REP
###
sub Logistic_Key{
	my $data=shift @_;
	# 基础测试
	my %variable=();
	foreach my $id(keys %{$data}){
		return if(not exists $data->{$id}{"DV"});
		return if(not exists $data->{$id}{"REP"});
		return if(exists $data->{$id}{"CV"});
		# 变量
		foreach my $v(keys %{$data->{$id}}){
			next if($v eq "DV" or $v eq "REP");
			$variable{$v}{$data->{$id}{$v}}=1;
		}
	}
	# 删除恒定变量
	foreach my $v(keys %variable){
		delete $variable{$v} if(keys %{$variable{$v}}<=1);
	}
	# 数据测试
	foreach my $v(keys %variable){
		foreach my $id(keys %{$data}){
			return if(not exists $data->{$id}{$v});
		}
	}
	$variable{"CV"}=1;
	return %variable;
}
###
# Logistic极大似然法
###
sub Logistic_Likelihood{
	my $Beita=shift @_;
	my $logisticKey=shift @_;
	my $data=shift @_;
	# 通用常数
	my %constant=();
	foreach my $id(keys %{$data}){
		my $beishuchu=0;
		foreach my $key1(keys %{$logisticKey}){
			my $x1=1;$x1=$data->{$id}{$key1} if(exists $data->{$id}{$key1});
			$beishuchu+=$Beita->{$key1}*$x1;
		}
		$constant{$id}{"beishuchu"}=exp($beishuchu);
		$constant{$id}{"chushu"}=1+exp($beishuchu);
	}
	# 二阶导数矩阵生成
	my %Derivative=();
	my $row=1;
	foreach my $key1(sort{$a cmp $b}keys %{$logisticKey}){
		my $col=1;
		foreach my $key2(sort{$a cmp $b}keys %{$logisticKey}){
			foreach my $id(keys %{$data}){
				my $x1=1;$x1=$data->{$id}{$key1} if(exists $data->{$id}{$key1});
				my $x2=1;$x2=$data->{$id}{$key2} if(exists $data->{$id}{$key2});
				my $beishuchu=$constant{$id}{"beishuchu"};
				my $chushu=$constant{$id}{"chushu"};
				# 转置的矩阵
				for my $rep(1..$data->{$id}{"REP"}){
					$Derivative{"SecondDerivative"}{$col}{$row}+=-$x1*$x2*$beishuchu/($chushu*$chushu);
				}
			}
			$col++;
		}
		$row++;
	}
	foreach my $id(keys %{$data}){
		my $beishuchu=$constant{$id}{"beishuchu"};
		my $chushu=$constant{$id}{"chushu"};
		foreach my $key1(keys %{$logisticKey}){
			my $x1=1;$x1=$data->{$id}{$key1} if(exists $data->{$id}{$key1});
			for my $rep(1..$data->{$id}{"REP"}){
				$Derivative{"Derivative"}{$key1}+=$x1*($data->{$id}{"DV"}-$beishuchu/$chushu);
			}
		}
	}
	# 矩阵求逆
	my $ArrayValue=package::Math::determinant::value($Derivative{"SecondDerivative"},0,0);
	return if($ArrayValue==0);
	$row=1;
	foreach my $key1(sort{$a cmp $b}keys %{$logisticKey}){
		my $col=1;
		foreach my $key2(sort{$a cmp $b}keys %{$logisticKey}){
			my $SE=package::Math::determinant::value($Derivative{"SecondDerivative"},$row,$col)/$ArrayValue;
			$Derivative{"Result"}{$key1}+=$SE*$Derivative{"Derivative"}{$key2};
			$Derivative{"SE"}{$key1}{$key2}=$SE;
			$col++;
		}
		$row++;
	}
	return %Derivative;
}




1