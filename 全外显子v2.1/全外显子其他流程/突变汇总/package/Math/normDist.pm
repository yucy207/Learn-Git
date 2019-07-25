package package::Math::normDist;
use strict;
use warnings;

###
# 正态分布
# 输入z值，得到单尾概率
###
sub pvalue{
	my $a=shift @_;
	my $pi = 3.14159265358979;
	my $p = 0.2316419;
	my $b1 = 0.31938153;
	my $b2 = -0.356563782;
	my $b3 = 1.781477937;
	my $b4 = -1.821255978;
	my $b5 = 1.330274429;
    my $x = abs($a);
    my $t = 1 / (1 + $p * $x);
    my $NormSDist = 1 - (1 / (sqrt(2 * $pi)) * exp(-package::Math::basic::Pow($a,2) / 2)) * ($b1 * $t + $b2 * package::Math::basic::Pow($t,2) + $b3 * package::Math::basic::Pow($t,3) + $b4 * package::Math::basic::Pow($t,4) + $b5 * package::Math::basic::Pow($t,5));
    $NormSDist = 1 - $NormSDist if($a<0);
	return $NormSDist;
}



1