package package::Math::p1p2p3p4;
use strict;
use warnings;
use Statistics::ChisqIndep;
use Text::NSP::Measures::2D::Fisher2::twotailed;

###
# 计算4种情况下case/control的卡方p值
# input case_HOMR,case_HET,case_HOMA,control_HOMR,control_HET,control_HOMA
###
sub run{
	my ($hc00,$hc10,$hc11,$hn00,$hn10,$hn11)=@_;
	my $case_refs=$hc00*2+$hc10;
	my $case_alts=$hc11*2+$hc10;
	my $control_refs=$hn00*2+$hn10;
	my $control_alts=$hn11*2+$hn10;
	my $h00=$hc00+$hn00;
	my $h10=$hc10+$hn10;
	my $h11=$hc11+$hn11;
	my (@alts,@a1,@a2,@a3)=();
	push(@alts,[$case_refs,$case_alts]);
	push(@alts,[$control_refs,$control_alts]);
	push(@a1,[$hc00,$hc10,$hc11]);
	push(@a1,[$hn00,$hn10,$hn11]);
	push(@a2,[$hc00,$hc10+$hc11]);
	push(@a2,[$hn00,$hn10+$hn11]);	
	push(@a3,[$hc00+$hc10,$hc11]);
	push(@a3,[$hn00+$hn10,$hn11]);	 
	my $chi4 = new Statistics::ChisqIndep;
	$chi4->load_data(\@alts);
	my $p4= $chi4->{p_value} ;
	if($alts[0][0]<5 || $alts[0][1]<5 || $alts[1][0]<5 || $alts[1][1]<5){
		my $n11= $alts[0][0];
		my $n1p =$alts[0][0]+$alts[0][1];
		my $np1 = $alts[0][0]+$alts[1][0];
		my $npp = $alts[0][0]+$alts[0][1]+$alts[1][0]+$alts[1][1];
		$p4= calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	}
	my $chi1 = new Statistics::ChisqIndep;
	$chi1->load_data(\@a1);
	my $p1= $chi1->{p_value} ;
	my $chi2 = new Statistics::ChisqIndep;
	$chi2->load_data(\@a2);
	my $p2= $chi2->{p_value} ;
	if($a2[0][0]<5 || $a2[0][1]<5 || $a2[1][0]<5 || $a2[1][1]<5) {
		my $n11= $a2[0][0];
		my $n1p =$a2[0][0]+$a2[0][1];
		my $np1 = $a2[0][0]+$a2[1][0];
		my $npp = $a2[0][0]+$a2[0][1]+$a2[1][0]+$a2[1][1];
		$p2=calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	}
	my $chi3 = new Statistics::ChisqIndep;
	$chi3->load_data(\@a3);
	my $p3= $chi3->{p_value} ; 
	if($a3[0][0]<5 || $a3[0][1]<5 || $a3[1][0]<5 || $a3[1][1]<5) {
		my $n11= $a3[0][0];
		my $n1p =$a3[0][0]+$a3[0][1];
		my $np1 = $a3[0][0]+$a3[1][0];
		my $npp = $a3[0][0]+$a3[0][1]+$a3[1][0]+$a3[1][1];
		$p3= calculateStatistic( n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
	}
	my $freq=0;$freq=($h10/2+$h11)/($h10+$h11+$h00) if($h10+$h11+$h00>0);
	if($freq>0.5){
		my	$p_t=$p3;
		$p3=$p2;
		$p2=$p_t;
	}
	return ($p1,$p2,$p3,$p4);
}


1