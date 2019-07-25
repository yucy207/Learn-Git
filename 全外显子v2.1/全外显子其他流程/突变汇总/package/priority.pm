package package::priority;
use strict;
use warnings;

sub Site{
	my ($hgmd,$sift,$g1000,$contperc,$esp6500,$snvq,$genoq,$homoHit,$geneskyHitFreq)=@_;
	my $Priority="Third";
	#Second
	if($snvq ne "L" and $genoq>=0.5){
		if($g1000<0.01 or $contperc<0.01){
			if($homoHit<3){
				$Priority="Second";
			}
		}
	}
	#First2
	if($g1000<0.001 or $contperc<0.01){
		if($esp6500<0.01){
			if($snvq ne "L" and $genoq>=0.5){
				if($homoHit==1 and $geneskyHitFreq<0.005){
					$Priority="First2";
				}
			}
		}
	}
	#First1
	$Priority="First1" if($hgmd==1);
	if($esp6500<0.01){
		if($snvq ne "L" and $genoq>=0.5){
			if($homoHit==1 and $geneskyHitFreq<0.005){
				if( ($g1000<0.001 or $contperc<0.01) and ($hgmd==1 or $sift==1) ){
					$Priority="First1";
				}
			}
		}
	}
	return $Priority;
}

sub Gene{
	my ($model,$tmppriority)=@_;
	my %hash_tmp=();
	if($model eq "Recessive"){
		foreach my $sample(keys %{$tmppriority}){
			my $Priority="";
			my $tmp_add=0;
			if(exists $tmppriority->{$sample}{"First1"}){$tmp_add+=$tmppriority->{$sample}{"First1"};}
			if($tmp_add>=2 and $Priority eq ""){$Priority="First1";}
			if(exists $tmppriority->{$sample}{"First2"}){$tmp_add+=$tmppriority->{$sample}{"First2"};}
			if($tmp_add>=2 and $Priority eq ""){$Priority="First2";}
			if(exists $tmppriority->{$sample}{"Second"}){$tmp_add+=$tmppriority->{$sample}{"Second"};}
			if($tmp_add>=2 and $Priority eq ""){$Priority="Second";}
			if(exists $tmppriority->{$sample}{"Third"}){$tmp_add+=$tmppriority->{$sample}{"Third"};}
			if($tmp_add>=2 and $Priority eq ""){$Priority="Third";}
			$hash_tmp{$Priority}++;
		}
	}
	else{
		foreach my $sample(keys %{$tmppriority}){
			if(exists $tmppriority->{$sample}{"First1"}){$hash_tmp{"First1"}++;}
			if(exists $tmppriority->{$sample}{"First2"}){$hash_tmp{"First2"}++;}
			if(exists $tmppriority->{$sample}{"Second"}){$hash_tmp{"Second"}++;}
			if(exists $tmppriority->{$sample}{"Third"}){$hash_tmp{"Third"}++;}
		}
	}
	my $new_Priority="";
	if(exists $hash_tmp{"First1"}){$new_Priority="First1";}
	elsif(exists $hash_tmp{"First2"}){$new_Priority="First2";}
	elsif(exists $hash_tmp{"Second"}){$new_Priority="Second";}
	elsif(exists $hash_tmp{"Third"}){$new_Priority="Third";}
	return $new_Priority;
}


1