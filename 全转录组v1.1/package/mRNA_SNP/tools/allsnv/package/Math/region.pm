package package::Math::region;
use strict;
use warnings;

# input $hash{chr}{pos}
sub chain{
	my $data=shift @_;
	my %hash=();
	foreach my $chr(keys %{$data}){
		my $pre="";
		foreach my $pos(sort{$a<=>$b}keys %{$data->{$chr}}){
			if($pre eq ""){
				$pre=$pos;
				next;
			}
			$hash{$chr}{$pre}=$pos;
			$pre=$pos;
		}
	}
	return %hash;
}

# input %{$hash{chr}{weiyi}}=("s"=>start,"e"=>end,"value"=>value)
sub annotation{
	my $chain=shift @_;
	my $data=shift @_;
	my %hash=();
	foreach my $chr(keys %{$data}){
		foreach my $weiyi(keys %{$data->{$chr}}){
			my $pos=$data->{$chr}{$weiyi}{"s"};
			while(exists $chain->{$chr}{$pos} and $chain->{$chr}{$pos}<=$data->{$chr}{$weiyi}{"e"}){
				$hash{$chr}{$pos}{$chain->{$chr}{$pos}}{$data->{$chr}{$weiyi}{"value"}}++;
				$pos=$chain->{$chr}{$pos};
			}
		}
	}
	return %hash;
}

# 合并blast不唯一的各个区域
sub blastMerge{
	my $annotation=shift @_;
	my %hash=();
	foreach my $chr(keys %{$annotation}){
		foreach my $start(sort{$a<=>$b}keys %{$annotation->{$chr}}){
			foreach my $end(keys %{$annotation->{$chr}{$start}}){
				my $count=0;
				foreach my $value(keys %{$annotation->{$chr}{$start}{$end}}){
					$count+=$annotation->{$chr}{$start}{$end}{$value};
				}
				$hash{$chr}{$start}{"end"}=$end;
				if($count>1){
					$hash{$chr}{$start}{"value"}="Repeat";
				}else{
					$hash{$chr}{$start}{"value"}="Unknown";
				}
			}
		}
	}
	my %del=();
	foreach my $chr(keys %hash){
		my ($preStart,$preEnd,$preValue)=("","","");
		foreach my $start(sort{$a<=>$b}keys %{$hash{$chr}}){
			if($hash{$chr}{$start}{"end"}<=$start){
				$del{$chr}{$start}=1;
				next;
			}
			if($preValue eq ""){
				$preValue=$hash{$chr}{$start}{"value"};
				$preStart=$start;
				$preEnd=$hash{$chr}{$start}{"end"};
				next;
			}
			if($preValue eq $hash{$chr}{$start}{"value"} and $preEnd+1>=$start){
				$hash{$chr}{$preStart}{"end"}=$hash{$chr}{$start}{"end"};
				$preEnd=$hash{$chr}{$start}{"end"};
				$del{$chr}{$start}=1;
				next;
			}
			$preValue=$hash{$chr}{$start}{"value"};
			$preStart=$start;
			$preEnd=$hash{$chr}{$start}{"end"};
		}
	}
	foreach my $chr(keys %del){
		foreach my $start(keys %{$del{$chr}}){
			delete $hash{$chr}{$start};
		}
	}
	return %hash;
}


1