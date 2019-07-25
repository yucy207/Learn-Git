package package::meancontrol;
use strict;
use warnings;
use Statistics::Descriptive;

sub run{
    my $config=shift @_;
    my $hashtag=shift @_;
	my @groups=keys %{$config->{'Group'}};
	my %hashcontrol=();
	print "Calculating Control Samples Mean normalized_coverage...";
	foreach my $group(@groups){
	    my ($case,$control)=split /;/,$config->{'Group'}{$group};
		my @controls=split /,/,$control;
		my @titles=keys %{$hashtag->{$controls[0]}};
		foreach my $title(@titles){
		    my @cover=();
		    foreach my $control(@controls){push @cover,$hashtag->{$control}{$title}{'normc'};}
			$hashcontrol{$group}{$title}=average(\@cover);# 均值
			# $hashcontrol{$group}{$title}=median(\@cover);# 中值
		}
	}    
	print "OK\n";
    return %hashcontrol;
}

sub average{
    my $data=shift @_;
	my $num=@$data;
	my $sum=0;
	foreach(@$data){$sum=$sum+$_;}
	my $mean=$sum/$num;
	return $mean;
}

sub median{
    my $data=shift @_;
	my $stat= Statistics::Descriptive::Full->new();
	$stat->add_data(@$data);
	my $median=$stat->median();
	return $median;
}

1