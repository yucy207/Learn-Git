package package::CNV;
use strict;
use warnings;
use Statistics::Descriptive;

sub run{
    my $config=shift @_;
    my $hashtag=shift @_;
	my %hashok=();###判断该样本的该位点是否要输出，主要是control样本是否>0.3
	my %hashCNV=();
	print "Calculating CNV...";
	my @groups=();@groups=keys %{$config->{'Group'}} if(exists($config->{'Group'}));
	my @shows=();@shows=keys %{$config->{'Show'}} if(exists($config->{'Show'}));
	my @counts=(@groups,@shows);
	foreach my $count(@counts){
	    my $info=$config->{'Group'}{$count} if(exists($config->{'Group'}{$count}));
		   $info=$config->{'Show'}{$count} if(exists($config->{'Show'}{$count}));
	    my ($case,$control)=split /;/,$info;
		my @cases=split /,/,$case;
		my @controls=split /,/,$control;
		my @titles= keys %{$hashtag->{$cases[0]}};
		foreach my $title(@titles){
		    my ($chr,$tmp)=split /\|/,$title;
		    foreach my $sample(@cases){
			    my @cnvs;
			    foreach my $control(@controls){
				    if($hashtag->{$control}{$title}{'normc'}>0.3){
					   my $cnvvalue=$hashtag->{$sample}{$title}{'normc'}/$hashtag->{$control}{$title}{'normc'};
					   if($chr eq 'X' or $chr eq 'Y'){
					      if($config->{'Sample'}{$sample}{'sex'} eq $config->{'Sample'}{$control}{'sex'}){
					         push @cnvs,$cnvvalue;
					      }elsif($config->{'Sample'}{$sample}{'sex'} eq 'male'){
					         push @cnvs,$cnvvalue*2;
					      }else{
					         push @cnvs,$cnvvalue/2;
					      }
					   }else{
					      push @cnvs,$cnvvalue;
					   }
					}
				}
				if(exists($cnvs[0])){
				   my $value=calculateCNV(\@cnvs);
				   $hashCNV{$sample}{$title}=$value;
				   $hashok{$sample}{$title}=1;
				}else{				   
				   $hashCNV{$sample}{$title}=1;####如果control样本均值为0，则将该数值设为10
				   $hashok{$sample}{$title}=0;
				}
			}
		}
		
	}  
	print "OK\n";
	my @results=();
	push @results,\%hashCNV;
	push @results,\%hashok;
    return @results;
}

sub calculateCNV{
    my $data=shift @_;
	my $stat= Statistics::Descriptive::Full->new();
	$stat->add_data(@$data);
	my $median=$stat->median();
	my $mean=$stat->mean();
	return $mean if(@$data<=2);
	return $median if(@$data>2);
}



1