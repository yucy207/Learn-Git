package package::readtag;
use strict;
use warnings;

sub run{
    my $config=shift @_;
	my @groups=();@groups=keys %{$config->{'Group'}} if(exists($config->{'Group'}));
	my @shows=();@shows=keys %{$config->{'Show'}} if(exists($config->{'Show'}));
    my %hash;
	foreach my $group(@groups){
	    my @tmps=split /[,;]/,$config->{'Group'}{$group};
		map{$hash{$_}++ if($_=~/\w/);}@tmps;
	}
	foreach my $show(@shows){
	    my @tmps=split /[,;]/,$config->{'Show'}{$show};
		map{$hash{$_}++ if($_=~/\w/);}@tmps;
	}
	my @samples=sort keys %hash;	
	my %hashtag=();
	foreach my $sample(@samples){
	    my $tagfile=$config->{'Sample'}{$sample}{'tag'};
		print "Reading $tagfile...";
		my $line=0;
		open TAG,"$tagfile";
		while(<TAG>){
		    $line++;
			next if($line==1);
			$_=~s/[\r\n]//g;
			my ($chr,$start,$end,$length,$name,$gc,$meanc,$normc,$tmp)=split /\t/,$_,9;
			my $title="$chr|$start|$end";
			$hashtag{$sample}{$title}{'length'}=$length;
		    $hashtag{$sample}{$title}{'name'}=$name;
			$hashtag{$sample}{$title}{'gc'}=$gc;
			$hashtag{$sample}{$title}{'meanc'}=$meanc;
			$hashtag{$sample}{$title}{'normc'}=$normc;
		}
		close TAG;
		print "OK\n";		
	}   
    return %hashtag;
}


1