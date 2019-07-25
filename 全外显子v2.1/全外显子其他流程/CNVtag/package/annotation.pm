package package::annotation;
use strict;
use warnings;

sub run{
    my $config=shift @_;
	my $hashtag=shift @_;
	my @samples=keys %$hashtag;
	my $sample=$samples[0];
	my $dir=$config->{'Report'}."/CNVtag";
	mkdir $dir if(not -e $dir);
	my $library=$dir."/library";
	my $geneanno=$dir."/library.variant_function";
	my $functionanno=$dir."/library.exonic_variant_function";
	my $AnnovarBuild = $config->{'AnnovarBuild'};
	my $AnnovarDIR   = $config->{'annovar'};
	
	my %hashAnnotationList = package::annotation_db::get_annotation(); # # 获取注释列表
	my @anno_models = keys %hashAnnotationList;
	
	if(not -e $geneanno or -s $geneanno==0){
	    my @titles=keys %{$hashtag->{$sample}};
	    open OUT,">$library";
	    foreach my $title(@titles){
	       my ($chr,$start,$end)=split /\|/,$title;
		   print OUT "$chr\t$start\t$end\t0\t-\n";
	    }
	    close OUT;
		
        foreach my $anno_model(@anno_models) {
			my $pid = fork;
			if($pid == 0) {
			    next if($hashAnnotationList{$anno_model}{'From'} ne 'Annovar');
                my $AnnovarDB    = $hashAnnotationList{$anno_model}{'DBDir'}{'Human_hg19'}; # # cnvtag没有物种配置，默认均为hg19
                system "$AnnovarDIR/$hashAnnotationList{$anno_model}{'Para'} --buildver $AnnovarBuild $library $AnnovarDB";
				exit;
			}
		}
		while(my $pid = waitpid(-1, 0)>0){};   	   
	}
	my %hashanno;
	open IN,$geneanno;
	while(<IN>){
	    $_=~s/[\r\n]//g;
		my ($region,$gene1,$chr,$start,$end,$tmp)=split /\t/,$_,6;
		my ($gene)=$gene1=~/^([\.\-\_\w]+)/;
		my $title="$chr|$start|$end";
		$hashanno{$title}{'Gene Region'}=$region;
		$hashanno{$title}{'Gene'}=$gene;
	}
    close IN;
    open IN,$functionanno;
    while(<IN>){
    	$_=~s/[\r\n]//g;
		my ($line,$function,$ppvs,$chr,$start,$end,$tmp)=split /\t/,$_,7;
		my @annos;
		foreach my $ppv((split /,/,$ppvs)){
			my ($gene,$NMnumber,$Exon,$tmp)=split /:/,$ppv,4;
			next if(!defined $Exon);
            push @annos,"$gene:$NMnumber:$Exon";
		}
		my $title="$chr|$start|$end";
		$hashanno{$title}{'Gene Exon'}=join ",",@annos;  	
    }
    close IN;
    my @db=("cytoBand","dgvFreq","iscaPathGainCum","iscaPathLossCum","iscaLikelyPathogenic","iscaPathogenic","iscaCuratedPathogenic","CNVD","DECIPHER");
    foreach(@db){
	    open FILE,"$library.".$config->{"AnnovarBuild"}."_$_";
	    while(my $line=<FILE>){
		    $line=~ s/[\r\n]//g;
		    my @split_line=split /\t/,$line;
		    $split_line[1]=~ s/^Name\=//;
		    my $value=$split_line[1];
		    if($_ eq "dgvFreq"){
			   my ($gainAdd,$lossAdd,$allAdd)=(0,0,0);
			   foreach my $key(split /,/,$split_line[1]){
			       my ($gain,$loss,$all)=split /\;/,$key;
			       $gainAdd+=$gain;
			       $lossAdd+=$loss;
			       $allAdd+=$all;
			   }
			   $value=0;$value=($gainAdd+$lossAdd)/$allAdd if($allAdd>0);
		    }
			my $title="$split_line[2]|$split_line[3]|$split_line[4]";
			$hashanno{$title}{$_}=$value;
	    }
	    close FILE;
    }
	return %hashanno;
}


1