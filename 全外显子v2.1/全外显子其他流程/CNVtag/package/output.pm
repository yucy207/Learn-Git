package package::output;
use strict;
use warnings;

sub run{
    my $config=shift @_;
    my $hashtag=shift @_;
    my $hashCNV=shift @_;
	my $hashok=shift @_;
	my $hashanno=shift @_;
	my $Reportdir=$config->{'Report'};
	my @groups=keys %{$config->{'Group'}};
	my @samples=sort keys %$hashCNV;
	my @controlsamples=getcontrolsample($config);
	my @titles=keys %{$hashCNV->{$samples[0]}};	
	my %hashsort=sorttitle(\@titles);###按染色体与位置排序
	my %hashCNVtype=CNVtype($hashCNV,$config);###把拷贝数转换为Gain与Loss
	my %hashoutput=tongji(\%hashCNVtype,\%hashsort,$hashanno,$config);# 输出过滤
	print "Output $Reportdir/document/CNVtag.xlsx\n";
	mkdir "$Reportdir/document/" if(not -e "$Reportdir/document/");
	mkdir "$Reportdir/document/3_CNV" if(not -e "$Reportdir/document/3_CNV");
	my $workbook=Excel::Writer::XLSX->new("$Reportdir/document/3_CNV/CNVtag.xlsx");
	my %format=package::format::addformat($workbook);
	my $report=$workbook->add_worksheet("CNV");	
	my @db=('Gene','Gene Region','Gene Exon',"cytoBand","dgvFreq","iscaPathGainCum","iscaPathLossCum","iscaLikelyPathogenic","iscaPathogenic","iscaCuratedPathogenic","CNVD","DECIPHER");
    my @cnvsamples=@samples;#给出拷贝数的样本
    my @tagsamples=(@samples,@controlsamples);# 给出原始tag文件中的数据
    my @title_cnvsamples;map{push @title_cnvsamples,"$_"."-CNV"}@cnvsamples;# CNV样本表头
    my @title_controls;map{push @title_controls,"$_"."-Control"}@controlsamples;# Control样本表头
    my @toptitles=('Chr','Start','End','Length','%gc',@db,@title_cnvsamples,@samples,@title_controls);
	my $row=0;
    foreach my $i(0..$#toptitles){
	    $report->write($row,$i,$toptitles[$i],$format{"title"});
	}
    $row++;
    @titles=gettitle(\%hashsort);
    foreach my $title(@titles){
	    my $mark=0;###判断该区域是否输出
		foreach my $sample(@samples){$mark=1 if(exists($hashoutput{$title}{$sample}));}
		next if($mark==0);
		my @OriginalDatas=();my @colors=();
		my ($chr,$start,$end)=split /\|/,$title;
		push @OriginalDatas,$chr;push @colors,'normal';
		push @OriginalDatas,$start;push @colors,'normal';
		push @OriginalDatas,$end;push @colors,'normal';
		push @OriginalDatas,$hashtag->{$samples[0]}{$title}{'length'};push @colors,'normal';
		push @OriginalDatas,$hashtag->{$samples[0]}{$title}{'gc'};push @colors,'normal';
		foreach(@db){
		    my $value=(exists($hashanno->{$title}{$_})) ? $hashanno->{$title}{$_} : " ";
			push @OriginalDatas,$value;push @colors,'normal';
		}
		foreach my $sample(@cnvsamples){
		    my $thiscolor='normal';
			my $CNV=$hashCNV->{$sample}{$title};
			my $gainloss=(exists($hashCNVtype{$sample}{$title})) ? $hashCNVtype{$sample}{$title} : "";
			$thiscolor='green' if($gainloss eq 'Gain');
			$thiscolor='red' if($gainloss eq 'Loss');
			$CNV="ControlError" if($hashok->{$sample}{$title}==0); # 对照样本有问题
			push @colors,$thiscolor;
            push @OriginalDatas,$CNV;		
		}
		foreach my $sample(@tagsamples){
		   push @OriginalDatas,$hashtag->{$sample}{$title}{'normc'};push @colors,'normal';
		}
		
        foreach my $i(0..$#OriginalDatas){
	        $report->write($row,$i,$OriginalDatas[$i],$format{$colors[$i]});
	    }
        $row++;		
	}
	my $Readme=$workbook->add_worksheet("READ ME");
    package::readme::run($Readme,\%format,"readme.txt");	
}

sub gettitle{
    my $hashsort=shift @_;
	my @tmpchrs=keys %$hashsort;
	my @chrs=();
	my @normalchrs=();
	my @otherchrs=();
	foreach my $chr(@tmpchrs){
	    if($chr=~/\d/){
		   push @normalchrs,$chr;
		}else{
		   push @otherchrs,$chr;
		}
	}
	@normalchrs=sort {$a<=>$b} @normalchrs;
	@otherchrs=sort @otherchrs;
	@chrs=(@normalchrs,@otherchrs);	
	my @titles=();
	foreach my $chr(@chrs){
	    my @thistitles=split /,/,$hashsort->{$chr};
		@titles=(@titles,@thistitles);
	}
	return @titles;
}

sub tongji{
    my $hashCNVtype=shift @_;
	my $hashsort=shift @_;
	my $hashanno=shift @_;
	my $config=shift @_;
	my @groups=keys %{$config->{'Group'}};
	my %hashoutput=();
	my @samples=sort keys %$hashCNVtype;
	my @tmpchrs=keys %$hashsort;
	my $singlecontinue=(exists($config->{'SingleContinue'}) and $config->{'SingleContinue'}=~/\d/) ? $config->{'SingleContinue'} : 2;
	my @chrs=();
	my @normalchrs=();
	my @otherchrs=();
	foreach my $chr(@tmpchrs){
	    if($chr=~/\d/){
		   push @normalchrs,$chr;
		}else{
		   push @otherchrs,$chr;
		}
	}
	@normalchrs=sort {$a<=>$b} @normalchrs;
	@otherchrs=sort @otherchrs;
	@chrs=(@normalchrs,@otherchrs);
	foreach my $group(@groups){
	    my ($case,$control)=split /;/,$config->{'Group'}{$group};
		my @cases=split /,/,$case;
		if(@cases==1){##单个样本，看基因上连续2个片段
		    CNVContinue($hashCNVtype,\%hashoutput,\@chrs,$hashsort,\@cases,$hashanno,$singlecontinue);
		}else{##多个样本，看多个样本是否均满足相同的变异
		    CNVMulti($hashCNVtype,\%hashoutput,\@chrs,$hashsort,\@cases);
		}
	}
    return %hashoutput;
}

sub CNVContinue{
    my $hashCNVtype=shift @_;
	my $hashoutput=shift @_;
	my $chrs=shift @_;
	my $hashsort=shift @_;
	my $samples=shift @_;
	my $hashanno=shift @_;
	my $continue=shift @_;# 连续相同类型的数量
	my $sample=$$samples[0];# 样本
	foreach my $chr(@$chrs){
	    my @titles=split /,/,$hashsort->{$chr};
		my $thistype="";###类型typetitles
		my $thisgene="";
		my @typetitles=();###连续target
		foreach my $title(@titles){
		    my $type=$hashCNVtype->{$sample}{$title};
			my $gene=" ";$gene=$hashanno->{$title}{'Gene'};
			if($type eq $thistype and $gene eq $thisgene){# 相同,列入候选
			   push @typetitles,$title;
			}else{ # 不同
			   #（1）检测是否需要保存，并保存
			   keepContinue(\@typetitles,$hashoutput,$sample) if(@typetitles>=$continue);
			   #（2）清空，并检测当前title是否可以导入候选
			   $thistype="";
			   $thisgene="";
			   @typetitles=();
			   if($type ne 'Normal'){
			      $thistype=$type;
				  $thisgene=$gene;
				  push @typetitles,$title;
			   }
			}
        }
        keepContinue(\@typetitles,$hashoutput,$sample) if(@typetitles>=$continue);		
	}		
}

sub keepContinue{
    my $typetitles=shift @_;
	my $hashoutput=shift @_;
	my $sample=shift @_;
	map{$hashoutput->{$_}{$sample}++;}@$typetitles;
}

sub CNVMulti{
    my $hashCNVtype=shift @_;
	my $hashoutput=shift @_;
	my $chrs=shift @_;
	my $hashsort=shift @_;
	my $samples=shift @_;
	foreach my $chr(@$chrs){
	    my @titles=split /,/,$hashsort->{$chr};
		foreach my $title(@titles){
		    my %hash;
			map{$hash{$hashCNVtype->{$_}{$title}}++;}@$samples;
			my @types=keys %hash;
			if(@types==1 and $types[0] ne 'Normal'){
			   map{$hashoutput->{$title}{$_}++;}@$samples;
			}		     
		}
	}
}

sub sorttitle{
    my $titles=shift @_;
	my %hashdata=();
	foreach my $title(@$titles){
	    my ($chr,$start,$end)=split /\|/,$title;
		$hashdata{$chr}{$start}{$end}=$title;
	}
	my %hashsort=();
	my @chrs=keys %hashdata;
	foreach my $chr(@chrs){
	    my @thistitles=();
	    my @starts=sort {$a<=>$b} keys %{$hashdata{$chr}};
		foreach my $start(@starts){
		    my @ends=sort {$a<=>$b} keys %{$hashdata{$chr}{$start}};
			foreach my $end(@ends){
			    push @thistitles,$hashdata{$chr}{$start}{$end}
			}
		}
		my $str=join ",",@thistitles;
		$hashsort{$chr}=$str;
	}
    return %hashsort;
}

sub CNVtype{
    my $hashCNV=shift @_;
    my $config=shift @_;
    my $Gain=(exists($config->{'Gain'}) and $config->{'Gain'}=~/\d/) ? $config->{'Gain'} : 1.8;
    my $Loss=(exists($config->{'Loss'}) and $config->{'Loss'}=~/\d/) ? $config->{'Loss'} : 0.6;
	my %hashCNVtype=();
	my @samples=sort keys %$hashCNV; 
	my @titles=keys %{$hashCNV->{$samples[0]}};	
	foreach my $title(@titles){
	    foreach my $sample(@samples){
		    my $value=$hashCNV->{$sample}{$title};
			if($value>$Gain){###大于1.8设为多拷贝
			   $hashCNVtype{$sample}{$title}='Gain';
			}elsif($value<$Loss){###小于0.6设为低拷贝
			   $hashCNVtype{$sample}{$title}='Loss';
			}else{
			   $hashCNVtype{$sample}{$title}='Normal';
			}
		}
	}
    return %hashCNVtype;
}

sub getcontrolsample{
    my $config=shift @_;
	my @groups=();@groups=keys %{$config->{'Group'}} if(exists($config->{'Group'}));
	my @shows=();@shows=keys %{$config->{'Show'}} if(exists($config->{'Show'}));
	my @counts=(@groups,@shows);
	my %hash;
	foreach my $count(@counts){
	    my $info=$config->{'Group'}{$count} if(exists($config->{'Group'}{$count}));
		   $info=$config->{'Show'}{$count} if(exists($config->{'Show'}{$count}));
	    my ($case,$control)=split /;/,$info;
		my @controls=split /,/,$control;
		map{$hash{$_}=1 if($_=~/\w/);}@controls;
	}
    my @samples=sort keys %hash;
    return @samples;	
}


1