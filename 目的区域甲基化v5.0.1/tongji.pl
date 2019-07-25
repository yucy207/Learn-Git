$|=1;
use strict;
use warnings;
use Encode;
use Excel::Writer::XLSX;
use Statistics::R;
my $tab="##########";
# 读取配置
die "Usage perl $0 ConfigFile\n" if(@ARGV!=1);
my $configFile=shift @ARGV;
my %config=ReadConfig($configFile);
my @samples=GetSample(\%config,'case','control');
# 开始流程
system "clear";
print "$tab Start TongJi Pipline ".ymd()." ".hms()." $tab\n";
die "[ERR] Loss Report\n" if(not exists $config{"Report"} or not -e $config{"Report"});
# Excel
mkdir $config{"Report"}."/document" if(not -e $config{"Report"}."/document");
mkdir $config{"Report"}."/document/01_QC" if(not -e $config{"Report"}."/document");
print $config{"Report"}."/document/01_QC/Data_statistics.xlsx\n";
mkdir $config{"Report"}."/document/01_QC";
my $workbook=Excel::Writer::XLSX->new($config{"Report"}."/document/01_QC/Data_statistics.xlsx");
my %format=Format($workbook);

## fastqc文件拷贝
prepare_QC_image(\%config, \@samples); 

## fastqc_clean 文件拷贝
prepare_QC_Clean_image(\%config, \@samples) if(-e "$config{'Report'}/fastqc_clean");

fastqc(\%config,$workbook,\%format,\@samples);# 下机数据质控

QCCleanData(\%config,$workbook,\%format,\@samples);# 过滤后的数据质控

TargetCoverage(\%config,$workbook,\%format,\@samples);# 目标区域数据量统计，PCR

Metrics(\%config,$workbook,\%format,\@samples);# 目标区域数据量统计,Picard

Tag(\%config,$workbook,\%format,\@samples);# tag文件输出

ReadMe($workbook,\%format,"readme/readme_tongji.txt");
########################## FUNCTION ########################
sub ReadMe{
	my ($workbook,$format,$file)=@_;
    # package路径
    my $package_dir = Cwd::abs_path(get_dirname(File::Spec->rel2abs(__FILE__)));
    return if(checkFileSizeOK("$package_dir/$file")==0);    
	my $sheet=$workbook->add_worksheet("Read Me");
	$sheet->set_row(0, 65);
	$sheet->set_column('A:A', 35);
	$sheet->set_column('B:B', 25);
	$sheet->set_column('C:C', 110);
	my $row=0;
	open FILE,"$package_dir/$file";
	while(<FILE>){
		$_=~ s/\^/\n/g;
		my @split_line=split /\t/,$_;
		my $col=0;
		foreach(@split_line)
		{
			my $text=decode("UTF-8",$_);
			if($row==0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme1'});}
			if($row==0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme2'});}
			if($row==0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme2tmp'});}
			if($row>0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme3'});}
			if($row>0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme4'});}
			if($row>0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme5'});}
			$col++;
		}
		$row++;
	}
	close FILE;
}
# RRBS结果统计
sub RRBSSummary{
    my $config   = shift @_;
    my $workbook = shift @_;
    my $format   = shift @_;
    my $samples  = shift @_;
    print "Process Summary...";
    my $Report=$config->{'Report'};
    my %hashCheckResult=CheckFile($samples,"#####","$Report/bismark/#####_PE_report.txt");# 检查文件
    my $isError=ShowError(\%hashCheckResult,"\n[WARN] Loss RRBS PE_report, Sample:");# 错误信息输出
    # return if($isError);
    # 读数据
    my %hashsummary;
    foreach my $sample(@samples){
        my $SummaryFile="$Report/bismark/$sample"."_PE_report.txt";
        next if(checkFileSizeOK($SummaryFile)==0);
        open SUMMARY,$SummaryFile;
        my $count=0;
        while(<SUMMARY>){
            $count++;
            $_=~s/[\r\n]//g;
            next if($_!~/\w/);
            my @data=split /\t/,$_;
            next if($#data==0);
            $hashsummary{$sample}{'Total pair reads'}=$data[1] if($data[0]=~/Sequence pairs analysed in total/);
            $hashsummary{$sample}{'Unique Best Hit'}=$data[1] if($data[0]=~/Number of paired-end alignments with a unique best hit:/);
            $hashsummary{$sample}{'mapping'}=$data[1] if($data[0]=~/Mapping efficiency:/);
            $hashsummary{$sample}{'No Alignments'}=$data[1] if($data[0]=~/Sequence pairs with no alignments/);
            $hashsummary{$sample}{'Not Map Uniquely'}=$data[1] if($data[0]=~/Sequence pairs did not map uniquely/);
            $hashsummary{$sample}{'totalC'}=$data[1] if($data[0]=~/Total number of C's analysed:/);
            $hashsummary{$sample}{'MCinCPG'}=$data[1] if($data[0]=~/Total methylated C's in CpG context:/);
            $hashsummary{$sample}{'MCinCHG'}=$data[1] if($data[0]=~/Total methylated C's in CHG context:/);
            $hashsummary{$sample}{'MCinCHH'}=$data[1] if($data[0]=~/Total methylated C's in CHH context:/);
            $hashsummary{$sample}{'unMCinCPG'}=$data[1] if($data[0]=~/Total unmethylated C's in CpG context:/);
            $hashsummary{$sample}{'unMCinCHG'}=$data[1] if($data[0]=~/Total unmethylated C's in CHG context:/);
            $hashsummary{$sample}{'unMCinCHH'}=$data[1] if($data[0]=~/Total unmethylated C's in CHH context:/);
            $hashsummary{$sample}{'percentinCPG'}=$data[1] if($data[0]=~/C methylated in CpG context:/);
            $hashsummary{$sample}{'percentinCHG'}=$data[1] if($data[0]=~/C methylated in CHG context:/);
            $hashsummary{$sample}{'percentinCHH'}=$data[1] if($data[0]=~/C methylated in CHH context:/);        
        }
        close SUMMARY;          
    }
    my $sheet=$workbook->add_worksheet("summary");
    $sheet->set_row(0, 60);
    $sheet->set_column(0,0,30);
    $sheet->set_column(1,scalar (@samples),12);
    my $row=0;
    my @titles=("Sample",@$samples);
    foreach my $col(0..$#titles){
        $sheet->write($row,$col,$titles[$col],$format{"title"});
    }
    $row++;
    my @rowtitles;
    push @rowtitles,"Total pair reads";
    push @rowtitles,"Unique Best Hit";
    push @rowtitles,"Mapping efficiency";   
    push @rowtitles,"No Alignments";    
    push @rowtitles,"Not Map Uniquely";
    push @rowtitles,"Total number of C's analysed";
    push @rowtitles,"Total methylated C's in CpG context";
    push @rowtitles,"Total methylated C's in CHG context";
    push @rowtitles,"Total methylated C's in CHH context";
    push @rowtitles,"Total unmethylated C's in CpG context";
    push @rowtitles,"Total unmethylated C's in CHG context";
    push @rowtitles,"Total unmethylated C's in CHH context";
    push @rowtitles,"C methylated in CpG context";
    push @rowtitles,"C methylated in CHG context";
    push @rowtitles,"C methylated in CHH context";
    foreach my $rowtitle(@rowtitles){
        my $col=0;
        $sheet->write($row,$col,$rowtitle,$format{"normal"});
        $col++;
        foreach my $sample(@samples){
            $sheet->write($row,$col,$hashsummary{$sample}{'Total pair reads'},$format{"normal"}) if($rowtitle eq "Total pair reads");
            $sheet->write($row,$col,$hashsummary{$sample}{'Unique Best Hit'},$format{"normal"}) if($rowtitle eq "Unique Best Hit");
            $sheet->write($row,$col,$hashsummary{$sample}{'mapping'},$format{"normal"}) if($rowtitle eq "Mapping efficiency");  
            $sheet->write($row,$col,$hashsummary{$sample}{'No Alignments'},$format{"normal"}) if($rowtitle eq "No Alignments");
            $sheet->write($row,$col,$hashsummary{$sample}{'Not Map Uniquely'},$format{"normal"}) if($rowtitle eq "Not Map Uniquely");
            $sheet->write($row,$col,$hashsummary{$sample}{'totalC'},$format{"normal"}) if($rowtitle eq "Total number of C's analysed");
            $sheet->write($row,$col,$hashsummary{$sample}{'MCinCPG'},$format{"normal"}) if($rowtitle eq "Total methylated C's in CpG context");
            $sheet->write($row,$col,$hashsummary{$sample}{'MCinCHG'},$format{"normal"}) if($rowtitle eq "Total methylated C's in CHG context");
            $sheet->write($row,$col,$hashsummary{$sample}{'MCinCHH'},$format{"normal"}) if($rowtitle eq "Total methylated C's in CHH context");
            $sheet->write($row,$col,$hashsummary{$sample}{'unMCinCPG'},$format{"normal"}) if($rowtitle eq "Total unmethylated C's in CpG context");
            $sheet->write($row,$col,$hashsummary{$sample}{'unMCinCHG'},$format{"normal"}) if($rowtitle eq "Total unmethylated C's in CHG context");
            $sheet->write($row,$col,$hashsummary{$sample}{'unMCinCHH'},$format{"normal"}) if($rowtitle eq "Total unmethylated C's in CHH context");
            $sheet->write($row,$col,$hashsummary{$sample}{'percentinCPG'},$format{"normal"}) if($rowtitle eq "C methylated in CpG context");
            $sheet->write($row,$col,$hashsummary{$sample}{'percentinCHG'},$format{"normal"}) if($rowtitle eq "C methylated in CHG context");
            $sheet->write($row,$col,$hashsummary{$sample}{'percentinCHH'},$format{"normal"}) if($rowtitle eq "C methylated in CHH context");    
            $col++;         
        }
        $row++;
    }
    print "OK\n";
}
# Tag文件输出
sub Tag{
    my $config   = shift @_;
    my $workbook = shift @_;
    my $format   = shift @_;
    my $samples  = shift @_;
    print "Process Tag...";
    print "\n\t[Option] Output Tag Info?[y/n]";
    my $input=input();
    return if($input ne 'y');
    my $Report=$config->{'Report'};
    my %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.tag");# 检查文件
    my $isError=ShowError(\%hashCheckResult,"\n[WARN] Loss Tag, Sample:");# 错误信息输出
    # return if($isError);
    print "\tProcessing ...";
    # 读数据
    my %hashTag;
    foreach my $sample(@$samples){
        my $TagFile="$Report/status/$sample.tag";
        next if(checkFileSizeOK($TagFile)==0);
        open TAG,$TagFile;
        my $count=0;
        while(my $line=<TAG>){
            $count++;
            next if($count==1);
            $line=~ s/[\r\n]//g;
            my @datas=split /\t/,$line;
            my $chr=$datas[0];
            my $start=$datas[1];
            my $end=$datas[2];
            my $length=$datas[3];
            my $name=$datas[4];
            my $gc=$datas[5];
            my $cover=$datas[6];
            $hashTag{"$chr\t$start\t$end\t$length\t$name\t$gc"}{'OutputSortOrder'}=$count;# 排序用
            $hashTag{"$chr\t$start\t$end\t$length\t$name\t$gc"}{$sample}=$cover;
            
        }
        close TAG;
    } 
    # 输出
    my $sheet=$workbook->add_worksheet("Tag");
    $sheet->set_row(0, 60);
    my $row=0;
    my @titles=("chrom","start","end","length","name","\%gc",@$samples);
    foreach my $col(0..$#titles){
        $sheet->write($row,$col,$titles[$col],$format->{'title'});
    }
    $row++;
    foreach my $target(sort {$hashTag{$a}{'OutputSortOrder'}<=>$hashTag{$b}{'OutputSortOrder'}} keys %hashTag){
        my @datas=split /\t/,$target;
        map{push @datas,$hashTag{$target}{$_};}@$samples;
        foreach my $col(0..$#datas){
            $sheet->write($row,$col,$datas[$col],$format->{'normal'});
        }
        $row++;
    }
    print "OK\n";
}
# 覆盖情况统计
sub Metrics{
    my $config   = shift @_;
    my $workbook = shift @_;
    my $format   = shift @_;
    my $samples  = shift @_;
    print "Process Metrics...";
    my $Report=$config->{'Report'};
    my %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.metrics");# 检查文件
    my $isError=ShowError(\%hashCheckResult,"\n[WARN] Loss Metrics, Sample:");# 错误信息输出
    # return if($isError);
    my @MetricsTitles=("BAIT_TERRITORY","TOTAL_READS","PF_READS","PF_UQ_READS_ALIGNED","PCT_PF_UQ_READS_ALIGNED","PF_BASES_ALIGNED","PF_UQ_BASES_ALIGNED","ON_BAIT_BASES","NEAR_BAIT_BASES","OFF_BAIT_BASES","PCT_SELECTED_BASES","PCT_OFF_BAIT","FOLD_ENRICHMENT","PCT_USABLE_BASES_ON_BAIT","MEAN_BAIT_COVERAGE","PCT_TARGET_BASES_1X","PCT_TARGET_BASES_2X","PCT_TARGET_BASES_10X","PCT_TARGET_BASES_20X","PCT_TARGET_BASES_30X","PCT_TARGET_BASES_40X","PCT_TARGET_BASES_50X","PCT_TARGET_BASES_100X");
    my @QTitles=("Q20","Q30");
    my @PCRDupTitles=("PERCENT_DUPLICATION");
    my @titles=("Sample",@MetricsTitles);
    my %hashStat;
    my %hashCoverage;# 覆盖深度，绘图用
    # 读取Metrics文件提取数据
    foreach my $sample(@$samples){
        my $MetricsFile="$Report/status/$sample.metrics";
        next if(checkFileSizeOK($MetricsFile)==0);
        my %hashtmp;
        open METRICS,$MetricsFile;
        while(my $title=<METRICS>){
            next if($title=~ /^\#/);
            next if($title!~ /\w/);
            my $line=<METRICS>;
            my @titles=split /\t/,$title;
            my @datas=split /\t/,$line;
            foreach my $col(0..$#titles){
                $hashtmp{$titles[$col]}=$datas[$col];
            }
            last;
        }
        close METRICS;
        foreach my $title(@MetricsTitles){
            $hashStat{$sample}{$title}=$hashtmp{$title};
        }
        $hashStat{$sample}{'Sample'}=$sample;
        $hashCoverage{$sample}=$hashtmp{'MEAN_BAIT_COVERAGE'};
    }
    # 读取Q20Q30数据
    %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.quality");# 检查文件
    if((keys %hashCheckResult)==0){
        push @titles,@QTitles;
        foreach my $sample(@$samples){
            my $QualityFile="$Report/status/$sample.quality";
            open QUALITY,$QualityFile;
            my $c20=0;
            my $c30=0;
            my $all=0;
            while(my $line=<QUALITY>){
                $line=~s/[\r\n]//g;
                next if($line=~/^\#/ or $line!~/\w/);
                my @datas=split /\t/,$line;
                next if(@datas!=2);
                next if($datas[0]!~ /\d/ or $datas[1]!~ /\d/);
                $c20+=$datas[1] if($datas[0]>=20);
                $c30+=$datas[1] if($datas[0]>=30);
                $all+=$datas[1];
            }
            close QUALITY;
            my ($q20,$q30)=(0,0);
            if($all>0){
                $q20=$c20/$all;
                $q30=$c30/$all;
            }
            $hashStat{$sample}{'Q20'}=$q20;         
            $hashStat{$sample}{'Q30'}=$q30;         
        }
    }
    # 读取Duplicate数据
    %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####_sort_dup.metrics");# 检查文件
    if((keys %hashCheckResult)==0){
        push @titles,@PCRDupTitles;
        foreach my $sample(@$samples){
            my $DupFile="$Report/status/$sample"."_sort_dup.metrics";
            my %hashtmp;
            open DUP,$DupFile;
            while(my $title=<DUP>){
                next if($title=~ /^\#/);
                next if($title!~ /\w/);
                my $line=<DUP>;
                my @titles=split /\t/,$title;
                my @datas=split /\t/,$line;
                foreach my $col(0..$#titles){
                    $hashtmp{$titles[$col]}=$datas[$col];
                }
                last;
            }
            close DUP;
            foreach my $title(@PCRDupTitles){
                $hashStat{$sample}{$title}=$hashtmp{$title};
            }                   
        }
    }
    # 输出
    my $sheet=$workbook->add_worksheet("Metrics");
    $sheet->set_row(0, 60);
    my $row=0;
    foreach my $col(0..$#titles){
        $sheet->write($row,$col,$titles[$col],$format->{'title'});
    }
    $row++;
    foreach my $sample(@$samples){
        my @datas;
        foreach my $title(@titles){
            push @datas,$hashStat{$sample}{$title};
        }
        foreach my $col(0..$#datas){
            $sheet->write($row,$col,$datas[$col],$format->{'normal'});
        }
        $row++;          
    }
    PlotCoverage(\%hashCoverage,$Report);
    print "OK\n";
}
# 覆盖情况统计
sub GenomeMetrics{
    my $config   = shift @_;
    my $workbook = shift @_;
    my $format   = shift @_;
    my $samples  = shift @_;
    print "Process Metrics...";
    my $Report=$config->{'Report'};
    my %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.collect_wgs_metrics");# 检查文件
    my $isError=ShowError(\%hashCheckResult,"\n[WARN] Loss Metrics, Sample:");# 错误信息输出
    # return if($isError);
    my @MetricsTitles=("GENOME_TERRITORY","MEAN_COVERAGE","SD_COVERAGE","MEDIAN_COVERAGE","MAD_COVERAGE","PCT_EXC_MAPQ","PCT_EXC_DUPE","PCT_EXC_UNPAIRED","PCT_EXC_BASEQ","PCT_EXC_CAPPED","PCT_1X","PCT_5X","PCT_10X","PCT_15X","PCT_20X","PCT_25X","PCT_30X","PCT_40X","PCT_50X","PCT_60X","PCT_70X","PCT_80X","PCT_90X","PCT_100X");
    my @QTitles=("Q20","Q30");
    my @PCRDupTitles=("PERCENT_DUPLICATION");
    my @titles=("Sample",@MetricsTitles);
    my %hashStat;
    my %hashCoverage;# 覆盖深度，绘图用
    # 读取Metrics文件提取数据
    foreach my $sample(@$samples){
        my $MetricsFile="$Report/status/$sample.collect_wgs_metrics";
        next if(checkFileSizeOK($MetricsFile)==0);
        my %hashtmp;
        open METRICS,$MetricsFile;
        while(my $title=<METRICS>){
            next if($title=~ /^\#/);
            next if($title!~ /\w/);
            my $line=<METRICS>;
            my @titles=split /\t/,$title;
            my @datas=split /\t/,$line;
            foreach my $col(0..$#titles){
                $hashtmp{$titles[$col]}=$datas[$col];
            }
            last;
        }
        close METRICS;
        foreach my $title(@MetricsTitles){
            $hashStat{$sample}{$title}=$hashtmp{$title};
        }
        $hashStat{$sample}{'Sample'}=$sample;
        $hashCoverage{$sample}=$hashtmp{'MEAN_COVERAGE'};
    }
    # 读取Q20Q30数据
    %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.quality");# 检查文件
    if((keys %hashCheckResult)==0){
        push @titles,@QTitles;
        foreach my $sample(@$samples){
            my $QualityFile="$Report/status/$sample.quality";
            open QUALITY,$QualityFile;
            my $c20=0;
            my $c30=0;
            my $all=0;
            while(my $line=<QUALITY>){
                $line=~s/[\r\n]//g;
                next if($line=~/^\#/ or $line!~/\w/);
                my @datas=split /\t/,$line;
                next if(@datas!=2);
                next if($datas[0]!~ /\d/ or $datas[1]!~ /\d/);
                $c20+=$datas[1] if($datas[0]>=20);
                $c30+=$datas[1] if($datas[0]>=30);
                $all+=$datas[1];
            }
            close QUALITY;
            my ($q20,$q30)=(0,0);
            if($all>0){
                $q20=$c20/$all;
                $q30=$c30/$all;
            }
            $hashStat{$sample}{'Q20'}=$q20;         
            $hashStat{$sample}{'Q30'}=$q30;         
        }
    }
    # 读取Duplicate数据
    %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####_sort_dup.metrics");# 检查文件
    if((keys %hashCheckResult)==0){
        push @titles,@PCRDupTitles;
        foreach my $sample(@$samples){
            my $DupFile="$Report/status/$sample"."_sort_dup.metrics";
            my %hashtmp;
            open DUP,$DupFile;
            while(my $title=<DUP>){
                next if($title=~ /^\#/);
                next if($title!~ /\w/);
                my $line=<DUP>;
                my @titles=split /\t/,$title;
                my @datas=split /\t/,$line;
                foreach my $col(0..$#titles){
                    $hashtmp{$titles[$col]}=$datas[$col];
                }
                last;
            }
            close DUP;
            foreach my $title(@PCRDupTitles){
                $hashStat{$sample}{$title}=$hashtmp{$title};
            }                   
        }
    }
    # 输出
    my $sheet=$workbook->add_worksheet("Metrics");
    $sheet->set_row(0, 60);
    my $row=0;
    foreach my $col(0..$#titles){
        $sheet->write($row,$col,$titles[$col],$format->{'title'});
    }
    $row++;
    foreach my $sample(@$samples){
        my @datas;
        foreach my $title(@titles){
            push @datas,$hashStat{$sample}{$title};
        }
        foreach my $col(0..$#datas){
            $sheet->write($row,$col,$datas[$col],$format->{'normal'});
        }
        $row++;          
    }
    PlotCoverage(\%hashCoverage,$Report);
    print "OK\n";
}
# 覆盖深度绘图
sub PlotCoverage{
    my $hashdata =  shift @_;
    my $reportdir = shift @_;
    my $imgdir="$reportdir/document/01_QC/Image";
    mkdir $imgdir if(not -e $imgdir);
    my @samples;
    my @covers;
    foreach my $sample(sort {$hashdata->{$a} <=> $hashdata->{$b}} keys %$hashdata){
        push @samples,"'$sample'";
        push @covers,$hashdata->{$sample};
    }
    my $coveragelist=join ",",@covers;
    my $samplenamelist=join ",",@samples;
    my $samplenum=@samples;
    my $pic="$imgdir/mapping_coverage.png";
    my $R = Statistics::R->new();
    $R->startR;
    $R->send(qq` myColor=c("#aa4643","#89a54e","#71588f","#4198af","#db843d") `);
    $R->send(qq` png(\"$pic\",width=800,height=600) `) ;
    $R->send(qq` barplot(c($coveragelist),col=myColor[1],names.arg=c($samplenamelist),xlab="Samples",ylab="Percentage",lwd=2) `) if($samplenum<=5);
    $R->send(qq` barplot(c($coveragelist),col=myColor[1],xlab="Samples",ylab="Percentage",lty=1,lwd=2) `) if($samplenum>5);
    $R->send(qq` legend("topleft",c("MEAN_TARGET_COVERAGE"),cex=1,lwd=2,col=myColor[1],lty=c(1)) `);
    $R->send(qq` dev.off() `);  
}
# 目标片段区域覆盖统计
sub TargetCoverage{
    my $config   = shift @_;
    my $workbook = shift @_;
    my $format   = shift @_;
    my $samples  = shift @_;
    print "Process TargetCoverage...";
    my $Report=$config->{'Report'};
    my %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.status");# 检查文件
    my $isError=ShowError(\%hashCheckResult,"\n[WARN] Loss Status, Sample:");# 错误信息输出
    # return if($isError);
    my %hashTargetReads;
    foreach my $sample(@$samples){
        my $statusFile="$Report/status/$sample.status";
        next if(checkFileSizeOK($statusFile)==0);
        open FILE,$statusFile;
        while(my $line=<FILE>){
            $line=~ s/[\r\n]//g;
            next if($line=~/^\#/);
            my $target=(split /[\t\|]/,$line)[0];
            $hashTargetReads{$target}{$sample}+=(split /\t/,$line)[1];
        }
    }
    # 统计
    my %hashStat;
    my @targets=sort keys %hashTargetReads;
    my $targetnum=@targets;
    my @nX=(2,10,20,30);###记录nX数量
    foreach my $sample(@$samples){
        my $sum=0;
        my %hashNX;
        foreach my $target(@targets){
            my $reads=(exists($hashTargetReads{$target}{$sample}) and $hashTargetReads{$target}{$sample}=~/\d/) ? $hashTargetReads{$target}{$sample} : 0;
            $sum+=$reads;
            foreach my $X(@nX){# 统计大于X的数量
                if($reads>$X){
                    $hashNX{$X}+=1;
                }else{
                    $hashNX{$X}+=0;
                }
            }
        }
        $hashStat{$sample}{'Coverage'}=sprintf "%0.3f",$sum/$targetnum;
        map{$hashStat{$sample}{">$_"."X"}=sprintf "%0.3f",$hashNX{$_}/$targetnum;}@nX;
    }
    # 准备输出
    my $sheet=$workbook->add_worksheet("TargetCoverage");
    $sheet->set_row(0, 60);
    my @titles=("Sample","Coverage");
    map{push @titles,">$_"."X"}@nX;
    @titles=(@titles,@targets);
    my $row=0;
    for my $col(0..@titles-1){
        $sheet->write($row,$col,$titles[$col],$format->{"title"});
    }
    $row++;
    foreach my $sample(@$samples){
        my @datas;
        my @colors;
        foreach my $title(@titles){
            my $value="";
               $value=$sample if($title eq "Sample");
               $value=$hashStat{$sample}{$title} if(exists($hashStat{$sample}{$title}));
               $value=$hashTargetReads{$title}{$sample} if(exists($hashTargetReads{$title}{$sample}));
            my $thiscolor=(exists($hashTargetReads{$title}{$sample}) and $value<10) ? "gray" : "normal";
            push @datas,$value;
            push @colors,$thiscolor;
        }
        foreach my $col(0..$#datas){
            $sheet->write($row,$col,$datas[$col],$format->{$colors[$col]});
        }
        $row++;        
    }
    print "OK\n";
}
# 过滤后的Reads质控
sub QCCleanData{
    my $config   = shift @_;
    my $workbook = shift @_;
    my $format   = shift @_;
    my $samples  = shift @_;
    print "Process QC CleanData...";
    my $Report=$config->{'Report'};
    my %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.status");# 检查文件
    my $isError=ShowError(\%hashCheckResult,"\n[WARN] Loss Status, Sample:");# 错误信息输出
    # return if($isError);
    my %hashReads;
    my $flag_qc=0;
    foreach my $sample(@$samples){
        my %hashqc=ReadFastqc("$Report/fastqc_clean/$sample"."_fastqc/fastqc_data.txt");
        my $statusFile="$Report/status/$sample.status";
        next if(checkFileSizeOK($statusFile)==0);
        open FILE,$statusFile;
        while(my $line=<FILE>){
            $line=~ s/[\r\n]//g;
            if($line=~ /^\#/){
                my ($read,$readAll)=$line=~/(\d+)\t(\d+)/;
                $hashReads{$sample}{"Clean Reads"}+=$read;
                $hashReads{$sample}{"Raw Reads"}+=$readAll;
                $hashReads{$sample}{"Clean Reads Ratio"}=($hashReads{$sample}{"Raw Reads"}>0) ? sprintf("%.2f", $hashReads{$sample}{"Clean Reads"}/$hashReads{$sample}{"Raw Reads"}) : 0;
            }
        }
        $hashReads{$sample}{"Sample"}=$sample;
        $hashReads{$sample}{'Q20'}=$hashqc{'Q20'} if(exists($hashqc{'Q20'}));
        $hashReads{$sample}{'Q30'}=$hashqc{'Q30'} if(exists($hashqc{'Q30'}));
        $flag_qc++ if(exists($hashqc{'Q20'}));
        close FILE;
    }
    my $Reads=$workbook->add_worksheet("QC_CleanData");
    $Reads->set_row(0, 60);
    $Reads->set_column(1,3,15);
    my $row=0;
    my @titles=("Sample","Clean Reads","Raw Reads","Clean Reads Ratio");
    @titles=(@titles,"Q20","Q30") if($flag_qc>0);
    for my $col(0..@titles-1){
        $Reads->write($row,$col,$titles[$col],$format->{"title"});
    }
    $row++;
    foreach my $sample(@$samples){
        my @datas;
        foreach my $title(@titles){
            my $value=(exists($hashReads{$sample}{$title})) ? $hashReads{$sample}{$title} : "";
            push @datas,$value;
        }
        foreach my $col(0..$#datas){
            $Reads->write($row,$col,$datas[$col],$format->{"normal"});
        }
        $row++;
    }
    print "OK\n";
}
# 读取QC文件
sub ReadFastqc{
    my $qcfile=shift @_;
    my %hashqc;
    return %hashqc if(checkFileSizeOK($qcfile)==0);
    open FILE,$qcfile;
    my %hash=();
    my $flag="";
    while(my $line=<FILE>){
        $line=~ s/[\r\n]//g;
        $flag="" if($line=~ /\>\>END\_MODULE/);
        if($line=~ /Total\s+Sequences/){
            my $v=$line;$v=~ s/Total\s+Sequences//;$v=~ s/[\s\t]//g;
            $hashqc{"Total Reads"}=$v;
            next;
        }
        if($line=~ /Sequence\s+length/){
            my $v=$line;$v=~ s/Sequence\s+length//;$v=~ s/[\s\t]//g;
            $hashqc{"Reads length"}=$v;
            next;
        }
        if($flag eq "Quality"){
            my @split_line=split /[\s\t]+/,$line;
            $hash{"all"}+=$split_line[1];
            $hash{"q20"}+=$split_line[1] if($split_line[0]>=20);
            $hash{"q30"}+=$split_line[1] if($split_line[0]>=30);
        }
        $flag="Quality" if($line=~ /\#Quality/);
    }
    close FILE;
    $hash{"Q20"}=0;$hash{"Q20"}=$hash{"q20"}/$hash{"all"} if(exists $hash{"q20"} and exists $hash{"all"} and $hash{"all"}>0);
    $hash{"Q30"}=0;$hash{"Q30"}=$hash{"q30"}/$hash{"all"} if(exists $hash{"q30"} and exists $hash{"all"} and $hash{"all"}>0);    
    $hashqc{'Q20'}=sprintf "%0.2f",100*$hash{"Q20"};$hashqc{"Q20"}.="%";
    $hashqc{'Q30'}=sprintf "%0.2f",100*$hash{"Q30"};$hashqc{"Q30"}.="%";
    die 
    return %hashqc;
}
# RawData质控
sub fastqc{
    my $config   = shift @_;
    my $workbook = shift @_;
    my $format   = shift @_;
    my $samples  = shift @_;
    print "Process FastQC...";
    my $Report=$config->{'Report'};
    my %hashCheckResult=CheckFile($samples,"#####","$Report/fastqc/#####_R1_fastqc/fastqc_data.txt","$Report/fastqc/#####_R2_fastqc/fastqc_data.txt");# 检查文件
    my $isError=ShowError(\%hashCheckResult,"\n[WARN] Loss FastQC, Sample:");# 错误信息输出
    # return if($isError);
    my $sheet=$workbook->add_worksheet("QC_RawData");
    my @array=("Total Reads","Reads length","Q20","Q30","BaseCount(M)");
    $sheet->merge_range(0,0,1,0,"Sample",$format->{"title"});
    $sheet->merge_range(0,1,0,1+@array-1,"R1",$format->{"title"});
    $sheet->merge_range(0,1+@array,0,1+@array*2-1,"R2",$format->{"title"});
    for my $i(0..@array-1){
        $sheet->write(1,1+$i,$array[$i],$format->{"title"});
        $sheet->write(1,1+@array+$i,$array[$i],$format->{"title"});
    }
    my $row=2;
    foreach my $sample(@$samples){
        $sheet->write_string($row,0,$sample,$format->{"normal"});
        for my $i(1..2){
            my $fastqcData = "$Report/fastqc/$sample"."_R".$i."_fastqc/fastqc_data.txt";
            next if(checkFileSizeOK($fastqcData)==0);
            open FILE,$fastqcData;
            my %hash=();
            my %flag=();
            while(my $line=<FILE>){
                $line=~ s/[\r\n]//g;
                %flag=() if($line=~ /\>\>END\_MODULE/);
                if($line=~ /Total\s+Sequences/){
                    my $v=$line;$v=~ s/Total\s+Sequences//;$v=~ s/[\s\t]//g;
                    $hash{"Total Reads"}=$v;
                }
                if($line=~ /Sequence\s+length/){
                    my $v=$line;$v=~ s/Sequence\s+length//;$v=~ s/[\s\t]//g;
                    $hash{"Reads length"}=$v;
                }                
                
                if(exists $flag{"Quality"}){
                    my @split_line=split /[\s\t]+/,$line;
                    $hash{"all"}+=$split_line[1];
                    $hash{"q20"}+=$split_line[1] if($split_line[0]>=20);
                    $hash{"q30"}+=$split_line[1] if($split_line[0]>=30);
                }
                $flag{"Quality"}=1 if($line=~ /\#Quality/);
            }
            close FILE;
            my $fastqcStat = "$Report/fastqc/$sample"."_R".$i."_fastqc/fastqc_stat.txt";
            next if(checkFileSizeOK($fastqcStat)==0);
            open FILE,$fastqcStat;
            while(my $line=<FILE>){
                $line=~ s/[\r\n]//g;
                if($line=~ /^Base/){
                    my $v=$line;$v=~ s/Base//;$v=~ s/[\s\t]//g;
                    $hash{"BaseCount(M)"}=sprintf "%0.2f",$v/1000000;
                }
                if($line=~ /^Q20/){
                    my $v=$line;$v=~ s/Q20//;$v=~ s/[\s\t]//g;
                    $hash{"Q20-Base"}=$v;
                }
                if($line=~ /^Q30/){
                    my $v=$line;$v=~ s/Q30//;$v=~ s/[\s\t]//g;
                    $hash{"Q30-Base"}=$v;
                }
            }
            close FILE;
            $hash{"Q20"}=0;$hash{"Q20"}=$hash{"q20"}/$hash{"all"} if(exists $hash{"q20"} and exists $hash{"all"} and $hash{"all"}>0);
            $hash{"Q30"}=0;$hash{"Q30"}=$hash{"q30"}/$hash{"all"} if(exists $hash{"q30"} and exists $hash{"all"} and $hash{"all"}>0);
            $hash{"Q20"}=sprintf "%0.2f",100*$hash{"Q20"};$hash{"Q20"}.="%";
            $hash{"Q30"}=sprintf "%0.2f",100*$hash{"Q30"};$hash{"Q30"}.="%";            
            for my $j(0..@array-1){
                my $value=(exists($hash{$array[$j]})) ? $hash{$array[$j]} : "";
                $sheet->write($row,1+@array*($i-1)+$j,$value,$format->{"normal"});
            }
        }
        $row++;
    }
    print "OK\n";
}

# fastqc拷贝
sub prepare_QC_image{
    my $hashConfig = shift @_;
    my $samples    = shift @_;
    my $report_dir     = $hashConfig->{'Report'};
    my $document_dir   = "$report_dir/document/01_QC";
    my $image_dir      = "$document_dir/Image";
    my $fastqc_cp_dir  = "$image_dir/fastqc";  # 拷贝路径
    my $fastqc_dir     = "$report_dir/fastqc"; # 原图路径

    mkdir $document_dir if(not -e $document_dir);
    mkdir $image_dir if(not -e $image_dir);
    mkdir $fastqc_cp_dir if(not -e $fastqc_cp_dir);

    print "Copy fastqc image ... ";
    foreach my $sample(@$samples){
        my $R1QC       = "$fastqc_dir/$sample\_R1_fastqc/Images/per_base_quality.png";
        my $R1ATCG     = "$fastqc_dir/$sample\_R1_fastqc/Images/per_base_sequence_content.png";
        my $R2QC       = "$fastqc_dir/$sample\_R2_fastqc/Images/per_base_quality.png";
        my $R2ATCG     = "$fastqc_dir/$sample\_R2_fastqc/Images/per_base_sequence_content.png";
        my $Errorate   = "$fastqc_dir/$sample\_ErrorRate.png";
        my $InsertSize = "$report_dir/status/$sample.insert.pdf";
        system("cp $R1QC       $fastqc_cp_dir/$sample\_R1_per_base_quality.png") if(-e $R1QC);
        system("cp $R2QC       $fastqc_cp_dir/$sample\_R2_per_base_quality.png") if(-e $R2QC);
        system("cp $R1ATCG     $fastqc_cp_dir/$sample\_R1_per_base_sequence_content.png") if(-e $R1ATCG);
        system("cp $R2ATCG     $fastqc_cp_dir/$sample\_R2_per_base_sequence_content.png") if(-e $R2ATCG);
        system("cp $Errorate   $fastqc_cp_dir/") if(-e $Errorate);
        system("cp $InsertSize $fastqc_cp_dir/") if(-e $InsertSize);
    }
    print "OK\n";
}

# fastqc clean 拷贝
sub prepare_QC_Clean_image{
    my $hashConfig = shift @_;
    my $samples    = shift @_;
    my $report_dir     = $hashConfig->{'Report'};
    my $document_dir   = "$report_dir/document/01_QC";
    my $image_dir      = "$document_dir/Image";
    my $fastqc_cp_dir  = "$image_dir/fastqc_clean";  # 拷贝路径
    my $fastqc_dir     = "$report_dir/fastqc_clean"; # 原图路径

    return if(not -e $fastqc_dir);

    mkdir $document_dir if(not -e $document_dir);
    mkdir $image_dir if(not -e $image_dir);
    mkdir $fastqc_cp_dir if(not -e $fastqc_cp_dir);

    print "Copy fastqc clean image ... ";
    foreach my $sample(@$samples){
        my $R1QC       = "$fastqc_dir/$sample\_fastqc/Images/per_base_quality.png";
        my $R1ATCG     = "$fastqc_dir/$sample\_fastqc/Images/per_base_sequence_content.png";
        my $R2QC       = "$fastqc_dir/$sample\_fastqc/Images/per_base_quality.png";
        my $R2ATCG     = "$fastqc_dir/$sample\_fastqc/Images/per_base_sequence_content.png";
        my $Errorate   = "$fastqc_dir/$sample\_ErrorRate.png";
        system("cp $R1QC       $fastqc_cp_dir/$sample\_R1_per_base_quality.png") if(-e $R1QC);
        system("cp $R2QC       $fastqc_cp_dir/$sample\_R2_per_base_quality.png") if(-e $R2QC);
        system("cp $R1ATCG     $fastqc_cp_dir/$sample\_R1_per_base_sequence_content.png") if(-e $R1ATCG);
        system("cp $R2ATCG     $fastqc_cp_dir/$sample\_R2_per_base_sequence_content.png") if(-e $R2ATCG);
        system("cp $Errorate   $fastqc_cp_dir/") if(-e $Errorate);
    }    
    print "OK\n";
}

sub ShowError{
    my $hashCheckResult = shift @_;
    my $showinfo        = shift @_;
    my @errors=sort keys %$hashCheckResult;
    if(@errors>0){
        my $string=join ",",@errors;
        print "$showinfo $string\n";
        return 1;
    }
    return 0;
}

sub CheckFile{
    my $samples     = shift @_;
    my $replacemark = shift @_;
    my @files       = @_;
    my %hashCheckResult;# 有问题的样本
    foreach my $sample(@$samples){
        my $filenum=@files;
        my @checkFiles=@files;
        my $oknum=0;
        foreach my $file(@checkFiles){
            $file=~s/$replacemark/$sample/;# 还原文件真实路径
            next if(checkFileSizeOK($file)==0);
            $oknum++;
        }
        $hashCheckResult{$sample}++ if($filenum!=$oknum);
    }
    return %hashCheckResult;
}
sub checkFileSizeOK{
    my $file=shift @_;
    if(-e $file and -s $file>0){
        return 1;
    }else{
        return 0;
    }
}
# 从config里提取样本名称
sub GetSample{
    my $config=shift @_;
    my @types=@_;
    my @samples;
    foreach my $type(@types){
        next if(!exists($config->{$type}));
        foreach my $sample(split /,/,$config->{$type}){
            push @samples,$sample if($sample=~ /\w/);
        }
    }
    return @samples;    
}
# 读取配置文件
sub ReadConfig{
    my $file=shift @_;
    my %hash=();
    open FILE,$file or die "[ERR] Loss Config File,$file\n";
    while(my $line=<FILE>){
        $line=~ s/[\s\t\r\n]//g;
        my @split_line=split /\=/,$line;
        next if(@split_line!=2);
        if(exists $hash{$split_line[0]} and ($split_line[0] eq "case" or $split_line[0] eq "control") ){
            $hash{$split_line[0]}.=",".$split_line[1];
        }else{
            $hash{$split_line[0]}=$split_line[1];
        }
    }
    close FILE;
    return %hash;
}
sub input{
    my $input=<STDIN>;
    $input=~ s/[\r\n]//g;
    return $input;
}
sub ymd{
    my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
    $year+=1900;
    $mon+=1;
    return "$year-$mon-$day";
}

sub hms{
    my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
    return "$hour:$min:$sec";
}

# 获取文件路径名
sub get_dirname{
    my $path = shift @_;
    my $path_curf = File::Spec->rel2abs($path);
    my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
    return $dirs;
}

# Excel 表格格式
sub Format{
    my ($workbook)=@_;
    my %format=();
    $format{'title'} = $workbook->add_format();
    $format{'title'} ->set_align('center');
    $format{'title'} ->set_align('vcenter');
    $format{'title'} ->set_size(12);
    $format{'title'} ->set_font("Times New Roman");
    $format{'title'} ->set_border();
    $format{'title'} ->set_bg_color("yellow");
    $format{'title'} ->set_color("black");
    
    $format{'normal'} = $workbook->add_format();
    $format{'normal'} ->set_align('center');
    $format{'normal'} ->set_align('vcenter');
    $format{'normal'} ->set_size(12);
    $format{'normal'} ->set_font("Times New Roman");
    $format{'normal'} ->set_border();
    
    $format{'small'} = $workbook->add_format();
    $format{'small'} ->set_align('vcenter');
    $format{'small'} ->set_size(10);
    $format{'small'} ->set_font("Times New Roman");
    $format{'small'} ->set_border();
    
    $format{'seq'} = $workbook->add_format();
    $format{'seq'} ->set_align('vcenter');
    $format{'seq'} ->set_size(11);
    $format{'seq'} ->set_font("Courier New");
    $format{'seq'} ->set_border();
    
    $format{'left'} = $workbook->add_format();
    $format{'left'} ->set_align('vcenter');
    $format{'left'} ->set_size(12);
    $format{'left'} ->set_font("Times New Roman");
    $format{'left'} ->set_border();
    
    $format{'orange'} = $workbook->add_format();
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();

    $format{'gray'} = $workbook->add_format();
    $format{'gray'} ->set_align('vcenter');
    $format{'gray'} ->set_size(12);
    $format{'gray'} ->set_font("Times New Roman");
    $format{'gray'} ->set_bg_color("#BFBFBF");
    $format{'gray'} ->set_border();

    $format{'skyblue'} = $workbook->add_format();
    $format{'skyblue'} ->set_align('vcenter');
    $format{'skyblue'} ->set_size(12);
    $format{'skyblue'} ->set_font("Times New Roman");
    $format{'skyblue'} ->set_bg_color("#538ed5");
    $format{'skyblue'} ->set_border();

    $format{'bold'} = $workbook->add_format( bold => 1 );
    $format{'blue'} = $workbook->add_format( color => "#538ed5" );
    $format{'redbold'} = $workbook->add_format( color => "#ff0000", bold => 1, );
    $format{'italic'} = $workbook->add_format( italic => 1 );
    $format{'boldblue'} = $workbook->add_format( bold => 1, color => "#538ed5" );
    $format{'bolditalic'} = $workbook->add_format( bold => 1, italic => 1 );
    $format{'blueitalic'} = $workbook->add_format( color => "#538ed5", italic => 1 );
    $format{'boldblueitalic'} = $workbook->add_format( bold => 1, color => "#538ed5", italic => 1 );
    
    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    
    $format{'readme1'} = $workbook->add_format();
    $format{'readme1'}->set_align('center');
    $format{'readme1'}->set_align('vcenter');
    $format{'readme1'}->set_bold();
    $format{'readme1'}->set_size(14);
    $format{'readme1'}->set_font("Times New Roman");
    $format{'readme1'}->set_border();

    $format{'readme2'} = $workbook->add_format();
    $format{'readme2'}->set_align('vcenter');
    $format{'readme2'}->set_bold();
    $format{'readme2'}->set_size(14);
    $format{'readme2'}->set_font("Times New Roman");

    $format{'readme2tmp'} = $workbook->add_format();
    $format{'readme2tmp'}->set_right();

    $format{'readme3'} = $workbook->add_format();
    $format{'readme3'}->set_align('center');
    $format{'readme3'}->set_align('vcenter');
    $format{'readme3'}->set_bold();
    $format{'readme3'}->set_size(11);
    $format{'readme3'}->set_font("Times New Roman");
    $format{'readme3'}->set_border();

    $format{'readme4'} = $workbook->add_format();
    $format{'readme4'}->set_align('vcenter');
    $format{'readme4'}->set_bold();
    $format{'readme4'}->set_size(11);
    $format{'readme4'}->set_font("Times New Roman");
    $format{'readme4'}->set_border();

    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    return %format;
}