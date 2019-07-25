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
my $report_dir   = $config{"Report"};
my $document_dir = "$report_dir/document";
my $qc_dir       = "$document_dir/1_QC";
mkdir $document_dir if(not -e $document_dir);
mkdir $qc_dir if(not -e $qc_dir);
print "$qc_dir/Data_statistics.xlsx\n";
my $workbook=Excel::Writer::XLSX->new("$qc_dir/Data_statistics.xlsx");
my %format=Format($workbook);
fastqc(\%config,$workbook,\%format,\@samples);# 下机数据质控
Metrics(\%config,$workbook,\%format,\@samples);# 目标区域数据量统计,Picard
Tag(\%config,$workbook,\%format,\@samples);# tag文件输出
ReadMe($workbook,\%format,"readme_tongji.txt");
########################## FUNCTION ########################
sub ReadMe{
    my ($workbook,$format,$file)=@_;
    return if(checkFileSizeOK($file)==0);
    my $sheet=$workbook->add_worksheet("Read Me");
    $sheet->set_row(0, 65);
    $sheet->set_column('A:A', 35);
    $sheet->set_column('B:B', 25);
    $sheet->set_column('C:C', 110);
    my $row=0;
    open FILE,$file;
    while(<FILE>){
        $_=~ s/\^/\n/g;
        my @split_line=split /\t/,$_;
        my $col=0;
        foreach(@split_line)
        {
            my $text=decode("gb2312",$_);
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
    my @titles=("chrom","start","end","length","name","\%gc");
    foreach my $col(0..$#titles){
        $sheet->write($row,$col,$titles[$col],$format->{'title'});
    }
    foreach my $col (0..$#$samples)
    {
        $sheet->write_string($row,$col+$#titles+1,$$samples[$col],$format->{'title'});      
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
    my %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.metrics");# 检查文件sample.metrics
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
    %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####.quality.metrics");# 检查文件
    if((keys %hashCheckResult)==0){
        push @titles,@QTitles;
        foreach my $sample(@$samples){
            my %hashtmp = ();
            my $QualityFile="$Report/status/$sample.quality.metrics";
            next if(checkFileSizeOK($QualityFile)==0);
            open QUALITY,$QualityFile;
            my $c20=0;
            my $c30=0;
            my $all=0;
            while(my $title=<QUALITY>){
                  next if($title=~ /^\#/);
                  next if($title!~ /\w/);
                  my $line=<QUALITY>;
                  my @titles=split /\t/,$title;
                  my @datas=split /\t/,$line;
                  foreach my $col(0..$#titles){
                      $hashtmp{$titles[$col]}=$datas[$col];
                  }
            }
            close QUALITY;
            my $sum = $hashtmp{"TOTAL_BASES"};
            foreach my $title(@QTitles){
                my $newtitle = "$title"."_BASES";
                $hashStat{$sample}{$title}= $hashtmp{$newtitle} / $sum;
            }        
        }
    }
    # 读取Duplicate数据
    %hashCheckResult=CheckFile($samples,"#####","$Report/status/#####_sort_dup.metrics");# 检查文件
    if((keys %hashCheckResult)==0){
        push @titles,@PCRDupTitles;
        foreach my $sample(@$samples){
            my $DupFile="$Report/status/$sample"."_sort_dup.metrics";
            next if(checkFileSizeOK($DupFile)==0);
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
            my $value=(exists($hashStat{$sample}{$title})) ? $hashStat{$sample}{$title} : "";
            push @datas,$value;            
        }
        $sheet->write_string($row,0,$sample,$format->{'normal'});
        foreach my $col(1..$#datas){
            $sheet->write($row,$col,$datas[$col],$format->{'normal'});
        }
        $row++;          
    }
    PlotCoverage(\%hashCoverage,$Report);
    # 行业内部质控标准检查
    my $string1 = "Warnings : PERCENT_DUPLICATION > 0.25 : ";
    my $string2 = "Warnings : 10X/20X Coverage < 0.95 : ";
    my $string3 = "Warnings : PCT_USABLE_BASES_ON_BAIT < 0.5 : ";
    my $string4 = "Warnings : PCT_PF_UQ_READS_ALIGNED < 0.95 : ";
    my $flag = 0;
    foreach my $title(@titles){
        foreach my $sample(@$samples){
            my $value=(exists($hashStat{$sample}{$title})) ? $hashStat{$sample}{$title} : "";
            if($title eq "PERCENT_DUPLICATION" and $value=~/\d/ and $value > 0.25){
                $string1.="$sample,";$flag=1;
            }
            if(($title eq "PCT_TARGET_BASES_10X" or $title eq "PCT_TARGET_BASES_20X") and $value=~/\d/ and $value<0.95){
                $string2.="$sample,";$flag=1;
            }
            if($title eq "PCT_USABLE_BASES_ON_BAIT" and $value=~/\d/ and $value<0.5){
                $string3.="$sample,";$flag=1;
            }
            if($title eq "PCT_PF_UQ_READS_ALIGNED" and $value=~/\d/ and $value<0.95){
                $string4.="$sample,";$flag=1;
            }
        }
    }
    print "$string1\n$string2\n$string3\n$string4\n" if($flag==1);
    print "OK\n";
}

# 覆盖深度绘图
sub PlotCoverage{
    my $hashdata  = shift @_;
    my $reportdir = shift @_;
    my $imgdir    = "$reportdir/document/1_QC/fastqc";
    mkdir $imgdir if(not -e $imgdir);
    my @samples;
    my @covers;
    foreach my $sample(sort {$hashdata->{$a} <=> $hashdata->{$b}} keys %$hashdata){
        push @samples,"'$sample'";
        push @covers,sprintf "%0.1f", $hashdata->{$sample};
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
    return %hashqc;
}
# RawData质控
#fastqc(\%config,$workbook,\%format,\@samples);
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
            my $fastqcData="$Report/fastqc/$sample"."_R".$i."_fastqc/fastqc_data.txt";
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
                if($line=~ /^BaseCount/){
                    my $v=$line;$v=~ s/BaseCount//;$v=~ s/[\s\t]//g;
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
                if(exists $flag{"Quality"}){
                    my @split_line=split /[\s\t]+/,$line;
                    $hash{"all"}+=$split_line[1];
                    $hash{"q20"}+=$split_line[1] if($split_line[0]>=20);
                    $hash{"q30"}+=$split_line[1] if($split_line[0]>=30);
                }
                $flag{"Quality"}=1 if($line=~ /\#Quality/);
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
sub fastqcBase{
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
    my @array=("Total Reads","Reads length","Q20","Q30","BaseCount(M)","Q20-Base","Q30-Base");
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
            my $fastqcData="$Report/fastqc/$sample"."_R".$i."_fastqc/fastqc_data.txt";
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
                if($line=~ /^BaseCount/){
                    my $v=$line;$v=~ s/BaseCount//;$v=~ s/[\s\t]//g;
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
                if(exists $flag{"Quality"}){
                    my @split_line=split /[\s\t]+/,$line;
                    $hash{"all"}+=$split_line[1];
                    $hash{"q20"}+=$split_line[1] if($split_line[0]>=20);
                    $hash{"q30"}+=$split_line[1] if($split_line[0]>=30);
                }
                $flag{"Quality"}=1 if($line=~ /\#Quality/);
            }
            close FILE;
            $hash{"Q20"}=0;$hash{"Q20"}=$hash{"q20"}/$hash{"all"} if(exists $hash{"q20"} and exists $hash{"all"} and $hash{"all"}>0);
            $hash{"Q30"}=0;$hash{"Q30"}=$hash{"q30"}/$hash{"all"} if(exists $hash{"q30"} and exists $hash{"all"} and $hash{"all"}>0);
            $hash{"Q20"}=sprintf "%0.2f",100*$hash{"Q20"};$hash{"Q20"}.="%";
            $hash{"Q30"}=sprintf "%0.2f",100*$hash{"Q30"};$hash{"Q30"}.="%";            
            for my $j(0..@array-1){
                $sheet->write($row,1+@array*($i-1)+$j,$hash{$array[$j]},$format->{"normal"});
            }
        }
        $row++;
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
# 没有去掉#
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
            #'case'=>'case,15A,18B'
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