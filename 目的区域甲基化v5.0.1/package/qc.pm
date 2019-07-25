package package::qc;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start QC ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::qc($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in QC run para, please check\n" if($isOK == 0);

    my $report_dir = $hashConfig->{'Report'}; 
    my $fastq_dir  = $hashConfig->{'Fastq'}; 
    my $fastqc_dir = "$report_dir/fastqc";
    package::utils::make_dir($fastqc_dir);

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control');
    my %hashCondition; 
    foreach my $sample(@samples)
    {
        my $qc_R1 = "$fastqc_dir/$sample\_R1_fastqc";
        my $qc_R2 = "$fastqc_dir/$sample\_R2_fastqc";
        my $fastq_R1 = "$fastq_dir/$sample\_R1.fastq.gz";
        my $fastq_R2 = "$fastq_dir/$sample\_R2.fastq.gz";
        # 已完成该步骤
        if(package::utils::is_dir_ok($qc_R1) and package::utils::is_dir_ok($qc_R2) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$qc_R1 $qc_R1";
            next;            
        }
        # 原始数据没问题
        if(package::utils::is_file_ok($fastq_R1) and package::utils::is_file_ok($fastq_R2)){
            $hashCondition{"Good2Run"}{$sample} = "$fastq_R1 $fastq_R2";
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$fastq_R1 $fastq_R2";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "QC", \%hashCondition); 
    die "Sample lost fastq :" . (join ",", sort keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));

    # 执行QC
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0)
    {
        my $threshold = exists $hashPara->{"Process"}{"QC"} ? $hashPara->{"Process"}{"QC"} : 10;
           $threshold = $hashConfig->{"Process_QC"} if(exists $hashConfig->{"Process_QC"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::utils::process_bar_array($sample, \@sample_runs)});# 进度条
        foreach my $sample(@sample_runs)
        {
            $pm->start($sample) and next;
            run_qc($hashPara, $hashConfig, $sample);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else
    {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}

sub run_qc{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $sample     = shift @_;
    my $report_dir = $hashConfig->{'Report'}; 
    my $fastq_dir  = $hashConfig->{'Fastq'}; 
    my $fastqc_dir = "$report_dir/fastqc";
    my $FastQC    = $hashPara->{"Soft"}{"FastQC"};
    my $FastqStat = $hashPara->{"Soft"}{"FastqStat"};
    my $tmp_dir   = $hashPara->{'Soft'}{'Tmp'};
    my $fastq_R1 = "$fastq_dir/$sample\_R1.fastq.gz";
    my $fastq_R2 = "$fastq_dir/$sample\_R2.fastq.gz"; 

    system("$FastQC --extract -q -o $fastqc_dir -d $tmp_dir $fastq_R1 $fastq_R2");
    system("rm $fastqc_dir/$sample\_R1_fastqc.zip $fastqc_dir/$sample\_R2_fastqc.zip");
    system "rm $fastqc_dir/$sample\_R1_fastqc.html $fastqc_dir/$sample\_R2_fastqc.html";
    
    my $R1QCData="$fastqc_dir/$sample\_R1_fastqc/fastqc_data.txt";
    my $R2QCData="$fastqc_dir/$sample\_R2_fastqc/fastqc_data.txt";
    QCplot($hashPara, $fastqc_dir, $sample, $R1QCData, $R2QCData) if(package::utils::is_file_ok($R1QCData) and package::utils::is_file_ok($R2QCData));
    
    # Q20,Q30统计
    QC2030($fastq_R1, "$fastqc_dir/$sample\_R1_fastqc/fastqc_stat.txt", $FastqStat);
    QC2030($fastq_R2, "$fastqc_dir/$sample\_R2_fastqc/fastqc_stat.txt", $FastqStat);
}

sub QC2030{
    my $fastq     = shift @_;
    my $qc_file   = shift @_;
    my $FastqStat = shift @_;
    my $info = `$FastqStat $fastq`;
    my ($BaseCount) = $info =~ /Base=(\d+)/;
    my ($BaseQ20Count) = $info =~ /BaseQ20=(\d+)/;
    my ($BaseQ30Count) = $info =~ /BaseQ30=(\d+)/;
    my $Q20 = sprintf "%0.2f", (100 * $BaseQ20Count) / $BaseCount;
    my $Q30 = sprintf "%0.2f", (100 * $BaseQ30Count) / $BaseCount;
    # $BaseCount=$BaseCount/1000000;# 单位为M
    open QC,">$qc_file";
    print QC "Base\t$BaseCount\n";
    print QC "Q20\t$Q20\n";
    print QC "Q30\t$Q30\n";
    close QC;
}

sub QCplot{
    my $hashPara   = shift @_;
    my $fastqc_dir = shift @_;
    my $sample     = shift @_;
    my $R1QCData   = shift @_;
    my $R2QCData   = shift @_;

    my $R    = $hashPara->{'Soft'}{'R'};
    my $RLib = $hashPara->{'Soft'}{'RLib'};

    my %hashquality;
    my %hashAGCTN;
    read_quality($R1QCData,\%hashquality,\%hashAGCTN,'R1');
    read_quality($R2QCData,\%hashquality,\%hashAGCTN,'R2');
    PlotError($R, $RLib, \%hashquality, "$fastqc_dir/$sample\_ErrorRate.png");   
}

sub PlotError{
    my $Rbin        = shift @_;
    my $RLib        = shift @_;
    my $hashquality = shift @_;
    my $out_png     = shift @_;
    my $count=0;
    my @positionsR1=sort {$a<=>$b} keys %{$hashquality->{'R1'}};
    my @positionsR2=sort {$a<=>$b} keys %{$hashquality->{'R2'}};
    my @row=();
    my @values=();
    foreach(@positionsR1){
        $count++;
        push @row,$count;push @values,$hashquality->{'R1'}{$_};
    }
    foreach(@positionsR2){
        $count++;
        push @row,$count;push @values,$hashquality->{'R2'}{$_};
    }
    my $rowstr=join ",",@row;
    my $valuestr=join ",",@values;
    #画图
 
    my $R = Statistics::R->new(bin => $Rbin);
    $R->startR;
    $R->send(qq` .libPaths("$RLib") \n`) ;
    $R->send(qq` library(ggplot2) \n`) ;
    $R->send(qq` Error=data.frame(Position = c($rowstr), ErrorRate = c($valuestr))\n`) ;
    $R->send(qq` png("$out_png",width=1600,height=1200)`) ;
    $R->send(qq` color <- "#32CD32"\n`) ; 
    $R->send(qq` fill_color <- "#00FF7F"\n`) ;
    $R->send(qq` p <- ggplot(data=Error,aes(x=Position,y=ErrorRate))\n`) ;
    $R->send(qq` p + geom_bar(fill=fill_color,stat="identity",color=color) + xlab("Position along reads") + ylab("\%Error rate")\n`) ;
    $R->send(qq` last_plot()+ theme(axis.title.x = element_text(size=20,colour = "black",face = "bold"),axis.title.y = element_text(size=20,colour = "black",face = "bold"),axis.text.y = element_text(size=20,colour = "black",face = "bold"),axis.text.x = element_text(size=20,colour = "black",face = "bold"))\n`) ;
    $R->send(qq` dev.off()`) ;
    $R->stop(); 
}

sub read_quality{
    my $QCData      = shift @_;
    my $hashquality = shift @_;
    my $hashAGCTN   = shift @_;
    my $type        = shift @_;
    open IN,$QCData;
    my $model = '';
    my @ATCGheads=();
    while(<IN>){
        $_=~s/[\r\n]//g;
        next if($_!~/\w/ or $_=~/^##/);
        if($_=~/^>>END_MODULE/)
        {
            $model = '';
            next;
        }
        if($_=~/^>>/)
        {
            my ($info, $tmp) = split /\t/, $_;
            $model = $info;
            $model =~s/^>>//;
            next;
        }
        if($model eq 'Per base sequence quality')
        {
            next if($_=~/^#/);
            my ($base, $quality, $tmp)=split /\t/, $_, 3;
            my @bases = split /-/, $base;
            my $errorrate = 10**(-$quality/10);#错误率
            foreach my $position($bases[0]..$bases[$#bases])
            {
                $hashquality->{$type}{$position}=100*$errorrate;
            }
        }

    }
    close IN;
}
1