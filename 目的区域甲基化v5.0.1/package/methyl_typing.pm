package package::methyl_typing;
use strict;
use warnings;
$|=1;


sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start MethylTyping ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::methyl_typing($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in MethylTyping run para, please check\n" if($isOK == 0);

    my $report_dir = $hashConfig->{'Report'}; 
    my $output_dir = $hashConfig->{'Output'};
    my $methyl_dir = "$report_dir/methylation";
    package::utils::make_dir($methyl_dir);

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control', 'Ycontrol');
    my %hashCondition; 
    foreach my $sample(@samples)
    {
        my $baseFilter   = "$output_dir/$sample/$sample.blast.base.filter.gz";
        my $methylResult = "$methyl_dir/$sample.base";# 最终分析结果
        # 已完成该步骤
        if(package::utils::is_file_ok($methylResult) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$methylResult";
            next;            
        }
        # 原始数据没问题
        if(package::utils::is_file_ok($baseFilter)){
            $hashCondition{"Good2Run"}{$sample} = "$baseFilter";
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$baseFilter";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "MethylTyping", \%hashCondition); 
    die "Sample lost blast/fastq/reads.belong :" . (join ",", sort keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));

    # 执行比对
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0)
    {   
        my %hashAnalysisC;# 分析的C位点
        my %hashOtherC;# 分析的C以外的位点
        my %hashAllPos;# 片段上的每一个位点
        find_analysis_position($hashConfig, \%hashAnalysisC, \%hashOtherC, \%hashAllPos);# 确认分析的位点

        my $threshold = exists $hashPara->{"Process"}{"MethylTyping"} ? $hashPara->{"Process"}{"MethylTyping"} : 10;
           $threshold = $hashConfig->{"Process_MethylTyping"} if(exists $hashConfig->{"Process_MethylTyping"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::utils::process_bar_array($sample, \@sample_runs)});# 进度条
        foreach my $sample(@sample_runs)
        {
            $pm->start($sample) and next;
            run_methyl_typing($hashConfig, $sample, \%hashAnalysisC, \%hashOtherC, \%hashAllPos);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else
    {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}

sub run_methyl_typing{
    my $hashConfig    = shift @_;
    my $sample        = shift @_;
    my $hashAnalysisC = shift @_;
    my $hashOtherC    = shift @_;
    my $hashAllPos    = shift @_;
    my $report_dir    = $hashConfig->{"Report"};
    my $output_dir    = $hashConfig->{"Output"};
    my $methyl_dir    = "$report_dir/methylation";

    my %hashMethylCount;# 记录需要分析的C的情况
    my $otherC_Methy = 0;# 非要分析的C的甲基化计数
    my $otherC_Cover = 0;# 非要分析的C的甲基化+非甲基化计数
    my %hashHaplotype;# 单倍型
    my %hashSNV;# 记录每一个位置上的碱基数量
    my $baseFilter = "$output_dir/$sample/$sample.blast.base.filter.gz";# reads每个位点识别结果
    open BASEFILTER,"gzip -dc $baseFilter |";
    while(<BASEFILTER>)
    {
        $_=~s/[\r\n]//g;
        my ($readsName, $target, $identifyBaseInfo) = split /\t/,$_;
        my %hashBase;# 该reads在target上每一个位置的碱基
        foreach my $identifyBase(split /;/, $identifyBaseInfo)
        {
            next if($identifyBase!~/\w/);
            my ($position, $base, $quality) = split /,/, $identifyBase;
            $hashBase{$position} = $base;
        }
        # 识别
        my $haplotype="";# 当前reads获得的单倍型
        foreach my $position(sort { $a <=> $b } keys %hashBase)
        {
            my $base = $hashBase{$position};
            $hashSNV{$target}{$position}{$base}++;
            if(exists($hashAnalysisC->{$target}{$position}))# 要分析的C
            {
                $hashMethylCount{$target}{$position}{'Methy'}++ if($base =~ /c/i);# 当前位置甲基化计数
                $hashMethylCount{$target}{$position}{'Cover'}++ if($base =~ /[ct]/i);# 当前位置甲基化+非甲基化计数
                $haplotype=$haplotype.$base;
            }
            elsif(exists($hashOtherC->{$target}{$position}))# 非要分析的C
            {
                $otherC_Methy++ if($base =~ /c/i);
                $otherC_Cover++ if($base =~ /[ct]/i);
            }
        }
        my $haplotype_length = keys %{$hashAnalysisC->{$target}};# 片段单倍型应有的长度
        next if($haplotype =~ /[^ct]/i or length($haplotype) != $haplotype_length);# 单倍型不满足条件(单倍型包含非ct碱基或长度不等于应有长度)
        $hashHaplotype{$target}{$haplotype}++;# 记录该单倍型
    }
    close BASEFILTER;
    
    # 输出
    # 甲基化位点
    open METH_BASE, ">$methyl_dir/$sample.base";# 甲基化文件
    foreach my $target (sort keys %hashMethylCount){
        foreach my $position (sort{ $a <=> $b }keys %{$hashMethylCount{$target}}){
            my $methy = exists($hashMethylCount{$target}{$position}{'Methy'}) ? $hashMethylCount{$target}{$position}{'Methy'} : 0;# 当前位点甲基化数目
            my $cover = $hashMethylCount{$target}{$position}{'Cover'};# 当前位点读到的次数
            my $Ctype = $hashAnalysisC->{$target}{$position}{'Type'};# 当前位点的甲基化类型
            print METH_BASE "$target\t$position\t$Ctype\t$methy\t$cover\n"; 
        }
    }
    print METH_BASE "summary(C that beyond CpG)\tMethyC=$otherC_Methy\tAllPositionCount=$otherC_Cover\n";
    close METH_BASE;

    # 单倍型
    open METH_HAP, ">$methyl_dir/$sample.haplotype";# 单倍型文件
    foreach my $target (sort keys %hashHaplotype){
        foreach my $haplotype (sort{ $hashHaplotype{$target}{$b} <=> $hashHaplotype{$target}{$a} } keys %{$hashHaplotype{$target}}){
            print METH_HAP "$target\t$haplotype\t$hashHaplotype{$target}{$haplotype}\n";    
        }
    }
    close METH_HAP;

    # 每个位点碱基识别
    open SNV,">$methyl_dir/$sample.baseCount";
    print SNV "Target\tTarget_Strand\tPosition\tRefBase\tC2TBase\tCheckBase\tChr\tGenomePosition\n";
    foreach my $target (sort keys %hashSNV){
        foreach my $position (sort{ $a <=> $b } keys %{$hashSNV{$target}}){
            my @bases          = sort keys %{$hashSNV{$target}{$position}};
            my @baseCount      = map{"$_:$hashSNV{$target}{$position}{$_}"}@bases;
            my $string         = join ";",@baseCount;
            my $target_Strand  = $hashAllPos->{$target}{$position}{'Target_Strand'};
            my $Original       = $hashAllPos->{$target}{$position}{'Original'};
            my $C2T            = $hashAllPos->{$target}{$position}{'C2T'};
            my $Chr            = $hashAllPos->{$target}{$position}{'Chr'};
            my $GenomePosition = $hashAllPos->{$target}{$position}{'GenomePosition'};
            print SNV "$target\t$target_Strand\t$position\t$Original\t$C2T\t$string\t$Chr\t$GenomePosition\n";
        }
    }
    close SNV;
}

sub find_analysis_position{
    my $hashConfig    = shift @_;
    my $hashAnalysisC = shift @_;
    my $hashOtherC    = shift @_;
    my $hashAllPos    = shift @_;
    my $report_dir    = $hashConfig->{'Report'};
    my $analysis_point = "$report_dir/c.point.analysis.txt";
    my $other_point    = "$report_dir/c.point.other.txt";

    open ANALYSIS_C, $analysis_point;
    <ANALYSIS_C>;
    while(<ANALYSIS_C>){
        $_=~s/[\r\n]//g;
        my($target, $position, $type, @tmps) = split/\t/,$_;
        $hashAnalysisC->{$target}{$position}{'Type'} = $type;
    }
    close ANALYSIS_C;

    open OTHER_C, $other_point;
    while(<OTHER_C>){
        $_=~s/[\r\n]//g;
        my($target, $position) = split/\t/,$_;
        $hashOtherC->{$target}{$position}++;
    }
    close OTHER_C;

    # 片段上所有位点的位置
    %$hashAllPos = read_ref_seq($hashConfig->{"TargetFasta"}, "$report_dir/methylation_fasta_db/methylation.fa");# 参照序列每个位置的碱基
    my %hashAnnoPos = package::utils::region_anno($hashAllPos, $hashConfig);
    foreach my $target (sort keys %$hashAllPos){
        foreach my $position (sort{ $a <=> $b } keys %{$hashAllPos->{$target}}){
            $hashAllPos->{$target}{$position}{'Chr'}            = $hashAnnoPos{$target}{$position}{'Chr'};
            $hashAllPos->{$target}{$position}{'GenomePosition'} = $hashAnnoPos{$target}{$position}{'GenomePosition'};
            $hashAllPos->{$target}{$position}{'Target_Strand'}  = $hashAnnoPos{$target}{$position}{'Target_Strand'};
        }
    }
}

sub read_ref_seq{
    my $originalseqfile = shift @_;
    my $C2Tseqfile      = shift @_;
    my %hashoriginal    = package::utils::read_fasta($originalseqfile);
    my %hashC2T         = package::utils::read_fasta($C2Tseqfile);
    my %hashAllPos;
    foreach my $target(keys %hashoriginal){
        my @originalbps = split //, $hashoriginal{$target};
        my @C2Tseqbps   = split //, $hashC2T{$target};   
        foreach my $position(1..length($hashoriginal{$target})){
            $hashAllPos{$target}{$position}{'Original'} = $originalbps[$position-1];
            $hashAllPos{$target}{$position}{'C2T'}      = $C2Tseqbps[$position-1];
        }       
    }
    return %hashAllPos;
}

1