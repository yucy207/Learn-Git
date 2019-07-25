package package::mapping_bwa;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start BWA Mapping ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::mapping_bwa($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in BWA Mapping run para, please check\n" if($isOK == 0);

    my $report_dir       = $hashConfig->{'Report'}; 
    my $output_dir       = $hashConfig->{'Output'}; 
    my $status_dir       = "$report_dir/status";
    package::utils::make_dir($status_dir);

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control', 'Ycontrol');
    my %hashCondition; 
    foreach my $sample(@samples)
    {
        my $bam_index = "$output_dir/$sample/$sample\_final.bai";
        my $fastq     = "$output_dir/$sample/$sample.fastq";
        # 已完成该步骤
        if(package::utils::is_file_ok($bam_index) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$bam_index";
            next;            
        }
        # 原始数据没问题
        if(package::utils::is_file_ok($fastq)){
            $hashCondition{"Good2Run"}{$sample} = "$fastq";
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$fastq";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "Mapping_BWA", \%hashCondition); 
    die "Sample lost fastq :" . (join ",", sort keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));

    # 执行比对
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0)
    {   
        prepareGenome($hashPara, $hashConfig);# 建立基因组index
        my $threshold = exists $hashPara->{"Process"}{"MappingBWA"} ? $hashPara->{"Process"}{"MappingBWA"} : 10;
           $threshold = $hashConfig->{"Process_MappingBWA"} if(exists $hashConfig->{"Process_MappingBWA"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::utils::process_bar_array($sample, \@sample_runs)});# 进度条
        foreach my $sample(@sample_runs)
        {
            $pm->start($sample) and next;
            run_mapping($hashPara, $hashConfig, $sample);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else
    {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}

sub run_mapping{
    my $hashPara    = shift @_;
    my $hashConfig  = shift @_;
    my $sample      = shift @_;
    my $Java        = $hashPara->{'Soft'}{'Java'};
    my $BWA         = $hashPara->{'Soft'}{'BWA'};
    my $sambamba    = $hashPara->{"Soft"}{"sambamba"};
    my $cpulimit    = $hashPara->{"Soft"}{"cpulimit"};
    my $samtools    = $hashPara->{'Soft'}{'SamTools'};
    my $Picard      = $hashPara->{'Soft'}{'Picard'};
    my $Tmp         = $hashPara->{"Soft"}{"Tmp"};
    my $report_dir  = $hashConfig->{'Report'};
    my $output_dir  = "$hashConfig->{'Output'}/$sample";
    my $status_dir  = "$report_dir/status";
    my $genome_dir  = "$report_dir/TargetFastaGenome";
    my $genome_file = "$genome_dir/seq.fa";     # 参考基因组
    my $target_bed  = "$genome_dir/target.bed"; # 目标区域
    # mapping
    my $bam          = "$output_dir/$sample.bam"; 
    my $mapped_fastq = "$output_dir/$sample.fastq"; 
    # system "$BWA mem -t 4 -k 12 -M -T 30 -B 8 -O 12 -L 10  -U 18 -R '\@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:illumina\tPU:barcode' $genome_file $mapped_fastq | $samtools view -bS -h - > $bam";
    system "$BWA mem -t 4 -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:illumina\\tPU:barcode' $genome_file $mapped_fastq | $samtools view -bS -h - > $bam";
    # sort
    my $bam_final     = "$output_dir/$sample\_final.bam"; 
    my $bam_final_idx = "$output_dir/$sample\_final.bai";
    system "$sambamba sort $bam -m 6G -t 4 -o $bam_final --tmpdir $Tmp";
    system "mv $bam_final\.bai $bam_final_idx";
    # system "$cpulimit -l 400 $Java -jar -Xmx6g $Picard SortSam INPUT=$bam OUTPUT=$bam_final SORT_ORDER=coordinate  VALIDATION_STRINGENCY=LENIENT TMP_DIR=$Tmp";
    # system "$cpulimit -l 400 $Java -jar -Xmx6g $Picard BuildBamIndex INPUT=$bam_final OUTPUT=$bam_final_idx  VALIDATION_STRINGENCY=LENIENT  TMP_DIR=$Tmp";
    # HsMetrics
    my $metrics     = "$status_dir/$sample.metrics"; 
    my $tag         = "$status_dir/$sample.tag";
    my $quality     = "$status_dir/$sample.quality";
    my $quality_pdf = "$status_dir/$sample.quality.pdf";
    system "$cpulimit -l 400 $Java -jar -Xmx4g $Picard CollectHsMetrics INPUT=$bam_final OUTPUT=$metrics R=$genome_file BI=$target_bed TI=$target_bed PER_TARGET_COVERAGE=$tag MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 VALIDATION_STRINGENCY=LENIENT TMP_DIR=$Tmp";
    system "$cpulimit -l 400 $Java -jar -Xmx4g $Picard QualityScoreDistribution INPUT=$bam_final OUTPUT=$quality R=$genome_file PF_READS_ONLY=true CHART_OUTPUT=$quality_pdf VALIDATION_STRINGENCY=LENIENT TMP_DIR=$Tmp";
}

sub prepareGenome{
    my $hashPara     = shift @_;
    my $hashConfig   = shift @_;
    my $report_dir   = $hashConfig->{'Report'};
    my $seq_region   = "$report_dir/seq.region";
    my $genome_dir   = "$report_dir/TargetFastaGenome";
    my $genome_file  = "$genome_dir/seq.fa";  # 参考基因组  
    mkdir $genome_dir if(not -e $genome_dir);
    system("cp $report_dir/methylation_fasta_db/methylation.fa $genome_file");# 把seq.fa单独拷贝至其他目录
    #####
    # index
    #####
    my $Java      = $hashPara->{'Soft'}{'Java'};
    my $BWA       = $hashPara->{'Soft'}{'BWA'};
    my $samtools  = $hashPara->{'Soft'}{'SamTools'};
    my $Picard    = $hashPara->{'Soft'}{'Picard'};
    my $dict      = $genome_file;
       $dict      =~ s/fa$/dict/;
    system("$BWA index $genome_file");# samtools index
    system("$samtools faidx $genome_file");# samtools index
    system("rm $dict") if(-e $dict);
    system("$Java -jar $Picard CreateSequenceDictionary R=$genome_file O=$dict");# picard dict
    #####
    # 生成bed文件,用于统计
    #####
    my %hashSeqRegion      = package::utils::read_seq_region($seq_region);# 读取引物序列
    my $target_bed         = "$genome_dir/target.bed";
    my $target_realign_bed = "$genome_dir/targetRealign.bed";
    system("cp $dict $target_bed");
    open BED,">>$target_bed";# 追加模式
    open REALIGNBED,">$target_realign_bed";
    open DICT, $dict;
    while(<DICT>){
        $_=~s/[\r\n]//g;
        next if($_=~/^\@HD/);
        my ($sq, $chr, $length, $tmp) = split /\t/,$_;
        $chr=~s/^SN://;
        $length=~s/^LN://;
        my $primerFSeq = $hashSeqRegion{$chr}{'PrimerF'};
        my $primerRSeq = $hashSeqRegion{$chr}{'PrimerR'};
        my $start = length($primerFSeq)+1;# 去掉开头结尾的引物区域
        my $end   = $length-length($primerRSeq);
        print BED "$chr\t$start\t$end\t+\t$chr\n";
        print REALIGNBED "$chr\t$start\t$end\t+\t$chr\n";
    } 
    close BED;
    close REALIGNBED; 
}

1