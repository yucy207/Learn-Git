package package::mapping;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start Mapping ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::mapping($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in Mapping run para, please check\n" if($isOK == 0);

    my $report_dir = $hashConfig->{'Report'}; 
    my $fastq_dir  = $hashConfig->{'Fastq'}; 
    my $output_dir = $hashConfig->{'Output'}; 
    my $status_dir = "$report_dir/status";
    package::utils::make_dir($status_dir);

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control');
    my %hashCondition; 
    foreach my $sample(@samples)
    {
        my $bam_final = "$output_dir/$sample/$sample\_final.bam";
        my $bai_final = "$output_dir/$sample/$sample\_final.bai";
        my $fastq_R1 = "$fastq_dir/$sample\_R1.fastq.gz";
        my $fastq_R2 = "$fastq_dir/$sample\_R2.fastq.gz";
        # 已完成该步骤
        if(package::utils::is_file_ok($bam_final) and package::utils::is_file_ok($bai_final) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$bam_final $bai_final";
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
    package::utils::write_log("$report_dir/run.log", "Mapping", \%hashCondition); 
    die "Sample lost fastq :" . (join ",", sort keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));
    
    # 执行比对
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0)
    {
        my $threshold = exists $hashPara->{"Process"}{"Mapping"} ? $hashPara->{"Process"}{"Mapping"} : 10;
           $threshold = $hashConfig->{"Process_Mapping"} if(exists $hashConfig->{"Process_Mapping"});
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
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $sample     = shift @_;
    my $report_dir = $hashConfig->{'Report'}; 
    my $fastq_dir  = $hashConfig->{'Fastq'}; 
    my $output_dir = "$hashConfig->{'Output'}/$sample"; 
    my $status_dir = "$report_dir/status";
    my $Species    = $hashConfig->{"Species"};
    my $TargetBed  = $hashConfig->{"TargetBed"};
    my $RealignBed = $hashConfig->{"RealignBed"};
    package::utils::make_dir($output_dir);

    my $tmp_dir   = $hashPara->{'Soft'}{'Tmp'};
    my $BWA       = $hashPara->{"Soft"}{"BWA"};
    my $SamTools  = $hashPara->{"Soft"}{"SamTools"};
    my $sambamba  = $hashPara->{"Soft"}{"sambamba"};
    my $cpulimit  = $hashPara->{"Soft"}{"cpulimit"};
    my $Picard    = $hashPara->{"Soft"}{"Picard"};
    my $Java      = $hashPara->{"Soft"}{"Java"};
    my $GATK3     = $hashPara->{"Soft"}{"GATK3"};
    my $GATK4_Loc = $hashPara->{"Soft"}{"GATK4_Loc"};
    my $Genome    = $hashPara->{$Species}{"Genome"};
    my $bwa_idx   = $hashPara->{$Species}{"bwa_idx"};
    my $knownDBsnp     = (exists $hashPara->{$Species}{"DBsnp"})     ? "-known $hashPara->{$Species}{'DBsnp'}"     :  "";
    my $knownInDel     = (exists $hashPara->{$Species}{"InDel"})     ? "-known $hashPara->{$Species}{'InDel'}"     :  "";
    my $knownInDelGold = (exists $hashPara->{$Species}{"InDelGold"}) ? "-known $hashPara->{$Species}{'InDelGold'}" :  "";

    my $knownSitesDBsnp     = (exists $hashPara->{$Species}{"DBsnp"})     ? "--known-sites $hashPara->{$Species}{'DBsnp'}"     :  "";
    my $knownSitesInDel     = (exists $hashPara->{$Species}{"InDel"})     ? "--known-sites $hashPara->{$Species}{'InDel'}"     :  "";
    my $knownSitesInDelGold = (exists $hashPara->{$Species}{"InDelGold"}) ? "--known-sites $hashPara->{$Species}{'InDelGold'}" :  "";
    
    my $fastq_R1 = "$fastq_dir/$sample\_R1.fastq.gz";
    my $fastq_R2 = "$fastq_dir/$sample\_R2.fastq.gz"; 

    # 比对   
    my $bam = "$output_dir/$sample.bam";
    system "$BWA mem -t 4 -k 32 -M -T 50 -B 8 -O 12 -L 10  -U 18 -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:illumina\\tPU:barcode' $bwa_idx $fastq_R1 $fastq_R2 | $SamTools view -bS -h - > $bam";             
    
    # 排序
    my $bam_sort     = "$output_dir/$sample\_sort.bam";
    system "$sambamba sort $bam -m 6G -t 4 -o $bam_sort --tmpdir $tmp_dir";
    # system "cpulimit -l 1000 $Java -jar -Xmx6g $Picard SortSam INPUT=$bam OUTPUT=$bam_sort CREATE_INDEX=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir";
    
    # 质量统计
    my $cover_metrics   = "$status_dir/$sample.metrics"; 
    my $cover_tag       = "$status_dir/$sample.tag";
    my $quality_metrics = "$status_dir/$sample.quality.metrics";
    system "$cpulimit -l 1000 $Java -jar -Xmx3g $Picard CollectHsMetrics            INPUT=$bam_sort OUTPUT=$cover_metrics R=$Genome BI=$TargetBed TI=$TargetBed PER_TARGET_COVERAGE=$cover_tag MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir";
    system "$cpulimit -l 1000 $Java -jar -Xmx3g $Picard CollectQualityYieldMetrics  INPUT=$bam_sort OUTPUT=$quality_metrics VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir";
    
    # 文库序列长度统计
    my $insert     = "$status_dir/$sample.insert";
    my $insert_pdf = "$status_dir/$sample.insert.pdf";
    system "$cpulimit -l 1000 $Java -jar -Xmx3g $Picard CollectInsertSizeMetrics INPUT=$bam_sort O=$insert H=$insert_pdf VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir";
    
    # 去重复
    my $bam_sort_dup         = "$output_dir/$sample\_sort_dup.bam";
    my $bam_sort_dup_metrics = "$status_dir/$sample\_sort_dup.metrics";
    system "$cpulimit -l 1000 $Java -jar -Xmx16g $Picard MarkDuplicates INPUT=$bam_sort OUTPUT=$bam_sort_dup METRICS_FILE=$bam_sort_dup_metrics CREATE_INDEX=true REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir";
   
    # 根据数据库，对已知插入缺失位点校正，同时，只保留目标区域的数据，减少数据量(如果后面分析突变时，只使用GATK的HaplotypeCaller方法，这一步就没有必要做)(实际上，非目标区域的reads也保留了一部分，原因未知)
    # 使用HaplotypeCaller 分析SNV，可以不用做该步骤
    my $bam_sort_dup_realign           ="$output_dir/$sample\_sort_dup_realign.bam";
    my $bam_sort_dup_realign_intervals ="$output_dir/$sample\_sort_dup_realign.intervals";
    system "$cpulimit -l 1000 $Java -Djava.io.tmpdir=$tmp_dir -jar -Xmx5g $GATK3 -l INFO -T RealignerTargetCreator -R $Genome -I $bam_sort_dup -L $RealignBed -o $bam_sort_dup_realign_intervals $knownInDel $knownInDelGold --validation_strictness LENIENT";
    system "$cpulimit -l 1000 $Java -Djava.io.tmpdir=$tmp_dir -jar -Xmx5g $GATK3 -l INFO -T IndelRealigner         -R $Genome -I $bam_sort_dup -L $RealignBed -o $bam_sort_dup_realign           $knownInDel $knownInDelGold  -targetIntervals $bam_sort_dup_realign_intervals --validation_strictness LENIENT";
    
    # 碱基质量校正 （与GATK3 相比，时间从7小时缩短至1小时）
    my $bam_sort_dup_realign_recalibrator_table = "$output_dir/$sample\_sort_dup_realign_recalibrator.table";
    my $bam_sort_dup_realign_recalibrator       = "$output_dir/$sample\_sort_dup_realign_recalibrator.bam";
    my $bam_sort_dup_realign_recalibrator_idx   = "$output_dir/$sample\_sort_dup_realign_recalibrator.bai";
    system("$cpulimit -l 1000 $Java -jar -Xmx5g $GATK4_Loc BaseRecalibrator -R $Genome -I $bam_sort_dup_realign  -O $bam_sort_dup_realign_recalibrator_table $knownSitesDBsnp $knownSitesInDelGold $knownSitesInDel --tmp-dir $tmp_dir") if($knownSitesDBsnp=~/\w/); # 缺失VCF时，无法进行该步骤
    system("$cpulimit -l 1000 $Java -jar -Xmx8g $GATK4_Loc ApplyBQSR        -R $Genome -I $bam_sort_dup_realign  -O $bam_sort_dup_realign_recalibrator --bqsr-recal-file $bam_sort_dup_realign_recalibrator_table  --tmp-dir $tmp_dir")  if($knownSitesDBsnp=~/\w/);
    
    # 最终比对bam
    my $bam_final_input     = ($knownSitesDBsnp=~/\w/) ? $bam_sort_dup_realign_recalibrator : $bam_sort_dup;
    my $bam_final_input_idx = $bam_final_input;
       $bam_final_input_idx =~s/\.bam$/\.bai/;
    # 重命名
    my $bam_final     = "$output_dir/$sample\_final.bam";
    my $bam_final_idx = "$output_dir/$sample\_final.bai";
    system("mv $bam_final_input     $bam_final");
    system("mv $bam_final_input_idx $bam_final_idx");

}

1