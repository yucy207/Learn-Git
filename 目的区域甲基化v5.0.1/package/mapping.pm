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

    my $report_dir       = $hashConfig->{'Report'}; 
    my $fastq_dir        = $hashConfig->{'Fastq'}; 
    my $output_dir       = $hashConfig->{'Output'}; 
    my $status_dir       = "$report_dir/status";
    my $fastqc_clean_dir = "$report_dir/fastqc_clean";
    package::utils::make_dir($status_dir);
    package::utils::make_dir($fastqc_clean_dir);

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control', 'Ycontrol');
    my %hashCondition; 
    foreach my $sample(@samples)
    {
        my $blast_out = "$output_dir/$sample/$sample.blast";
        my $fastq_R1  = "$fastq_dir/$sample\_R1.fastq.gz";
        my $fastq_R2  = "$fastq_dir/$sample\_R2.fastq.gz";
        # 已完成该步骤
        if(package::utils::is_file_ok($blast_out) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$blast_out";
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

    # 目录
    my $report_dir       = $hashConfig->{'Report'}; 
    my $fastq_dir        = $hashConfig->{'Fastq'}; 
    my $output_dir       = "$hashConfig->{'Output'}/$sample"; 
    my $status_dir       = "$report_dir/status";
    my $fastqc_clean_dir = "$report_dir/fastqc_clean";
    package::utils::make_dir($output_dir);

    # 软件
    my $FastX      = $hashPara->{"Soft"}{"FastX"};
    my $BlastPlus  = "$hashPara->{'Soft'}{'BlastPlus'}/blastn";
    my $FLASH      = $hashPara->{"Soft"}{"FLASH"};
    my $FastQC     = $hashPara->{"Soft"}{"FastQC"};
    my $FastqStat  = $hashPara->{"Soft"}{"FastqStat"};


    # 基础输入文件
    my $target_fasta_db = "$report_dir/methylation_fasta_db/methylation.fa";
    my $fastq_R1        = "$fastq_dir/$sample\_R1.fastq.gz";
    my $fastq_R2        = "$fastq_dir/$sample\_R2.fastq.gz"; 

    # Merge
    print "1.Merge Fastq : $sample\n";
    my $FLASH_Output_dir = "$output_dir/flash";
    package::utils::make_dir($FLASH_Output_dir);
    system "$FLASH -t 4 --allow-outies -m 15 -x 0.1 -d $FLASH_Output_dir -o $sample $fastq_R1 $fastq_R2";

    # 数据过滤
    print "2.Filter Fastq : $sample \n";
    my $merge_fastq      = "$FLASH_Output_dir/$sample.extendedFrags.fastq";
    my $R1_nomerge_fastq = "$FLASH_Output_dir/$sample.notCombined_1.fastq";
    my $R2_nomerge_fastq = "$FLASH_Output_dir/$sample.notCombined_2.fastq";
    my %hashTargetFasta  = package::utils::read_fasta($target_fasta_db);
    my %hashNoMergeResult = nomerge_fastq_analysis($BlastPlus, $FastX, $target_fasta_db, \%hashTargetFasta, $FLASH_Output_dir, $sample, $hashConfig); # 合并失败reads过滤
    my %hashMergeResult   = merge_fastq_analysis($BlastPlus, $FastX, $target_fasta_db, \%hashTargetFasta, $FLASH_Output_dir, $sample); # 合并成功reads过滤
    
    # 过滤结果输出 
    print "3.Write Filter Fastq : $sample \n";
    write_nomerge_data(\%hashNoMergeResult, $FLASH_Output_dir, $sample, "UNMAPPED");
    write_nomerge_data(\%hashNoMergeResult, $FLASH_Output_dir, $sample, "MAPPED"); 
    write_merge_data(\%hashMergeResult, $FLASH_Output_dir, $sample, "UNMAPPED");
    write_merge_data(\%hashMergeResult, $FLASH_Output_dir, $sample, "MAPPED");

    # 统计结果输出
    my $Status       = "$status_dir/$sample.status";
    my $reads_belong = "$output_dir/$sample.reads.belong";
    write_status(\%hashTargetFasta, \%hashMergeResult, \%hashNoMergeResult, $Status);
    write_reads_belong(\%hashMergeResult, \%hashNoMergeResult, $reads_belong); # 合格reads比对到的target信息输出

    # 质控过滤后，合格数据合并
    print "4.Output final Fastq : $sample \n";
    my $merge_MAPPED_fastq      = "$FLASH_Output_dir/$sample.merge_MAPPED.fastq";
    my $nomerge_r1_MAPPED_fastq = "$FLASH_Output_dir/$sample.r1_MAPPED.fastq";
    my $nomerge_r2_MAPPED_fastq = "$FLASH_Output_dir/$sample.r2_MAPPED.fastq";
    my $final_fastq             = "$output_dir/$sample.fastq"; # 最终fastq
    system("cat $merge_MAPPED_fastq $nomerge_r1_MAPPED_fastq $nomerge_r2_MAPPED_fastq > $final_fastq");

    # 最终fastq比对
    print "5.blast final Fastq : $sample \n";
    my $final_fasta = "$output_dir/$sample.fasta";
    my $final_blast = "$output_dir/$sample.blast";
    system "$FastX/fastq_to_fasta -Q 33 -n -i $final_fastq -o $final_fasta";
    system("$BlastPlus -task blastn -query $final_fasta -evalue 0.00001 -db $target_fasta_db -out $final_blast -num_threads 4 -dust no");

    # 合格的reads再质控
    print "6.fastqc final Fastq : $sample \n";
    system "$FastQC --extract -q -o $fastqc_clean_dir $final_fastq";
    system "rm $fastqc_clean_dir/$sample\_fastqc.zip";
    package::qc::QC2030($final_fastq, "$fastqc_clean_dir/$sample\_fastqc/fastqc_stat.txt", $FastqStat);
}

sub write_status{
    my $hashTargetFasta   = shift @_;
    my $hashMergeResult   = shift @_;
    my $hashNoMergeResult = shift @_;
    my $Status            = shift @_;
    open SAVE,">$Status";
    my $merge_bad    = (exists $hashMergeResult->{"MergeBad"})     ? $hashMergeResult->{"MergeBad"}     : 0;
    my $no_merge_bad = (exists $hashNoMergeResult->{"NoMergeBad"}) ? $hashNoMergeResult->{"NoMergeBad"} : 0;

    my $bad  = $merge_bad + $no_merge_bad;
    my $good = 0;
    foreach my $target_name(keys %{$hashTargetFasta})
    {
        my $merge_good    = (exists $hashMergeResult->{"Merge"}{$target_name})       ? $hashMergeResult->{"Merge"}{$target_name}     : 0;
        my $no_merge_good = (exists $hashNoMergeResult->{"NoMerge"}{$target_name})   ? $hashNoMergeResult->{"NoMerge"}{$target_name} : 0;
        my $good_tmp      = $merge_good + $no_merge_good;

        print SAVE "$target_name\t$good_tmp\n";
        $good += $good_tmp;
    }
    print SAVE "###$good\t".($good + $bad)."\n";
    close SAVE;
}
sub write_reads_belong{
    my $hashMergeResult   = shift @_;
    my $hashNoMergeResult = shift @_;
    my $reads_belong      = shift @_;
    open BELONG,">$reads_belong";
    foreach my $reads_name(sort keys %{$hashMergeResult->{'MergeBelong'}})
    {
        print BELONG "$reads_name\t$hashMergeResult->{'MergeBelong'}{$reads_name}\n";
    }
    foreach my $reads_name(sort keys %{$hashNoMergeResult->{'NoMergeBelong'}})
    {
        print BELONG "$reads_name\t$hashNoMergeResult->{'NoMergeBelong'}{$reads_name}\n";
    }
    close BELONG;
}
# 输出merge失败的reads过滤后的结果
sub write_nomerge_data{
    my $hashNoMergeResult = shift @_;
    my $FLASH_Output_dir  = shift @_;
    my $sample            = shift @_;
    my @out_types         = @_; # "MAPPED", "UNMAPPED"
    foreach my $out_type(@out_types)
    {
        open SAVE_R1, ">$FLASH_Output_dir/$sample.r1_$out_type.fastq";
        open SAVE_R2, ">$FLASH_Output_dir/$sample.r2_$out_type.fastq";
        foreach my $reads_name(keys %{$hashNoMergeResult->{$out_type}})
        {
            for my $i(1..4)
            {   
                # 注意：有可能会强制处理R1/R2,此时MAPPED/UNMAPPED会出现只有单端的情况
                print SAVE_R1 $hashNoMergeResult->{$out_type}{$reads_name}{"R1"}{$i} . "\n" if(exists $hashNoMergeResult->{$out_type}{$reads_name}{"R1"});
                print SAVE_R2 $hashNoMergeResult->{$out_type}{$reads_name}{"R2"}{$i} . "\n" if(exists $hashNoMergeResult->{$out_type}{$reads_name}{"R2"});
            }
        }
        close SAVE_R1;
        close SAVE_R2;
    }
}
# 输出merge成功的reads过滤后的结果
sub write_merge_data{
    my $hashMergeResult  = shift @_;
    my $FLASH_Output_dir = shift @_;
    my $sample           = shift @_;
    my @out_types        = @_; # "MAPPED", "UNMAPPED"
    foreach my $out_type(@out_types)
    {
        open SAVE, ">$FLASH_Output_dir/$sample.merge_$out_type.fastq";
        foreach my $reads_name(keys %{$hashMergeResult->{$out_type}})
        {
            for my $i(1..4)
            {
                print SAVE $hashMergeResult->{$out_type}{$reads_name}{$i} . "\n";
            }
        }
        close SAVE;
    }
}
# 合并成功reads过滤
sub merge_fastq_analysis{
    my $BlastPlus        = shift @_; 
    my $FastX            = shift @_;
    my $target_fasta_db  = shift @_;
    my $hashTargetFasta  = shift @_;
    my $FLASH_Output_dir = shift @_;
    my $sample           = shift @_;

    # fastq 转 fasta
    my $merge_fastq = "$FLASH_Output_dir/$sample.extendedFrags.fastq";
    my $merge_fasta = "$FLASH_Output_dir/$sample.extendedFrags.fa";   
    system "$FastX/fastq_to_fasta -Q 33 -n -i $merge_fastq -o $merge_fasta";

    # blast
    my $merge_blast = "$FLASH_Output_dir/$sample.extendedFrags.blast";
    system("$BlastPlus -task blastn -query $merge_fasta -outfmt 6 -evalue 0.00001 -db $target_fasta_db -out $merge_blast -num_threads 4 -max_target_seqs 5 -dust no");# max_target_seqs 只给出最先比配到的5个数据库染色体序列的结果，一个染色体可能有多个结果。注意：具有随机性，可能导致无法得到最佳匹配结果
    
    # 筛选
    my %hashmergeFastq = read_fastq($merge_fastq);
    my %hashmergeBlast = read_blast(\%hashmergeFastq, $hashTargetFasta, $merge_blast, 'Merge');
    my %hashMergeResult=();
    foreach my $reads_name(keys %hashmergeFastq)
    {
        if(exists $hashmergeBlast{$reads_name})
        {
            my @split_info = split /\t/, $hashmergeBlast{$reads_name};

            $hashmergeFastq{$reads_name}{1} = "\@$reads_name:MERGE";

            $hashMergeResult{"MAPPED"}{$reads_name} = $hashmergeFastq{$reads_name};
            $hashMergeResult{"Merge"}{$split_info[1]}++;
            $hashMergeResult{"MergeBelong"}{"$reads_name:MERGE"} = $split_info[1];# 记录当前reads比对到了哪一个片段
        }
        else
        {
            $hashMergeResult{"UNMAPPED"}{$reads_name} = $hashmergeFastq{$reads_name};
            $hashMergeResult{"MergeBad"}++;
        }
    }
    return %hashMergeResult;
}


# 合并失败reads过滤
sub nomerge_fastq_analysis{
    my $BlastPlus        = shift @_; 
    my $FastX            = shift @_;
    my $target_fasta_db  = shift @_;
    my $hashTargetFasta  = shift @_;
    my $FLASH_Output_dir = shift @_;
    my $sample           = shift @_;
    my $hashConfig       = shift @_;

    # fastq 转 fasta
    my $R1_nomerge_fastq = "$FLASH_Output_dir/$sample.notCombined_1.fastq";
    my $R2_nomerge_fastq = "$FLASH_Output_dir/$sample.notCombined_2.fastq";
    my $R1_nomerge_fasta = "$FLASH_Output_dir/$sample.notCombined_1.fa";
    my $R2_nomerge_fasta = "$FLASH_Output_dir/$sample.notCombined_2.fa";    
    system "$FastX/fastq_to_fasta -Q 33 -n -i $R1_nomerge_fastq -o $R1_nomerge_fasta";
    system "$FastX/fastq_to_fasta -Q 33 -n -i $R2_nomerge_fastq -o $R2_nomerge_fasta";

    # blast
    my $R1_nomerge_blast = "$FLASH_Output_dir/$sample.notCombined_1.blast";
    my $R2_nomerge_blast = "$FLASH_Output_dir/$sample.notCombined_2.blast";
    system("$BlastPlus -task blastn -query $R1_nomerge_fasta -outfmt 6 -evalue 0.00001 -db $target_fasta_db -out $R1_nomerge_blast -num_threads 4 -max_target_seqs 5 -dust no");
    system("$BlastPlus -task blastn -query $R2_nomerge_fasta -outfmt 6 -evalue 0.00001 -db $target_fasta_db -out $R2_nomerge_blast -num_threads 4 -max_target_seqs 5 -dust no");
            
    # 筛选
    my %hashR1NomergeFastq = read_fastq($R1_nomerge_fastq);
    my %hashR2NomergeFastq = read_fastq($R2_nomerge_fastq);
    my %hashR1NomergeBlast = read_blast(\%hashR1NomergeFastq, $hashTargetFasta, $R1_nomerge_blast, 'NoMerge');
    my %hashR2NomergeBlast = read_blast(\%hashR2NomergeFastq, $hashTargetFasta, $R2_nomerge_blast, 'NoMerge');
    my %hashNoMergeResult=();
    foreach my $reads_name(keys %hashR1NomergeFastq)
    {
        my @split_info1;
        my @split_info2;
           @split_info1 = split /\t/, $hashR1NomergeBlast{$reads_name} if(exists $hashR1NomergeBlast{$reads_name});
           @split_info2 = split /\t/, $hashR2NomergeBlast{$reads_name} if(exists $hashR2NomergeBlast{$reads_name});
        # 1 R1/R2都有最佳比对结果
        if(exists $hashR1NomergeBlast{$reads_name} and exists $hashR2NomergeBlast{$reads_name})
        {
            # 1.1 R1/R2最佳比对结果相同
            if($split_info1[1] eq $split_info2[1])
            {
                # 更换序列名称
                $hashR1NomergeFastq{$reads_name}{1} = "\@$reads_name:R1";
                $hashR2NomergeFastq{$reads_name}{1} = "\@$reads_name:R2";
                
                # 保留数据
                $hashNoMergeResult{"MAPPED"}{$reads_name}{"R1"} = $hashR1NomergeFastq{$reads_name};
                $hashNoMergeResult{"MAPPED"}{$reads_name}{"R2"} = $hashR2NomergeFastq{$reads_name};
                $hashNoMergeResult{"NoMerge"}{$split_info1[1]}++;
                $hashNoMergeResult{"NoMergeBelong"}{"$reads_name:R1"} = $split_info1[1];# 记录当前reads比对到了哪一个片段
                $hashNoMergeResult{"NoMergeBelong"}{"$reads_name:R2"} = $split_info1[1];# 记录当前reads比对到了哪一个片段
            }
            else
            {

                $hashNoMergeResult{"UNMAPPED"}{$reads_name}{"R1"}=$hashR1NomergeFastq{$reads_name};
                $hashNoMergeResult{"UNMAPPED"}{$reads_name}{"R2"}=$hashR2NomergeFastq{$reads_name};
                $hashNoMergeResult{"NoMergeBad"}++;
            }
        }
        # 2 R1/R2 存在没有比对上的
        else
        {
            ############# 
            ## 强制处理部分片段R1,R2（当测序质量较差，只有R1或R2端能够完美匹配）
            #############
            # R1
            my $force = 0;
            if(exists($split_info1[1]) and exists($hashConfig->{'Force'}{"$split_info1[1]:R1"}))
            {
                $hashR1NomergeFastq{$reads_name}{1}             = "\@$reads_name:R1";
                $hashNoMergeResult{"MAPPED"}{$reads_name}{"R1"} = $hashR1NomergeFastq{$reads_name};
                $hashNoMergeResult{"NoMerge"}{$split_info1[1]}++;
                $hashNoMergeResult{"NoMergeBelong"}{"$reads_name:R1"} = $split_info1[1];# 记录当前reads比对到了哪一个片段
                $hashNoMergeResult{"UNMAPPED"}{$reads_name}{"R2"}     = $hashR2NomergeFastq{$reads_name};
                $force++;
            }
            # R2
            if(exists($split_info2[1]) and exists($hashConfig->{'Force'}{"$split_info2[1]:R2"}))
            {
                $hashR2NomergeFastq{$reads_name}{1}             = "\@$reads_name:R2";
                $hashNoMergeResult{"MAPPED"}{$reads_name}{"R2"} = $hashR2NomergeFastq{$reads_name};
                $hashNoMergeResult{"NoMerge"}{$split_info2[1]}++;
                $hashNoMergeResult{"NoMergeBelong"}{"$reads_name:R2"} = $split_info2[1];# 记录当前reads比对到了哪一个片段
                $hashNoMergeResult{"UNMAPPED"}{$reads_name}{"R1"}     = $hashR1NomergeFastq{$reads_name};
                $force++;                
            }
            next if($force > 0);

            $hashNoMergeResult{"UNMAPPED"}{$reads_name}{"R1"} = $hashR1NomergeFastq{$reads_name};
            $hashNoMergeResult{"UNMAPPED"}{$reads_name}{"R2"} = $hashR2NomergeFastq{$reads_name};
            $hashNoMergeResult{"NoMergeBad"}++;
        }
    }
    return %hashNoMergeResult;
}

# 读取fastq文件
sub read_fastq{
    my $fastq = shift @_;
    my %hashFastq;
    open FASTQ,$fastq;
    while(my $line1 = <FASTQ>)
    {
        my $line2 = <FASTQ>;
        my $line3 = <FASTQ>;
        my $line4 = <FASTQ>;
        $line1=~ s/[\r\n]//g;
        $line2=~ s/[\r\n]//g;
        $line3=~ s/[\r\n]//g;
        $line4=~ s/[\r\n]//g;
        my $name=(split /\s/,$line1)[0];
           $name=~ s/^\@//;
        $hashFastq{$name}{1} = $line1;
        $hashFastq{$name}{2} = $line2;
        $hashFastq{$name}{3} = $line3;
        $hashFastq{$name}{4} = $line4;
    }
    close FASTQ;
    return %hashFastq;
}

# 读取blast结果
sub read_blast{
    my $hashFastq       = shift @_; 
    my $hashTargetFasta = shift @_;
    my $blast           = shift @_;
    my $type            = shift @_; # Merge表示对merge的reads进行比对。NoMerge表示对没有merge上的R1或R2的reads进行比对。
    
    # 读取blast比对信息，并进行基础过滤
    my %hashBlast;
    open BLAST,$blast;
    while(<BLAST>){
        $_=~ s/[\r\n]//g;
        my ($reads_name, $target_name, $identify_perc, $identify_len, $mismatch, $gap, $reads_start, $reads_end, $target_start, $target_end, $evalue, $score) = split /\t/, $_;
        my @array=split(/\t/,$_);
        my $reads_length  = length( $hashFastq->{$reads_name}{2} );
        my $target_length = length( $hashTargetFasta->{$target_name} );
        ($target_start, $target_end) = ($target_end, $target_start) if($target_end < $target_start);
 
        my $judge_value = $identify_len; # 默认用匹配长度作为后续候选的筛选规则
        #
        # 情况一：Reads覆盖整个目标区域
        # 条件，覆盖超过90%的区域
        # 条件，覆盖起始位置小于开头加10，结束位置大于末尾减10
        #
        if($type eq 'Merge' and $identify_len > $target_length * 0.9 and $target_start <= 10 and $target_end > $target_length - 10){
            $hashBlast{$reads_name}{$identify_len} = "$reads_name\t$target_name\t$reads_start\t$reads_end\t$target_start\t$target_end\t$identify_len\n";
        }
 
        # 情况二：Reads在目标区域的一端
        # 条件，Reads长度覆盖超过90%
        # 条件，覆盖起始位置小于开头加10，或结束位置大于末尾减10
        #
        if($type eq 'NoMerge' and $identify_len > $reads_length * 0.9 and ($target_start <= 10 or $target_end > $target_length - 10) ){
            $hashBlast{$reads_name}{$identify_len} = "$reads_name\t$target_name\t$reads_start\t$reads_end\t$target_start\t$target_end\t$identify_len\n";
        }
    }
    close BLAST;

    # 保留比对最长的
    my %hashBest=();
    foreach my $title(keys %hashBlast){
        foreach my $len (sort {$b <=> $a} keys %{$hashBlast{$title}}){
            $hashBest{$title} = $hashBlast{$title}{$len};
            last;
        }
    }       
    return %hashBest;
}



1