package package::target_typing;
use strict;
use warnings;
$|=1;


sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start TargetTyping ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::target_typing($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in TargetTyping run para, please check\n" if($isOK == 0);


    my $report_dir       = $hashConfig->{'Report'}; 
    my $output_dir       = $hashConfig->{'Output'}; 

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control', 'Ycontrol');
    my %hashCondition; 
    foreach my $sample(@samples)
    {

        my $blast       = "$output_dir/$sample/$sample.blast";# 比对结果
        my $fastq       = "$output_dir/$sample/$sample.fastq";# 比对fastq
        my $readsBelong = "$output_dir/$sample/$sample.reads.belong";# reads所属片段
        my $baseFilter  = "$blast.base.filter";# 最终分析结果
        # 已完成该步骤
        if(package::utils::is_file_ok($baseFilter) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$baseFilter";
            next;            
        }
        # 原始数据没问题
        if(package::utils::is_file_ok($blast) and package::utils::is_file_ok($fastq) and package::utils::is_file_ok($readsBelong)){
            $hashCondition{"Good2Run"}{$sample} = "$blast $fastq $readsBelong";
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$blast $fastq $readsBelong";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "TargetTyping", \%hashCondition); 
    die "Sample lost blast/fastq/reads.belong :" . (join ",", sort keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));

    # 执行比对
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0)
    {   
        my $threshold = exists $hashPara->{"Process"}{"TargetTyping"} ? $hashPara->{"Process"}{"TargetTyping"} : 10;
           $threshold = $hashConfig->{"Process_TargetTyping"} if(exists $hashConfig->{"Process_TargetTyping"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::utils::process_bar_array($sample, \@sample_runs)});# 进度条
        foreach my $sample(@sample_runs)
        {
            $pm->start($sample) and next;
            run_target_typing($hashPara, $hashConfig, $sample);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else
    {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}


sub run_target_typing{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $sample     = shift @_;

    my $report_dir   = $hashConfig->{'Report'};
    my $output_dir   = "$hashConfig->{'Output'}/$sample";
    my $seq_region   = "$report_dir/seq.region";
    my $blast             = "$output_dir/$sample.blast";# 比对结果
    my $fastq             = "$output_dir/$sample.fastq";# 比对fastq
    my $reads_belong      = "$output_dir/$sample.reads.belong";# reads所属片段
    my $reads_base        = "$blast.base.gz"; # 识别比对位置上的每一个碱基（由“位置，碱基，质量”构成，其中位置是参考序列上的位置）
    my $reads_base_filter = "$blast.base.filter.gz";# R1,R2重叠碱基合并，取高质量结果
    
    #####
    # 1.识别比对位置上的每一个碱基
    #####
    my %hashReadsBelong = read_reads_belong($reads_belong);# 经过上一步的质控，提取reads属于哪一个参考片段
    my %hashFastq       = read_fastq_quality($fastq); # fastq数据ASCII质量序列
    my %hashSeqRegion   = package::utils::read_seq_region($seq_region);
    my %hashFinishReads = (); # 保存分析过的reads名，用于过滤
    
    print "Output $reads_base start\n";
    open READSBASE, "|gzip >$reads_base";
    my $in = Bio::SearchIO->new(-file => $blast, -format => 'blast');
    while( my $r = $in->next_result ) 
    {
        while( my $h = $r->next_hit ) 
        {  
            while( my $hsp = $h->next_hsp ) 
            {  
                my $target_name   =  $h->name;
                   $target_name   =~ s/\s+//g; # 片段名称
                my $reads_name    =  $r->query_name; # reads名称
                last if(exists($hashFinishReads{$reads_name})); # 该reads已经分析完毕
                last if($target_name ne $hashReadsBelong{$reads_name}); # 该reads质控时对应的参考片段不是当前的参考片段,解决片段重叠问题
                
                my $target_string               = lc($hsp->hit_string);   # 比对结果，参考序列
                my $reads_string                = lc($hsp->query_string); # 比对结果，fastq序列
                my ($target_start, $target_end) = $hsp->range('hit');   # 始终是从小到大
                my ($reads_start, $reads_end)   = $hsp->range('query'); # 始终是从小到大
                my $hit_length                  = length($target_string); # 比对序列长度
                my $strand                      = $hsp->strand('hit'); # 比对方向
                if($strand == -1) # 方向问题,需要保证片段始终是与参照片段方向一致，即从小到大
                {    
                   $reads_string  = reverse $reads_string;
                   $reads_string  =~ tr/ATCGatcg/TAGCtagc/;
                   $target_string =  reverse $target_string;    
                   $target_string =~ tr/ATCGatcg/TAGCtagc/;                                            
                }
                # 识别
                my @identify_bases; # 识别的碱基
                my $target_position = $target_start; # 参照序列初始位置
                my $reads_position  = ($strand == -1) ? $reads_end : $reads_start; # reads序列初始位置
                for(my $i=1; $i <= $hit_length; $i++) # 遍历参考序列
                {
                    my $target_base = substr($target_string, $i - 1, 1); # 当前位置参照碱基
                    my $reads_base  = substr($reads_string, $i - 1, 1);  # 当前位置reads碱基 
                    my $ref_base    = substr($hashSeqRegion{$target_name}{'TargetSeq'}, $target_position - 1, 1); # 当前位置参照碱基（C未转化未T）
                    next if($i == $hit_length and $ref_base =~ /c/i); # 当前位置是该比对reads的最后一个碱基，且该位置在参考序列上是C。这么处理的原因在于，我们的参考序列所有的C->T，而我们是用blastn进行比对的，这就导致如果reads尾部是C的情况下，该碱基比对不上，如果是T，则能够比对上，产生T的偏好性，使得结果有问题
                    if($target_base =~ /[atcg]/i and $reads_base =~ /[actg]/i) # 参照与reads都是碱基
                    {
                        my $base_quality = substr($hashFastq{$reads_name}, $reads_position - 1, 1); # 当前碱基质量
                        push @identify_bases, "$target_position,$reads_base,$base_quality";
                    }
                    $target_position++ if($target_base =~ /[atcgn]/i);
                    $reads_position--  if($strand == -1 and $reads_base =~ /[atcgn]/i);
                    $reads_position++  if($strand == 1  and $reads_base =~ /[atcgn]/i);
                }
                $hashFinishReads{$reads_name}++;# reads已读过
                print READSBASE "$reads_name\t$target_name\t" . (join ";", @identify_bases) . "\n";
            }
        }
    }
    close READSBASE;
    print "Output $reads_base finish\n";
    undef %hashFastq; # 释放内存，内存会被标记为reusable。内存不会被立即释放，只有当系统内存不足时，这部分内存才会被释放。
    undef %hashFinishReads;
    undef %hashReadsBelong;

    #####
    # 2.R1,R2重叠碱基合并，取高质量结果
    #####
    my %hashReadsBase;# 每一条reads的碱基信息
    print "Read $reads_base start\n";
    open READSBASE, "gzip -cd $reads_base|";
    my $row = 0;
    while(<READSBASE>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        $row++;
        my ($reads_name, $target_name, $identify_base_string) = split /\t/, $_;
        my @tmps = split /:/, $reads_name;
        $reads_name = join ":", @tmps[0..$#tmps-1]; # reads名称
        my $merge_type = $tmps[$#tmps];# merge类型，存在3中值：MERGE,R1,R2
        $hashReadsBase{$reads_name}{$merge_type}   = $identify_base_string;
        $hashReadsBase{$reads_name}{'ReadsBelong'} = $target_name;
        $hashReadsBase{$reads_name}{'SortOrder'}   = $row if(not exists $hashReadsBase{$reads_name}{'SortOrder'} );
    }
    close READSBASE;
    print "Read $reads_base finish\n";

    # 重叠合并，输出
    print "Output $reads_base_filter start\n";
    open READSBASEFILTER, "|gzip >$reads_base_filter";
    foreach my $reads_name(sort {$hashReadsBase{$a}{'SortOrder'} <=> $hashReadsBase{$b}{'SortOrder'}} keys %hashReadsBase)
    {   
        # MERGE 不需要再处理
        if(exists $hashReadsBase{$reads_name}{'MERGE'})
        {
            print READSBASEFILTER "$reads_name:MERGE\t$hashReadsBase{$reads_name}{'ReadsBelong'}\t$hashReadsBase{$reads_name}{'MERGE'}\n";
            next;
        }

        # 处理R1/R2
        my $identifyBaseInfo = '';# 碱基合并
        my $newMergeType     = '';# reads合并类型
        my $readsBelong      = '';# reads所属片段
        foreach my $mergeType(('R1','R2'))# 顺序不能变
        {
            next if(not exists $hashReadsBase{$reads_name}{$mergeType} );
            $identifyBaseInfo = "$identifyBaseInfo;$hashReadsBase{$reads_name}{$mergeType}"; # 存在先后顺序。后面的分析中，同质量的情况下，优先保留已记录的结果。
            $newMergeType     = "$newMergeType$mergeType";
            $readsBelong      = $hashReadsBase{$reads_name}{'ReadsBelong'} if($readsBelong eq '');
        }
        # 合并
        my %hashBase;
        foreach my $identifyBase(split /;/, $identifyBaseInfo)
        {
            next if($identifyBase !~ /\w/);
            my ($position, $base, $quality) = split /,/, $identifyBase;
            $quality = ord($quality);
            if( (not exists $hashBase{$position}) or  $quality > $hashBase{$position}{'Quality'} )# 没有保存过 或者保存过，但新的质量更高一些
            {
                $hashBase{$position}{'Base'}    = $base;
                $hashBase{$position}{'Quality'} = $quality;             
            }
        }

        my @identifyBases = map{ my $quality = chr($hashBase{$_}{'Quality'}); "$_,$hashBase{$_}{'Base'},$quality" } sort {$a <=> $b} keys %hashBase;
        print READSBASEFILTER "$reads_name:$newMergeType\t$readsBelong\t" . (join ";", @identifyBases) . "\n";
    }
    close READSBASEFILTER;
    print "Output $reads_base_filter finish\n";
}
# 读fastq质量
sub read_fastq_quality{
    my $fastq = shift @_;
    print "Read $fastq start\n";
    my %hashFastq;
    my $in     = Bio::SeqIO->new(-file => $fastq, -format=>'Fastq') or die "Could not open up file $fastq: $!";
    while(my $inSeq = $in->next_seq)
    {
        my $id   = $inSeq->display_id;
        my $seq  = $inSeq->seq;
        my $qual = join('', map{chr(33+$_)} @{$inSeq->qual} );# 质量值
        $hashFastq{$id} = $qual;
    }
    print "Read $fastq finish\n";
    return %hashFastq;
}
# 读fasta序列
sub readFasta{
    my $fasta = shift @_;
    my %hashFasta;
    my $in     = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    while(my $inSeq = $in->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        $hashFasta{$id}=$seq;
    }
    return %hashFasta;  
}
# 读reads所属
sub read_reads_belong{
    my $reads_belong = shift @_;
    my %hashReadsBelong;
    open BELONG, $reads_belong;
    while(<BELONG>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/);
        my ($reads_name, $target_name) = split /\t/, $_;
        $hashReadsBelong{$reads_name}  = $target_name;
    }
    close BELONG;
    return %hashReadsBelong;
}

1