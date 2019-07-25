package package::check_methylation;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start check_methylation ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::check_methylation($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in check_methylation run para, please check\n" if($isOK == 0);

    my $report_dir = $hashConfig->{'Report'}; 

    # 1 建立甲基化序列索引
    build_methylation_db($hashConfig, $hashPara);      # 建立甲基化参考序列库

    # 2 检测fasta序列在参考基因组上的位置信息
    check_fasta_position($hashConfig, $hashPara);      # 检测参考序列基因组信息

    # 3 检查fasta序列 CG/CHG/CHH位点信息
    check_fasta_methyl_position($hashConfig, $hashPara);

    
}

sub check_fasta_methyl_position{
    my $hashConfig    = shift @_;
    my $hashPara      = shift @_;

    print "\n3: Check Fasta methylation Position and generate [c.point.analysis.txt] / [c.point.other.txt] file\n";


    my $species      = package::utils::get_species($hashConfig);
    my $special_mark = exists $hashPara->{$species}{'SpecialMark'}  ? $hashPara->{$species}{'SpecialMark'}  : "";

    my $report_dir     = $hashConfig->{'Report'};
    my $seq_region     = "$report_dir/seq.region";
    my %hashSeqRegion  = package::utils::read_seq_region($seq_region); 

    # 序列默认前后去掉的长度
    my $cut_F_default  = (exists $hashConfig->{'TargetCut'} ) ? $hashConfig->{'TargetCut'} : 25;
    my $cut_R_default  = $cut_F_default;

    # 项目分析的C类型,默认CG
    my $analysis_type = ( exists $hashConfig->{'AnalysisType'} ) ? $hashConfig->{'AnalysisType'} : 'CG';

    # (1)识别序列中的C
    print "    check methylation point\n";
    my %hashAnalysisC; # 要分析的C
    my %hashOtherC;    # 其他的C,主要用于检测盐化效率
    foreach my $target_name(sort keys %hashSeqRegion)
    {   
        # 参考序列
        my $seq        = $hashSeqRegion{$target_name}{'TargetSeq'}; 
        my $seq_length = $hashSeqRegion{$target_name}{'Length'}; 

        # 前后去掉的长度
        my $cut_F = (exists $hashSeqRegion{$target_name}{'PrimerF'}) ? length($hashSeqRegion{$target_name}{'PrimerF'}) : $cut_F_default;
        my $cut_R = (exists $hashSeqRegion{$target_name}{'PrimerR'}) ? length($hashSeqRegion{$target_name}{'PrimerR'}) : $cut_R_default;

        # 开始识别
        # CpG位点
        if($analysis_type =~ /CG/){
           while($seq =~ /cg/ig){
               my $position = pos($seq) - 1;
               if($position > $cut_F and $position <= ($seq_length - $cut_R)){
                  $hashAnalysisC{$target_name}{$position}{'Type'} = "CG";
               }
           }        
        }
        # CHG位点
        if($analysis_type =~ /CHG/){
           while($seq =~ /c[atc]g/ig){
               my $position = pos($seq) - 2;
               pos($seq)    = $position;
               if($position > $cut_F and $position <= ($seq_length - $cut_R)){
                  next if(exists $hashAnalysisC{$target_name}{$position} ); # 第二优先级CHG
                  $hashAnalysisC{$target_name}{$position}{'Type'} = "CHG";
               }
           }        
        }
        # CHHG位点
        if($analysis_type =~ /CHH/){
           while($seq =~ /c[atc][atc]g/ig){
               my $position = pos($seq) - 3;
               pos($seq)    = $position;
               if($position > $cut_F and $position <= ($seq_length - $cut_R)){
                  next if(exists $hashAnalysisC{$target_name}{$position} ); # 第三优先级CHHG
                  $hashAnalysisC{$target_name}{$position}{'Type'} = "CHH";
               }
           }        
        }
        # 上述3种位点外的所有位点
        if($analysis_type=~/AllOther/){
           while($seq =~ /c/ig  ){
               my $position = pos($seq);
               if($position > $cut_F and $position <= ($seq_length - $cut_R)){
                  next if(exists $hashAnalysisC{$target_name}{$position} ); # 第四优先级CN
                  my $nextallele = substr($seq, $position, 1);
                  next if($nextallele !~ /[atcg]/i);
                  $hashAnalysisC{$target_name}{$position}{'Type'} = "C$nextallele";
               }
           }        
        }
        # 上面需要分析的C之外的C位点
        while($seq =~ /c/ig){
            my $position = pos($seq);
            if($position > $cut_F  and  $position <= ($seq_length - $cut_R)){
                next if(exists $hashAnalysisC{$target_name}{$position} );           
                $hashOtherC{$target_name}{$position}++;
            }
        }         
    }

    # (2) 检查捕获序列是否有问题
    my @errors;
    foreach my $target_name(keys %hashSeqRegion)
    {
        next if(exists $hashAnalysisC{$target_name});
        push @errors, $target_name;
    }
    package::utils::is_continue2("    [Error] Target can not find AnalysisType C : @errors")  if(@errors > 0);
    
    # (3)确定每一个甲基化位置在参考基因组上的位置
    print "    check methylation point genome position\n";
    foreach my $target_name(sort keys %hashAnalysisC)
    {   
        # 参考序列情况
        my $chr           = $hashSeqRegion{$target_name}{'Chr'};
        my $tss           = $hashSeqRegion{$target_name}{'TSS'};
        my $mrna_strand   = $hashSeqRegion{$target_name}{'mRNA_Strand'};
        my $target_start  = $hashSeqRegion{$target_name}{'Start'};
        my $target_strand = $hashSeqRegion{$target_name}{'Target_Strand'};
        next if($target_start !~ /\d/); # 参考序列没有参考基因组信息

        foreach my $position(sort {$a <=> $b} keys %{$hashAnalysisC{$target_name}})
        {
            my $genome_position = ($target_strand eq '+') ? $target_start + $position - 1 : $target_start - $position + 1;
            my $distance_2tss   = '';
               $distance_2tss   = $genome_position - $tss if($mrna_strand eq '+');
               $distance_2tss   = $tss - $genome_position if($mrna_strand eq '-');

            $hashAnalysisC{$target_name}{$position}{'Chr'}            = $chr;
            $hashAnalysisC{$target_name}{$position}{'GenomePosition'} = $genome_position;
            $hashAnalysisC{$target_name}{$position}{'Distance2TSS'}   = $distance_2tss;
        }
    } 

    # (4)如果需要特殊标记，在Type符号上增加标记
    #    人的话，通常是增加高频rs编号
    if($special_mark ne '' and package::utils::is_file_ok($special_mark) == 1)
    {   
        print "    mark methylation point with special mark (rs id for human)\n";
        # 获取需要标记的基因组位置
        my %hashNeed;
        foreach my $target_name(keys %hashAnalysisC)
        {
            foreach my $position(keys %{$hashAnalysisC{$target_name}})
            {
                my $chr             = $hashAnalysisC{$target_name}{$position}{'Chr'};
                my $genome_position = $hashAnalysisC{$target_name}{$position}{'GenomePosition'};
                $hashNeed{"$chr,$genome_position"}++;
            }
        }
        # 读入
        my %hashSpecialMark;
        open SPECIALMARK, $special_mark;
        while(<SPECIALMARK>)
        {
            $_ =~ s/[\r\n]//g;
            my ($chr, $genome_position, $ref, $alt, $gene, $snp_id) = split /\t/, $_;
            next if(not exists $hashNeed{"$chr,$genome_position"}); # 减少内存消耗，全部读入的话，需要内存2.5G左右
            $hashSpecialMark{"$chr,$genome_position"} = $snp_id;
        }
        close SPECIALMARK;
        # 标注
        foreach my $target_name(keys %hashAnalysisC)
        {
            foreach my $position(keys %{$hashAnalysisC{$target_name}})
            {
                my $chr             = $hashAnalysisC{$target_name}{$position}{'Chr'};
                my $genome_position = $hashAnalysisC{$target_name}{$position}{'GenomePosition'};
                next if(not exists $hashSpecialMark{"$chr,$genome_position"});
                
                my $snp_id          = $hashSpecialMark{"$chr,$genome_position"};
                $hashAnalysisC{$target_name}{$position}{'Type'} .= "/$snp_id";
            }
        }
    }

    # (5)final 信息输出
    print "    output\n";
    my $analysis_c_point = "$report_dir/c.point.analysis.txt";
    my $other_c_point    = "$report_dir/c.point.other.txt";

    my @heads = ('Target', 
        'Position', 
        'Type', 
        'Chr', 
        'GenomePosition', 
        'Distance2TSS',
        );
    open ANALYSIS_C_POINT, ">$analysis_c_point";
    print ANALYSIS_C_POINT (join "\t", @heads) . "\n";
    foreach my $target_name(sort keys %hashAnalysisC)
    {   
        foreach my $position(sort {$a <=> $b} keys %{$hashAnalysisC{$target_name}})
        {
            my @values;
            foreach my $head(@heads)
            {
                my $value = '';
                $value = $hashAnalysisC{$target_name}{$position}{$head} if(exists $hashAnalysisC{$target_name}{$position}{$head});
                $value = $target_name if($head eq 'Target'); 
                $value = $position    if($head eq 'Position'); 
                push @values, $value;
            }
            print ANALYSIS_C_POINT (join "\t", @values) . "\n";
        }
    }
    close ANALYSIS_C_POINT;

    open OTHER_C_POINT, ">$other_c_point";
    print OTHER_C_POINT "Target\tPosition\n";
    foreach my $target_name(sort keys %hashOtherC)
    {   
        foreach my $position(sort {$a <=> $b} keys %{$hashOtherC{$target_name}})
        {
            print OTHER_C_POINT "$target_name\t$position\n";
        }
    }
    close OTHER_C_POINT;
}

sub check_fasta_position{
    my $hashConfig   = shift @_;
    my $hashPara     = shift @_;

    print "\n2: Check Fasta Genome Position and generate [seq.region] file\n";

    # 参数
    my $report_dir     = $hashConfig->{'Report'};
    my $target_fasta   = $hashConfig->{'TargetFasta'};
    my $primer         = $hashConfig->{'Primer'};
    my $blast_plus_dir = $hashPara->{'Soft'}{'BlastPlus'};
    my $species        = package::utils::get_species($hashConfig);
    my $genome_db      = exists $hashPara->{$species}{'Genome'}  ? $hashPara->{$species}{'Genome'}  : "";
    my $ref_gene       = exists $hashPara->{$species}{'RefGene'} ? $hashPara->{$species}{'RefGene'} : "";
    my $drop_UNnormal_chr = exists $hashConfig->{'Drop_UNnormal_chr'} ? $hashConfig->{'Drop_UNnormal_chr'} : "TRUE"; # 是否排除非常规染色体

    # 最终结果
    my $seq_region     = "$report_dir/seq.region"; 

    # （1）读取参考序列
    my %hashPrimer       = read_primer($primer);
    my %hashTargetFasta  = package::utils::read_fasta($target_fasta);

    # （2）与参考基因组比对，确定参考序列位置
    my $blast_out = "$report_dir/seq.fa.blast";
    if($genome_db ne "" and package::utils::is_file_ok($blast_out) == 0)
    {
        print "    Blast\n";
        system("$blast_plus_dir/blastn -task blastn -query $target_fasta -outfmt 6 -evalue 0.00001 -db $genome_db   -out $blast_out  -num_threads 4 -dust no");
    }else{
        print "    [Warnings] do not blast again as existing file $blast_out\n";
    }
    
    my %hashBlast = read_blast_result($blast_out, \%hashTargetFasta, $drop_UNnormal_chr); # 读取参考序列比对信息，锁定基因组位置

    # （3）读取转录本信息，确定基因最近的转录本
    my %hashGeneList = map{ ($hashBlast{$_}{'Gene_lowcase'}, 1)} keys %hashBlast;
    my %hashRefGene  = read_refgene($ref_gene, \%hashGeneList, $drop_UNnormal_chr); 

    # 确定距离基因最近的转录本
    my %hashGeneCatchPosition; # 基因捕获片段范围
    my %hashNearTranscript; # 最近的转录本
    foreach my $target_name(keys %hashBlast)
    {   
        my $gene_lowcase = $hashBlast{$target_name}{'Gene_lowcase'};
        my $start        = $hashBlast{$target_name}{'Start'};
        my $end          = $hashBlast{$target_name}{'End'};
        $hashGeneCatchPosition{$gene_lowcase}{$start}++;
        $hashGeneCatchPosition{$gene_lowcase}{$end}++;
    }
    foreach my $gene_lowcase(keys %hashGeneCatchPosition)
    {   
        my @positions = sort {$a<=>$b} keys %{$hashGeneCatchPosition{$gene_lowcase}};
        my ($start, $end) = ($positions[0], $positions[$#positions]);
        
        my $distance        = ""; # 最近的距离
        my $near_transcript = ""; # 最近的转录本
        foreach my $transcript (sort keys %{$hashRefGene{$gene_lowcase}})
        {
            
            my $tss = $hashRefGene{$gene_lowcase}{$transcript}{'TSS'};
            my $distance1 = abs($start - $tss);
            my $distance2 = abs($end   - $tss);

            if($distance eq "" or $distance1 < $distance or $distance2 < $distance )
            {
                $near_transcript = $transcript;
                $distance        = ($distance1 < $distance2) ? $distance1 : $distance2;
            }
        }
        next if($near_transcript eq "");

        $hashNearTranscript{$gene_lowcase} = $near_transcript;
    }

    # （4）根据最近转录本信息，确定相对TSS的距离
    foreach my $target_name(keys %hashTargetFasta)
    {
        my $gene_lowcase    = exists $hashBlast{$target_name}{'Gene_lowcase'} ? $hashBlast{$target_name}{'Gene_lowcase'} : "";
        my $near_transcript = exists $hashNearTranscript{$gene_lowcase}       ? $hashNearTranscript{$gene_lowcase}       : "";
        next if($gene_lowcase eq "" or $near_transcript eq "");

        my $target_start = $hashBlast{$target_name}{'Start'};
        my $tss          = $hashRefGene{$gene_lowcase}{$near_transcript}{'TSS'};
        my $mrna_strand  = $hashRefGene{$gene_lowcase}{$near_transcript}{'mRNA_Strand'};

        # 捕获序列起点距tss的距离
        $hashBlast{$target_name}{'Distance2TSS'} = ($mrna_strand eq '+') ? $target_start - $tss : $tss - $target_start;
        map{ $hashBlast{$target_name}{$_} = $hashRefGene{$gene_lowcase}{$near_transcript}{$_} if(not exists $hashBlast{$target_name}{$_}); } keys %{$hashRefGene{$gene_lowcase}{$near_transcript}};# 把转录本的信息放入Blast里
    }

    # (5) 输出
    open OUT, ">$seq_region";
    my @heads = ('Target',
        'Chr',
        'Gene',
        'mRNA',
        'mRNA_Strand',
        'TSS',
        'TES',
        'Target_Strand',
        'Start',
        'End',
        'Length',
        'Distance2TSS',
        'PrimerF',
        'PrimerR',
        'TargetSeq',
        );
    print OUT (join "\t", @heads) . "\n";
 
    foreach my $target_name(sort keys %hashTargetFasta)
    {   
        die "lost primer: $target_name" if(not exists $hashPrimer{$target_name}{'PrimerF'} or not exists $hashPrimer{$target_name}{'PrimerR'}); # 丢失引物信息
        my @values;
        foreach my $head(@heads)
        {
            my $value = '';
               $value = $target_name                            if($head eq 'Target');
               $value = $hashBlast{$target_name}{$head}         if(exists $hashBlast{$target_name}{$head});
               $value = $hashPrimer{$target_name}{$head}        if(exists $hashPrimer{$target_name}{$head});
               $value = length($hashTargetFasta{$target_name})  if($head eq 'Length');
               $value = $hashTargetFasta{$target_name}          if($head eq 'TargetSeq');
            push @values, $value; 
        }
        print OUT (join "\t", @values) . "\n";
    }
    close OUT;
}


# 读取引物
sub read_primer{
    my $primer = shift @_;

    my %hashPrimer;
    return %hashPrimer if(package::utils::is_file_ok($primer) == 0);

    open PRIMER, $primer;
    while(<PRIMER>){
        $_ =~ s/[\r\n]//g;
        next if($_!~/\w/);
        my ($target_name, $seq) = split /\t/, $_;
        die "PrimerError $_" if(!defined $seq );
        
        if($target_name =~ /F$/){
           $target_name =~ s/F$//;
           $hashPrimer{$target_name}{'PrimerF'}        = $seq;
           $hashPrimer{$target_name}{'PrimerF_length'} = length($seq); 

        }elsif($target_name =~ /R$/){
           $target_name =~ s/R$//;
           $hashPrimer{$target_name}{'PrimerR'}        = $seq;      
           $hashPrimer{$target_name}{'PrimerR_length'} = length($seq);        
        }
    }
    close PRIMER;

    return %hashPrimer;
}

# 读取转录本信息
sub read_refgene{
    my $ref_gene          = shift @_;
    my $hashGeneList      = shift @_;
    my $drop_UNnormal_chr = shift @_;

    my %hashRefGene;
    return %hashRefGene if(package::utils::is_file_ok($ref_gene) == 0);
       
    open REFGENE, $ref_gene;
    my $count = 0;
    while(<REFGENE>){
        $count++;
        my ($id, $mrna, $chr, $fangxiang, $g_s, $g_e, $cds_s, $cds_e, $n, $pos1, $pos2, $tp1, $gene, $tp2) = split /\t/, $_, 14;
        my $gene_lowcase = lc($gene);

        next if($drop_UNnormal_chr eq 'TRUE' and $chr =~ /_/); # 非常规染色体
        next if(!exists $hashGeneList->{$gene_lowcase} ); # 不是需要的基因

        $hashRefGene{$gene_lowcase}{"$mrna:$count"}{'Chr'}         = $chr; #"$mrna:$count",存在相同Mrna的情况
        $hashRefGene{$gene_lowcase}{"$mrna:$count"}{'Gene'}        = $gene;
        $hashRefGene{$gene_lowcase}{"$mrna:$count"}{'mRNA'}        = $mrna;
        $hashRefGene{$gene_lowcase}{"$mrna:$count"}{'mRNA_Strand'} = $fangxiang;
        $hashRefGene{$gene_lowcase}{"$mrna:$count"}{'TSS'}         = $fangxiang eq '+' ? $g_s : $g_e; # 转录起始终止
        $hashRefGene{$gene_lowcase}{"$mrna:$count"}{'TES'}         = $fangxiang eq '+' ? $g_e : $g_s;
    }
    close REFGENE;

    return %hashRefGene;
}

# 读取参考基因组比对结果
sub read_blast_result{
    my $blast_out         = shift @_;
    my $hashTargetFasta   = shift @_;
    my $drop_UNnormal_chr = shift @_;

    my %hashBlast;
    return %hashBlast if(package::utils::is_file_ok($blast_out) == 0);

    open BLAST, $blast_out;
    while(<BLAST>){
        $_=~ s/[\r\n]//g;
        my ($target_name, $chr, $identify_perc, $identify_len, $mismatch, $gap, $target_start, $target_end, $genome_start, $genome_end, $evalue, $score) = split /\t/, $_;
        next if($drop_UNnormal_chr eq 'TRUE' and $chr =~ /_/);

        my $target_length = length($hashTargetFasta->{$target_name});
        my ($gene_name, $tmp) = split /_/, $target_name;
        next if($mismatch > 0 or $gap > 0 or ($target_end - $target_start + 1) != $target_length); # 必须完美匹配

        $hashBlast{$target_name}{'Target'}        = $target_name;
        $hashBlast{$target_name}{'Chr'}           = $chr;
        $hashBlast{$target_name}{'Start'}         = $genome_start;
        $hashBlast{$target_name}{'End'}           = $genome_end;
        $hashBlast{$target_name}{'Target_Strand'} = ($genome_start < $genome_end) ? "+" : "-";
        $hashBlast{$target_name}{'Length'}        = $target_length;
        $hashBlast{$target_name}{'Gene_lowcase'}  = lc($gene_name);
    }
    close BLAST;
    return %hashBlast; 
}

 
# 建立甲基化序列库
sub build_methylation_db{
    my $hashConfig   = shift @_;
    my $hashPara     = shift @_;

    print "\n1:Build Methylation DB\n";
    # 参数
    my $report_dir   = $hashConfig->{'Report'};
    my $target_fasta = $hashConfig->{'TargetFasta'};
    my $blast_plus_dir = $hashPara->{'Soft'}{'BlastPlus'};

    # 数据库
    my $methylation_fasta_db_dir = "$report_dir/methylation_fasta_db";
    my $methylation_fasta        = "$methylation_fasta_db_dir/methylation.fa";
    package::utils::make_dir($methylation_fasta_db_dir); 

    # 甲基化处理
    open FASTA_IN, $target_fasta;
    open FASTA_ME, ">$methylation_fasta";
    while( my $target_name = <FASTA_IN>)
    {
        my $target_seq = <FASTA_IN>;
        $target_seq =~ s/C/T/ig; # C全部替换成T,甲基化
        print FASTA_ME "$target_name$target_seq";
    }
    close FASTA_IN;
    close FASTA_ME;

    # 建索引
    system("$blast_plus_dir/makeblastdb -in $methylation_fasta -dbtype nucl");

}


1
