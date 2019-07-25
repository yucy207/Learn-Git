package package::utils;
use strict;
use warnings;
$|=1;

# 读取配置文件
sub read_config{
    my $config_file = shift @_;
    my %hashConfig;

    open CONFIG, $config_file or die "[ERR] Loss Config File,$config_file\n";
    while(my $_=<CONFIG>){
        $_=~s/[\s\t\r\n]//g;
        next if($_ =~/^#/ or $_!~/\w/);
        my @datas = split /=/, $_;
        my $name  = $datas[0];
        my $value = $datas[1];
        next if(!defined $value);

        if($name eq 'case' or $name eq 'control')
        {
            $hashConfig{$name} .= "$value,";
        }
        elsif($name eq 'Force_Run' and $value =~ /\w/) 
        {
            map{ $hashConfig{'Force_Run'}{$_}++; } split /,/, $value;
        }
        elsif($name eq 'Gene' and $value =~ /\w/ and $datas[2] =~ /\w/) 
        {
            my @targets = split /,/, $datas[2];
            map{$hashConfig{"Gene"}{$value}{$_}++;}@targets;
        }
        else
        {
            $hashConfig{$name} = $value;
        }
    }
    close CONFIG;
    check_sample_duplicate(\%hashConfig); # 检查重名问题
    return %hashConfig;
}

sub check_sample_duplicate{
    my $hashConfig   = shift @_;
    my @samples      = get_sample($hashConfig, 'case', 'control'); # 样本
    my @ycontrols    = get_sample($hashConfig, 'Ycontrol');        # 内参对照

    # Ycontrol样本不能出现在case/control中
    foreach my $ycontrol(@ycontrols)
    {
        die "Ycontrol sample $ycontrol exists in case/control\n" if($ycontrol ~~ @samples);
    }
}

# 获取默认参数
sub get_para{
    my %hashPara;

    # 并行线程数配置
    $hashPara{'Process'}{'QC'}           = 20;
    $hashPara{'Process'}{'Mapping'}      = 20;
    $hashPara{'Process'}{'TargetTyping'} = 20;
    $hashPara{'Process'}{'MethylTyping'} = 20;
    $hashPara{'Process'}{'MappingBWA'}   = 20;
    $hashPara{'Process'}{'GATK_VCF'}     = 10;    

    # 软件配置
    $hashPara{'Soft'}{'Tmp'}        = "/home/tmp";   
    $hashPara{'Soft'}{'FastQC'}     = "/home/genesky/software/fastqc/0.11.5/fastqc";     
    $hashPara{'Soft'}{'FastqStat'}  = "/home/genesky/software/fastq_stat/fastq_stat"; 
    $hashPara{'Soft'}{'R'}          = "/home/genesky/software/r/3.5.1/bin/R"; 
    $hashPara{'Soft'}{'RLib'}       = "/home/genesky/software/r/3.5.1/lib64/R/library";
    $hashPara{'Soft'}{'BlastPlus'}  = "/home/genesky/software/blast+/2.7.1/bin"; 
    $hashPara{'Soft'}{'FastX'}      = "/home/genesky/software/fastx_toolkit/0.0.14/bin"; 
    $hashPara{'Soft'}{'FLASH'}      = "/home/genesky/software/flash/1.2.11/flash";

    $hashPara{'Soft'}{'BWA'}        = "/home/genesky/software/bwa/0.7.17/bwa";   
    $hashPara{'Soft'}{'SamTools'}   = "/home/genesky/software/samtools/1.9/samtools";
    $hashPara{'Soft'}{'sambamba'}   = "/home/genesky/software/sambamba/0.6.7/sambamba";
    $hashPara{'Soft'}{'cpulimit'}   = "/home/genesky/software/cpulimit/0.2/cpulimit";
    $hashPara{'Soft'}{'Picard'}     = "/home/genesky/software/picard/2.18.29/picard.jar";   
    $hashPara{'Soft'}{'GATK4_Loc'}  = "/home/genesky/software/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar"; 
    $hashPara{'Soft'}{'bcftools'}   = "/home/genesky/software/bcftools/1.9/bin/bcftools";   
    $hashPara{'Soft'}{'Bgzip'}      = "/home/genesky/software/htslib/1.9/bin/bgzip";   
    $hashPara{'Soft'}{'Java'}       = "/home/genesky/software/java/1.8.0_181/bin/java";
    $hashPara{'Soft'}{'AnnovarDIR'} = "/home/genesky/software/annovar/2018Apr16";

    # 参考基因组配置
    $hashPara{'Human_hg19'}{'Genome'}        = "/home/genesky/database/ucsc/hg19_modify/genome/blast+_idx/hg19_modify.fa";  # 参考基因组
    $hashPara{'Human_hg19'}{'RefGene'}       = "/home/genesky/database/ucsc/hg19_modify/gene/hg19_modify_refGene.txt";  # 基因信息
    $hashPara{'Human_hg19'}{'SpecialMark'}   = "/home/genesky/database/self_build_database/methylation_target_1000g_high_freq_snps/1000g_high_freq_snps.txt";  # 高频RS信息
    $hashPara{'Human_hg19'}{'AnnovarBuild'}  = "hg19";  # annovar 注释数据库前缀

    $hashPara{'Human_hg38'}{'Genome'}        = "/home/genesky/database/ucsc/hg38/genome/blast+_idx/hg38.fa";
    $hashPara{'Human_hg38'}{'RefGene'}       = "/home/genesky/database/ucsc/hg38/gene/hg38_refGene.txt";
    $hashPara{'Human_hg38'}{'AnnovarBuild'}  = "hg38";

    $hashPara{'Mouse_10'}{'Genome'}          = "/home/pub/database/Mouse/mm10/mm10.fa";
    $hashPara{'Mouse_10'}{'RefGene'}         = "/home/pub/database/Mouse/mm10/Annotation/mm10_refGene.txt";

    $hashPara{'Rat_rn6'}{'Genome'}           = "/home/pub/database/Rat/genome/rat_rn6.fa";
    $hashPara{'Rat_rn6'}{'RefGene'}          = "/home/pub/database/Rat/Annotation/rat_rn6_refGene.txt";

    return %hashPara;
}

# 展示配置内容
sub show_config{
    my $hashConfig = shift @_;
    my @case_samples     = get_sample($hashConfig, "case");
    my @control_samples  = get_sample($hashConfig, "control");
    my @Ycontrol_samples = get_sample($hashConfig, "Ycontrol");

    my $case_num      = @case_samples;
    my $control_num   = @control_samples;
    my $ycontrol_num  = @Ycontrol_samples;
    my $force_sample  = '';
       $force_sample  = join ",", keys %{$hashConfig->{'Force_Run'}} if (exists $hashConfig->{"Force_Run"});

    my $gender_file   = (exists $hashConfig->{"GenderFile"})   ? $hashConfig->{"GenderFile"}   : "[NA]";
    my $fastq_dir     = (exists $hashConfig->{"Fastq"})        ? $hashConfig->{"Fastq"}        : "[NA]";
    my $output_dir    = (exists $hashConfig->{"Output"})       ? $hashConfig->{"Output"}       : "[NA]";
    my $report_dir    = (exists $hashConfig->{"Report"})       ? $hashConfig->{"Report"}       : "[NA]";
    my $target_fasta  = (exists $hashConfig->{"TargetFasta"})  ? $hashConfig->{"TargetFasta"}  : "[NA]";
    my $primer        = (exists $hashConfig->{"Primer"})       ? $hashConfig->{"Primer"}       : "[NA]";
    my $analysis_type = (exists $hashConfig->{"AnalysisType"}) ? $hashConfig->{"AnalysisType"} : "[NA]";
    my $vcf_ploidy    = (exists $hashConfig->{"VCF_PLOIDY"})   ? $hashConfig->{"VCF_PLOIDY"}   : "[NA]";
    my $force         = (exists $hashConfig->{"Force"})        ? $hashConfig->{"Force"}        : "[NA]";
    my $species       = (exists $hashConfig->{"Species"})      ? $hashConfig->{"Species"}      : "[NA]";
 
    print "case:         $case_num\n";
    print "control:      $control_num\n";
    print "Ycontrol:     $ycontrol_num\n";
    print "Fastq:        $fastq_dir\n";
    print "Output:       $output_dir\n";
    print "Report:       $report_dir\n";
    print "TargetFasta:  $target_fasta\n";
    print "Primer:       $primer\n";
    print "AnalysisType: $analysis_type\n";
    print "VCF_PLOIDY:   $vcf_ploidy\n";
    print "Force:        $force\n";
    print "Species:      $species\n";
    system "echo -e '\\033[41;37mForce_Run : $force_sample \\033[0m'"  if($force_sample ne ''); # 警告 红底白字
}

# 获取当前项目分析的样本
sub get_sample{
    my $hashConfig = shift @_;
    my @need_types = @_;
    my @samples;

    foreach my $need_type(@need_types)
    {
        next if(!exists($hashConfig->{$need_type}));
        foreach my $sample(split /,/, $hashConfig->{$need_type})
        {
            push @samples, $sample if($sample=~/\w/);
        }
    }

    return @samples;
}

# 获取配对样本名
sub get_pair_sample{
    my $hashConfig = shift @_;
    my @samples = ();
    return @samples if(not exists $hashConfig->{'pairSamples'} or $hashConfig->{'pairSamples'} !~ /\w/);

    @samples = split /[,;]/, $hashConfig->{'pairSamples'};

    return @samples;
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 创建目录
sub make_dir{
    my @dir = @_;
    foreach(@dir){
        mkdir $_ if(not -e $_);
    }
    
}
# 是否继续
sub is_continue{
    print "\n[Option] Confirm and Start?[y/n]";
    my $input = get_input();
    exit if($input ne "y" and $input ne 'Y');
}

# 是否继续
sub is_continue2{
    my $info = shift @_;
    print "$info\n";
    print "\n[Option] Continue?[y/n]";
    my $input = get_input();
    exit if($input ne "y" and $input ne 'Y');
}

# 从屏幕获得输入
sub get_input{
    my $input=<STDIN>;
    $input=~ s/[\r\n]//g;
    return $input;
}

# 获取时间
sub get_time{
    my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
    $year+=1900;
    $mon+=1;
    return "$year-$mon-$day $hour:$min:$sec";
}

# 写运行日志
sub write_log{
    my $log_file = shift @_;
    my $title    = shift @_;
    my $hashCondition = shift @_;

    my $tab="##########";
    open LOG,">>$log_file";
    print LOG "$tab Start $title ".get_time()." $tab\n" if($title=~/\w/);
    # Finish
    my @finish_samples = exists $hashCondition->{'Finish'} ? sort keys %{$hashCondition->{'Finish'}} : ();
    print LOG 'Finish   : '. (join ",", @finish_samples) . "\n"; 
    # Good2Run
    my @good2run_samples = exists $hashCondition->{'Good2Run'} ? sort keys %{$hashCondition->{'Good2Run'}} : (); 
    print LOG 'Good2Run : '. (join ",", @good2run_samples) . "\n"; 
    # Error
    my @error_samples = exists $hashCondition->{'Error'} ? sort keys %{$hashCondition->{'Error'}} : ();
    print LOG 'Error    : '. (join ",", @error_samples) . "\n"; 

    print LOG "\n";
    close LOG;
}

# 写运行日志
sub write_log_simple{
    my $log_file = shift @_;
    my $title    = shift @_;

    my $tab="##########";
    open LOG,">>$log_file";
    print LOG "$tab Start $title ".get_time()." $tab\n" if($title=~/\w/);
    print LOG "\n";
    close LOG;
}

# 进度条， 根据输入向量计算
sub process_bar_array{
    my $process_name   = shift @_; # 需要分析的对象
    my $process_arrays = shift @_; # 总列表
    my $process_all_count = @$process_arrays;# 需要分析的总列表数量
    my ($pos) = grep $$process_arrays[$_] eq $process_name, 0..$process_all_count-1;
    my $process_count = $pos+ 1; # 当前对象的位置
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count; # 进度
    my $windows = 100; # 窗口宽度
    my $finished = $windows * int($process_perc) / 100; # 完成
    my $unfinished = $windows - $finished; # 未完成
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\n" if($process_count == $process_all_count);   
}

# 获取路径文件名
sub get_basename{
    my $path = shift @_;
    my $path_curf = File::Spec->rel2abs($path);
    my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
    return $file;
}

# 获取文件路径名
sub get_dirname{
    my $path = shift @_;
    my $path_curf = File::Spec->rel2abs($path);
    my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
    return $dirs;
}

# 检查VCF标题是否缺失样本
sub check_vcf_sample_isOK{
    my $vcf_file = shift @_;
    my $samples = shift @_;
    my $is_gzip = ($vcf_file=~/\.gz$/) ? 1 : 0;# VCF是否是压缩的

    my @losts;
    open VCF, "gzip -cd $vcf_file|" if($is_gzip == 1);
    open VCF, $vcf_file             if($is_gzip == 0);
    while(<VCF>)
    {
        next if($_=~/^##/);
        if($_=~/^#CHROM/)
        {
            $_=~s/[\r\n]//g;
            my @vcf_heads = split /\t/, $_;
            foreach my $sample(@$samples)
            {
                next if($sample ~~ @vcf_heads);
                push @losts, $sample;
            }
            last;
        }
    }
    close VCF;
    return @losts;
}

# 获取物种信息
sub get_species{
    my $hashConfig = shift @_;
    my $species = (exists $hashConfig->{'Species'}) ? $hashConfig->{'Species'} : 'Lost Species';
    return $species;
}

# 提取fasta序列
sub read_fasta{
    my $fasta = shift @_;
    my %hashFasta;
    print "read $fasta ... ";
    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    while(my $inSeq = $IN->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        $hashFasta{$id}= $seq;
    }
    print "OK\n";
    return %hashFasta;
}

# 读取甲基化参考序列信息
sub read_seq_region{
    my $seq_region = shift @_;
    my %hashSeqRegion;

    open SEQ_REGION, $seq_region;
    my $line1 = <SEQ_REGION>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    while(<SEQ_REGION>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        my %hashTmp = map{ ($heads[$_], $datas[$_]) } (0..$#heads);
        my $target = $hashTmp{'Target'};

        map{ $hashSeqRegion{$target}{$_} = $hashTmp{$_}; }keys %hashTmp;
    }
    close SEQ_REGION;

    return %hashSeqRegion;
}

sub region_anno{
    my $hashCPos      = shift @_;
    my $hashConfig    = shift @_;
    my $seq_region    = "$hashConfig->{'Report'}/seq.region";
    my %hashmark      = read_mark($hashConfig);# 需要显著标出的位置  
    my %hashSeqRegion = read_seq_region($seq_region);
    # 位点注释
    my %hashanno;
    foreach my $target(keys %$hashCPos){
        my $start  = exists($hashSeqRegion{$target}{'Start'}) ? $hashSeqRegion{$target}{'Start'}:"";
        my $chr    = exists($hashSeqRegion{$target}{'Chr'}) ? $hashSeqRegion{$target}{'Chr'}:"";
        my $strand = exists($hashSeqRegion{$target}{'Target_Strand'}) ? $hashSeqRegion{$target}{'Target_Strand'}:"";
        my $TSS    = exists($hashSeqRegion{$target}{'TSS'}) ? $hashSeqRegion{$target}{'TSS'}:"";
        my $mRNA_Strand = exists($hashSeqRegion{$target}{'mRNA_Strand'}) ? $hashSeqRegion{$target}{'mRNA_Strand'}:"";
        foreach my $position(keys %{$hashCPos->{$target}}){
            if($chr=~/\w/){
               $hashanno{$target}{$position}{'Chr'} = $chr;
               $hashanno{$target}{$position}{'GenomePosition'} = $strand eq '+'? $start + $position - 1 : $start - $position + 1;
               if($TSS=~/\w/){
                  my $genomeposition = $hashanno{$target}{$position}{'GenomePosition'};
                  $hashanno{$target}{$position}{'Distance2TSS'} = $mRNA_Strand eq '+' ? $genomeposition - $TSS : $TSS - $genomeposition;
               }
            }
            # 标注
            my $genomeposition=exists($hashanno{$target}{$position}{'GenomePosition'}) ? $hashanno{$target}{$position}{'GenomePosition'}:"";
            $hashanno{$target}{$position}{'Mark'}           = (exists($hashmark{'Target'}{$target}{$position}) or exists($hashmark{'Chr'}{$chr}{$genomeposition})) ? "checked":'';
            $hashanno{$target}{$position}{'Chr'}            = "" if(!exists($hashanno{$target}{$position}{'Chr'}));
            $hashanno{$target}{$position}{'GenomePosition'} = "" if(!exists($hashanno{$target}{$position}{'GenomePosition'}));
            $hashanno{$target}{$position}{'Distance2TSS'}   = "" if(!exists($hashanno{$target}{$position}{'Distance2TSS'}));
            $hashanno{$target}{$position}{'Target_Strand'}  = $strand;
        }
    }
    return %hashanno;
}

sub read_mark{
    my $hashConfig = shift @_;
    my $mark_file  = (exists($hashConfig->{'Mark'}) and -f $hashConfig->{'Mark'}) ? $hashConfig->{'Mark'}:"";
    my %hashmark;
    return %hashmark if($mark_file!~/\w/);
    open MARK, $mark_file;
    while(<MARK>){
        $_=~s/[\r\n\s]//g;
        next if($_!~/\w/);
        my($type,$name,$position) = split /=/,$_;
        $hashmark{$type}{$name}{$position} = 0;
    }
    close MARK;
    return %hashmark;   
}

# 表格格式模版
sub sheet_format{
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

    $format{'red'} = $workbook->add_format();
    $format{'red'} ->set_align('center');
    $format{'red'} ->set_align('vcenter');
    $format{'red'} ->set_size(12);
    $format{'red'} ->set_font("Times New Roman");
    $format{'red'} ->set_bg_color("#ff0000");
    $format{'red'} ->set_border();
    
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
    $format{'orange'} ->set_align('center');
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();

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
1
