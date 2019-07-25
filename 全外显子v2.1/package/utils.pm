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
        my ($name, $value) = split /=/, $_;
        next if(!defined $value);

        if($name eq 'case' or $name eq 'control')
        {
            $hashConfig{$name} .= "$value,";
        }
        elsif($name eq 'Force_Run' and $value =~ /\w/) 
        {
            map{ $hashConfig{'Force_Run'}{$_}++; } split /,/, $value;
        }
        else
        {
            $hashConfig{$name} = $value;
        }
    }
    close CONFIG;
    return %hashConfig;
}

# 获取默认参数
sub get_para{
    my %hashPara;

    # 并行线程数配置
    $hashPara{'Process'}{'QC'}         = 20;
    $hashPara{'Process'}{'Mapping'}    = 20;
    $hashPara{'Process'}{'GATK_GVCF'}  = 20;
    $hashPara{'Process'}{'GATK_VCF'}   = 10;
    $hashPara{'Process'}{'Annotation'} = 20;

    # 软件配置
    $hashPara{'Soft'}{'Tmp'}        = "/home/tmp";   
    $hashPara{'Soft'}{'FastQC'}     = "/home/genesky/software/fastqc/0.11.5/fastqc";     
    $hashPara{'Soft'}{'FastqStat'}  = "/home/genesky/software/fastq_stat/fastq_stat"; 
    $hashPara{'Soft'}{'BWA'}        = "/home/genesky/software/bwa/0.7.17/bwa";   
    $hashPara{'Soft'}{'SamTools'}   = "/home/genesky/software/samtools/1.9/samtools";
    $hashPara{'Soft'}{'sambamba'}   = "/home/genesky/software/sambamba/0.6.7/sambamba";
    $hashPara{'Soft'}{'cpulimit'}   = "/home/genesky/software/cpulimit/0.2/cpulimit";  
    $hashPara{'Soft'}{'Picard'}     = "/home/genesky/software/picard/2.18.29/picard.jar"; 
    $hashPara{'Soft'}{'GATK3'}      = "/home/genesky/software/gatk/3.5/GenomeAnalysisTK.jar";   
    $hashPara{'Soft'}{'GATK4_Loc'}  = "/home/genesky/software/gatk/4.1.2.0/gatk-package-4.1.2.0-local.jar";   
    $hashPara{'Soft'}{'GATK4_Spa'}  = "/home/genesky/software/gatk/4.1.2.0/gatk-package-4.1.2.0-spark.jar";
    $hashPara{'Soft'}{'bcftools'}   = "/home/genesky/software/bcftools/1.9/bin/bcftools";   
    $hashPara{'Soft'}{'Java'}       = "/home/genesky/software/java/1.8.0_181/bin/java"; 
    $hashPara{'Soft'}{'R'}          = "/home/genesky/software/r/3.5.1/bin/R"; 
    $hashPara{'Soft'}{'RLib'}       = "/home/genesky/software/r/3.5.1/lib64/R/library"; 
    $hashPara{'Soft'}{'AnnovarDir'} = "/home/genesky/software/annovar/2018Apr16"; 
    $hashPara{'Soft'}{'MegaBlast'}  = "/home/genesky/software/blast/20110130/bin/megablast"; 
    $hashPara{'Soft'}{'Python27'}   = "/home/genesky/software/python/2.7.13/bin/python"; 
    $hashPara{'Soft'}{'InterVar'}   = "/home/genesky/software/intervar/2.1.2"; 
    $hashPara{'Soft'}{'snpEFF'}     = "/home/genesky/software/snpeff/4_3t/snpEff/snpEff.jar";
    $hashPara{'Soft'}{'plink'}      = "/home/genesky/software/plink/1.90_beta/plink";
    $hashPara{'Soft'}{'king'}       = "/home/genesky/software/king/2.0/king";
    $hashPara{'Soft'}{'fasttree'}   = "/home/genesky/software/fasttree/2.1.10/FastTree";

    
    # 参考基因组配置
    $hashPara{'Human_hg19'}{'Genome'}            = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.fa";  # 参考基因组
    $hashPara{'Human_hg19'}{'blast_idx'}         = "/home/genesky/database/ucsc/hg19_modify/genome/blast+_idx/hg19_modify.fa";  
    $hashPara{'Human_hg19'}{'bwa_idx'}           = "/home/genesky/database/ucsc/hg19_modify/genome/bwa_idx/hg19_modify.fa";  
    $hashPara{'Human_hg19'}{'samtools_idx'}      = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.fa";  
    $hashPara{'Human_hg19'}{'DBsnp'}             = "/home/genesky/database/gatk/hg19_modify/dbsnp/dbsnp_138.hg19.vcf";  # DBSNP 用于GATK校正
    $hashPara{'Human_hg19'}{'InDel'}             = "/home/genesky/database/gatk/hg19_modify/1000G_phase1_indels/1000G_phase1.indels.hg19.vcf";  # 1000g indel 用于GATK校正
    $hashPara{'Human_hg19'}{'InDelGold'}         = "/home/genesky/database/gatk/hg19_modify/mills_and_1000g_gold_standard_indels/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"; # 1000g indel 用于GATK校正
    $hashPara{'Human_hg19'}{'VQSR_HAPMAP'}       = "/home/genesky/database/gatk/hg19_modify/hapmap/hapmap_3.3.hg19.sites.vcf.gz";                        # VQSR
    $hashPara{'Human_hg19'}{'VQSR_1000G_OMNI'}   = "/home/genesky/database/gatk/hg19_modify/1000G_omni/1000G_omni2.5.hg19.sites.vcf.gz";                     # VQSR
    $hashPara{'Human_hg19'}{'VQSR_1000G_SNP_HC'} = "/home/genesky/database/gatk/hg19_modify/1000G_phase1_hc/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz"; # VQSR    
    $hashPara{'Human_hg19'}{'AnnovarBuild'}      = "hg19";  # annovar 注释数据库前缀
    $hashPara{'Human_hg19'}{'snpEFF_config'}     = "/home/genesky/database/snpeff/4_3t/snpEff.config";    # snpEFF 数据库路径
    $hashPara{'Human_hg19'}{'snpEFF_build'}      = "hg19";    # snpEFF 数据库前缀
    $hashPara{'Human_hg19'}{'exclude_96_bed'}    = "/home/genesky/database/self_build_database/chip_capture_bed/hg19_exclude_96_sites.bed";    # GVCF生成时排除96位点所在区域

    $hashPara{'Human_hg19_HLA'}{'Genome'}        = "/home/genesky/database/ucsc/hg19_chr6/genome/hg19_chr6.fa";   
    $hashPara{'Human_hg19_HLA'}{'blast_idx'}     = "/home/genesky/database/ucsc/hg19_chr6/genome/blast+_idx/hg19_chr6.fa";   
    $hashPara{'Human_hg19_HLA'}{'bwa_idx'}       = "/home/genesky/database/ucsc/hg19_chr6/genome/bwa_idx/hg19_chr6.fa";   
    $hashPara{'Human_hg19_HLA'}{'samtools_idx'}  = "/home/genesky/database/ucsc/hg19_chr6/genome/hg19_chr6.fa";   
    $hashPara{'Human_hg19_HLA'}{'DBsnp'}         = "/home/genesky/database/gatk/hg19_modify/dbsnp/dbsnp_138.hg19.vcf";   
    $hashPara{'Human_hg19_HLA'}{'InDel'}         = "/home/genesky/database/gatk/hg19_modify/1000G_phase1_indels/1000G_phase1.indels.hg19.vcf";   
    $hashPara{'Human_hg19_HLA'}{'InDelGold'}     = "/home/genesky/database/gatk/hg19_modify/mills_and_1000g_gold_standard_indels/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"; 
    $hashPara{'Human_hg19_HLA'}{'AnnovarBuild'}  = "hg19";   
    $hashPara{'Human_hg19_HLA'}{'snpEFF_config'} = "/home/genesky/database/snpeff/4_3t/snpEff.config"; 
    $hashPara{'Human_hg19_HLA'}{'snpEFF_build'}  = "hg19";
    $hashPara{'Human_hg19_HLA'}{'exclude_96_bed'}= "/home/genesky/database/self_build_database/chip_capture_bed/hg19_exclude_96_sites.bed";    # GVCF生成时排除96位点所在区域

    $hashPara{'Human_hg38'}{'Genome'}            = "/home/genesky/database/ucsc/hg38_modify/genome/hg38_modify.fa";   
    $hashPara{'Human_hg38'}{'blast_idx'}         = "/home/genesky/database/ucsc/hg38_modify/genome/blast+_idx/hg38_modify.fa";   
    $hashPara{'Human_hg38'}{'bwa_idx'}           = "/home/genesky/database/ucsc/hg38_modify/genome/bwa_idx/hg38_modify.fa";   
    $hashPara{'Human_hg38'}{'samtools_idx'}      = "/home/genesky/database/ucsc/hg38_modify/genome/hg38_modify.fa";   
    $hashPara{'Human_hg38'}{'DBsnp'}             = "/home/genesky/database/gatk/hg38/dbsnp/dbsnp_146.hg38.vcf.gz";     
    $hashPara{'Human_hg38'}{'InDelGold'}         = "/home/genesky/database/gatk/hg38/mills_and_1000g_gold_standard_indels/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"; 
    $hashPara{'Human_hg38'}{'VQSR_HAPMAP'}       = "/home/genesky/database/gatk/hg38/hapmap/hapmap_3.3.hg38.vcf.gz";                        # VQSR
    $hashPara{'Human_hg38'}{'VQSR_1000G_OMNI'}   = "/home/genesky/database/gatk/hg38/1000G_omni/1000G_omni2.5.hg38.vcf.gz";                     # VQSR
    $hashPara{'Human_hg38'}{'VQSR_1000G_SNP_HC'} = "/home/genesky/database/gatk/hg38/1000G_phase1_hc/1000G_phase1.snps.high_confidence.hg38.vcf.gz"; # VQSR    
    $hashPara{'Human_hg38'}{'AnnovarBuild'}      = "hg38";   
    $hashPara{'Human_hg38'}{'snpEFF_config'}     = "/home/genesky/database/snpeff/4_3t/snpEff.config";     
    $hashPara{'Human_hg38'}{'snpEFF_build'}      = "hg38";  

    return %hashPara;
}

# 获取注释列表
sub get_annotation_list{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    my $species     = package::utils::get_species($hashConfig);
    my $report_dir  = $hashConfig->{'Report'}; 
    my $library_dir = "$report_dir/library";
    my $anno_log    = "$library_dir/anno.log";
    my %hashKeyWord = package::annotation_db::get_annotation_key_word($hashPara, $hashConfig); # 获取关键字哈希

    open ANNOLOG, ">$anno_log";
    my %hashAnnotation_all = package::annotation_db::get_annotation(); # 所有注释数据库
    my %hashAnnotation; # 当前物种最终可以注释的数据库
    foreach my $anno_model(sort keys %hashAnnotation_all)
    {
        my $model_species = $hashAnnotation_all{$anno_model}{'Species'};

        # 1 物种检测
        next if($model_species ne '[all]' and $model_species !~/\[$species\]/); # 过滤，非当前物种

        # 2 分析必要文件检测
        my $result = check_annotation_necessary_file($hashAnnotation_all{$anno_model}, \%hashKeyWord, $species);

        # 结论
        if($result eq '') # 正常
        {
            print ANNOLOG "annotation $anno_model is good to run\n";
            $hashAnnotation{$anno_model} = $hashAnnotation_all{$anno_model};
        }
        else # 有错误
        {
            print ANNOLOG "annotation $anno_model is Error: $result\n";
        }
    }
    return %hashAnnotation; 
}

# 检测必要文件
sub check_annotation_necessary_file{
    my $hashAnnotation = shift @_;
    my $hashKeyWord    = shift @_;
    my $species        = shift @_;

    my $isOK = 1;
    my @errors = ();
    foreach my $file_index(sort keys %{$hashAnnotation->{'NecessaryFile'}})
    {
        my $anno_dir = (exists $hashAnnotation->{'DBDir'}{$species}) ? $hashAnnotation->{'DBDir'}{$species}."/": "";
        my %hashInfo = split_annotation_info($hashAnnotation->{'NecessaryFile'}{$file_index});
        my $prefix = (exists($hashKeyWord->{$hashInfo{'prefix'}})) ? $hashKeyWord->{$hashInfo{'prefix'}} : $hashInfo{'prefix'};
        my $suffix = $hashInfo{'suffix'};
        my $need_file = $anno_dir.$prefix.$suffix;
        my $condition = (package::utils::is_file_ok($need_file)==1) ? 'OK' : 'Lost';
        next if($condition eq 'OK');
        push @errors, "Lost $need_file";
    }
    
    return (join ",", @errors);
}

# 拆分注释信息
sub split_annotation_info{
    my $info = shift @_;
    my %hashInfo;
    foreach(split /\t/, $info)
    {
        my ($name, $value) = split /=/, $_;
        $hashInfo{$name} = $value;
    }
    return %hashInfo;
}

# 展示配置内容
sub show_config{
    my $hashConfig = shift @_;
    my @case_samples    = get_sample($hashConfig, "case");
    my @control_samples = get_sample($hashConfig, "control");

    my $case_num      = @case_samples;
    my $control_num   = @control_samples;
    my $force_sample  = '';
       $force_sample  = join ",", keys %{$hashConfig->{'Force_Run'}} if (exists $hashConfig->{"Force_Run"});
    my $gender_file   = (exists $hashConfig->{"GenderFile"}) ? $hashConfig->{"GenderFile"} : "[NA]";
    my $fastq_dir     = (exists $hashConfig->{"Fastq"})      ? $hashConfig->{"Fastq"}      : "[NA]";
    my $output_dir    = (exists $hashConfig->{"Output"})     ? $hashConfig->{"Output"}     : "[NA]";
    my $report_dir    = (exists $hashConfig->{"Report"})     ? $hashConfig->{"Report"}     : "[NA]";
    my $target_bed    = (exists $hashConfig->{"TargetBed"})  ? $hashConfig->{"TargetBed"}  : "[NA]";
    my $realign_bed   = (exists $hashConfig->{"RealignBed"}) ? $hashConfig->{"RealignBed"} : "[NA]";
    my $species       = (exists $hashConfig->{"Species"})    ? $hashConfig->{"Species"}    : "[NA]";
 
    print "case:       $case_num\n";
    print "control:    $control_num\n";
    print "GenderFile: $gender_file\n";
    print "Fastq:      $fastq_dir\n";
    print "Output:     $output_dir\n";
    print "Report:     $report_dir\n";
    print "TargetBed:  $target_bed\n";
    print "RealignBed: $realign_bed\n";
    print "Species:    $species\n";
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
# 检查项目是否进行了强制运行
sub is_force_sample{
    my $hashConfig = shift @_;
    my $is_force_sample = 0;
    if(exists $hashConfig->{'Force_Run'})
    {
        my @force_samples = keys %{$hashConfig->{'Force_Run'}};
        $is_force_sample = 1 if(@force_samples > 0);
    }
    return $is_force_sample;
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
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
}
# 是否继续
sub is_continue{
    print "\n[Option] Confirm and Start?[y/n]";
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
    print "\n";
    return %hashFasta;
}

# 读取通用SNP注释
sub read_general_snp_annotation{
    my $hashAnnotation = shift @_;
    my $hashAnnotationList = shift @_;
    my $hashKeyWord        = shift @_;
    my $need_model         = shift @_;

    foreach my $out_index(keys %{$hashAnnotationList->{'Output'}})
    {
        my %hashFile_Info = split_annotation_info($hashAnnotationList->{'Output'}{$out_index}{'File'}); # 文件信息
        my %hashNeed_Info = split_general_snp_need_info($hashAnnotationList->{'Output'}{$out_index}{'Need_Info'}); # 需要的信息坐标
        # 文件检测
        my $prefix = (exists($hashKeyWord->{$hashFile_Info{'prefix'}})) ? $hashKeyWord->{$hashFile_Info{'prefix'}} : $hashFile_Info{'prefix'};
        my $suffix = $hashFile_Info{'suffix'};
        my $output_file = $prefix.$suffix; # 当前注释内容的输出文件
        if(is_file_ok($output_file)==0)
        {
            print "Lost $need_model : $output_file\n";
            next;
        }

        # 文件信息读取
        print "Process $need_model : $output_file ... ";
        open OUTPUTFILE, $output_file;
        while(<OUTPUTFILE>)
        {
            $_=~s/[\r\n]//g;
            next if($_=~/^#/);
            my @datas = split /\t/, $_;
            foreach my $info_name(keys %hashNeed_Info)
            {
                my $chr   = $datas[$hashNeed_Info{$info_name}{'chr_index'}];
                my $start = $datas[$hashNeed_Info{$info_name}{'start_index'}];
                my $end   = $datas[$hashNeed_Info{$info_name}{'end_index'}];
                my $ref   = $datas[$hashNeed_Info{$info_name}{'ref_index'}];
                my $alt   = $datas[$hashNeed_Info{$info_name}{'alt_index'}];
                my $snp_title   = "$chr|$start|$end|$ref|$alt"; # SNP表头
                my $value_index = $hashNeed_Info{$info_name}{'value_index'}; # 取值index
                my $value       = '';

                if($value_index=~/\//) # 复合式结果，需要拆分提取
                {
                    my ($index1, $split_mark, $index2) = split /\//, $value_index;
                    my $value1 = exists($datas[$index1]) ? $datas[$index1] : '' ;
                    my @datas2 = split /$split_mark/, $value1;

                    my @values = map{ my $value_tmp = exists($datas2[$_]) ? $datas2[$_] : '' ; $value_tmp} (split /,/, $index2); # 应对需要获得多个结果的情况
                    $value = join "|", @values;
                }
                else
                {   
                    my @values = map{ my $value_tmp = exists($datas[$_]) ? $datas[$_] : '' ; $value_tmp} (split /,/, $value_index); # 应对需要获得多个结果的情况
                    $value = join "|", @values;
                }
                
                $value =~ s/\\x2c/,/g; # annovar 注释内容会把‘，’或‘#’替换为特殊字符，防止冲突，所以使用时需要再替换回来
                $value =~ s/\\x23/#/g;

                $hashAnnotation->{$snp_title}{$info_name} = $value if($value  ne "");
            }
        }
        print "OK\n";
    }
 
}

# 读取通用GENE注释
sub read_general_gene_annotation{
    my $hashAnnotation     = shift @_;
    my $hashAnnotationList = shift @_;
    my $hashKeyWord        = shift @_;
    my $need_model         = shift @_;
    my $species            = shift @_;
    foreach my $out_index(keys %{$hashAnnotationList->{'NecessaryFile'}})
    {
        my %hashFile_Info = split_annotation_info($hashAnnotationList->{'NecessaryFile'}{$out_index}); # 文件信息
        my %hashNeed_Info = split_general_gene_need_info($hashAnnotationList->{'Output'}{$out_index}{'Need_Info'}); # 需要的信息坐标

        # 文件检测
        my $anno_dir = (exists $hashAnnotationList->{'DBDir'}{$species}) ? $hashAnnotationList->{'DBDir'}{$species}."/": "";
        my $prefix = (exists($hashKeyWord->{$hashFile_Info{'prefix'}})) ? $hashKeyWord->{$hashFile_Info{'prefix'}} : $hashFile_Info{'prefix'};
        my $suffix = $hashFile_Info{'suffix'};
        my $output_file = $anno_dir.$prefix.$suffix; # 当前注释内容的输出文件
        if(is_file_ok($output_file)==0)
        {
            print "Lost $need_model : $output_file\n";
            next;
        }

        # 文件信息读取
        print "Process $need_model : $output_file ... ";
        open OUTPUTFILE, $output_file;
        while(<OUTPUTFILE>)
        {
            $_=~s/[\r\n]//g;
            next if($_=~/^#/);
            my @datas = split /\t/, $_;
            foreach my $info_name(keys %hashNeed_Info)
            {
                my $gene        = $datas[$hashNeed_Info{$info_name}{'gene_index'}];
                my $value_index = $hashNeed_Info{$info_name}{'value_index'}; # 取值index
                my $value       = '';

                if($value_index=~/\//) # 复合式结果，需要拆分提取
                {
                    my ($index1, $split_mark, $index2) = split /\//, $value_index;
                    my $value1 = exists($datas[$index1]) ? $datas[$index1] : '' ;
                    my @datas2 = split /$split_mark/, $value1;

                    my @values = map{ my $value_tmp = exists($datas2[$_]) ? $datas2[$_] : '' ; $value_tmp} (split /,/, $index2); # 应对需要获得多个结果的情况
                    $value = join "|", @values;
                }
                else
                {   
                    my @values = map{ my $value_tmp = exists($datas[$_]) ? $datas[$_] : '' ; $value_tmp} (split /,/, $value_index); # 应对需要获得多个结果的情况
                    $value = join "|", @values;
                }

                $hashAnnotation->{$gene}{$info_name} = $value if($value ne "");
            }

        }
        print "OK\n";
    }
 
}

sub split_general_snp_need_info{
    my $need_info = shift @_;
    my %hashNeed_Info;
    foreach my $info_detail(split /\t/, $need_info)
    {
        my ($info_name, $value_index, $snp_title_index) = split /\=/, $info_detail;
        my ($chr_index, $start_index, $end_index, $ref_index, $alt_index) = split /,/, $snp_title_index;
        $hashNeed_Info{$info_name}{'value_index'} = $value_index;
        $hashNeed_Info{$info_name}{'chr_index'}   = $chr_index;
        $hashNeed_Info{$info_name}{'start_index'} = $start_index;
        $hashNeed_Info{$info_name}{'end_index'}   = $end_index;
        $hashNeed_Info{$info_name}{'ref_index'}   = $ref_index;
        $hashNeed_Info{$info_name}{'alt_index'}   = $alt_index;
    }
    return %hashNeed_Info;
}

sub split_general_gene_need_info{
    my $need_info = shift @_;
    my %hashNeed_Info;
    foreach my $info_detail(split /\t/, $need_info)
    {
        my ($info_name, $value_index, $gene_index) = split /\=/, $info_detail;
        $hashNeed_Info{$info_name}{'value_index'} = $value_index;
        $hashNeed_Info{$info_name}{'gene_index'}   = $gene_index;
    }
    return %hashNeed_Info;
}

1