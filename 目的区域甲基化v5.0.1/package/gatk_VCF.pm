package package::gatk_VCF;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start GATK VCF ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::gatk_VCF($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in GATK VCF run para, please check\n" if($isOK == 0);

    # 路径
    my $report_dir = $hashConfig->{'Report'};
    my $output_dir = $hashConfig->{'Output'};
    my $vcf_dir    = "$report_dir/vcf";
    my $gvcf_dir   = "$report_dir/vcf/gvcf";
    package::utils::make_dir($vcf_dir);
    package::utils::make_dir($gvcf_dir);

    # 文件
    my $genome_file         = "$report_dir/TargetFastaGenome/seq.fa";            # 参考基因组
    my $target_realign_bed  = "$report_dir/TargetFastaGenome/targetRealign.bed"; # 目标区域

    # 软件
    my $GATK4_Loc = $hashPara->{"Soft"}{"GATK4_Loc"};
    my $Java      = $hashPara->{"Soft"}{"Java"};
    my $Tmp       = $hashPara->{"Soft"}{"Tmp"};
    my $bgzip     = $hashPara->{"Soft"}{"Bgzip"};
    my $bcftools  = $hashPara->{"Soft"}{"bcftools"};

    # 参数
    my $vcf_ploidy= (exists($hashConfig->{"VCF_PLOIDY"}) and $hashConfig->{"VCF_PLOIDY"}=~/^\d+$/) ? $hashConfig->{"VCF_PLOIDY"} : 2;

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control');
    my %hashCondition; 
    foreach my $sample(@samples)
    {
        my $bam_index = "$output_dir/$sample/$sample\_final.bai";
        my $g_vcf  = "$gvcf_dir/$sample.g.vcf.idx";
        # 已完成该步骤
        if(package::utils::is_file_ok($g_vcf) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$g_vcf";
            next;            
        }
        # 原始数据没问题
        if(package::utils::is_file_ok($bam_index)){
            $hashCondition{"Good2Run"}{$sample} = "$bam_index";
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$bam_index";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "GATK_VCF", \%hashCondition); 
    die "Sample lost final bam :" . (join ",", sort keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));

    # 执行比对
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0)
    {   
        my $threshold = exists $hashPara->{"Process"}{"GATK_VCF"} ? $hashPara->{"Process"}{"GATK_VCF"} : 10;
           $threshold = $hashConfig->{"Process_GATK_VCF"} if(exists $hashConfig->{"Process_GATK_VCF"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::utils::process_bar_array($sample, \@sample_runs)});# 进度条
        foreach my $sample(@sample_runs)
        {
            $pm->start($sample) and next;
            my $bam_final = "$output_dir/$sample/$sample\_final.bam";
            system "$Java -jar -Xmx5g $GATK4_Loc HaplotypeCaller -I $bam_final -O $gvcf_dir/$sample.g.vcf.gz -ERC GVCF -R $genome_file --min-base-quality-score 20 -L $target_realign_bed --tmp-dir $Tmp --kmer-size 10  --kmer-size 20 --kmer-size 25 --kmer-size 30 --kmer-size 35 --kmer-size 40 --max-num-haplotypes-in-population 512 --max-assembly-region-size 200";
            $pm->finish;    
        }
        $pm->wait_all_children;

        #####
        # GVCF全部准备完成后，合并生成突变vcf
        #####
        print "Combining *.g.vcf.gz \n";
        my @okGVCF;   # gvcf生成的样本参数
        my @lostGVCF; # 样本丢失vcf或者无法生成gvcf
        foreach my $sample(@samples){
            my $gvcf = "$gvcf_dir/$sample.g.vcf.gz";
            if(package::utils::is_file_ok($gvcf)==1){
                push @okGVCF, "--variant $gvcf";
            }else{
                push @lostGVCF, $sample;
            }
        }
        print "Sample Lost GVCF:\t@lostGVCF\n" if(@lostGVCF > 0);
        return if(@okGVCF==0);
        my $gvcf_para     = join " ",@okGVCF;
        my $gvcf_combined = "$vcf_dir/allSample.combine.g.vcf.gz";
        system ("$Java  -jar -Xmx5g $GATK4_Loc CombineGVCFs  -L $target_realign_bed -R $genome_file $gvcf_para -O $gvcf_combined --tmp-dir $Tmp");
        system ("$Java  -jar -Xmx10g $GATK4_Loc GenotypeGVCFs -L $target_realign_bed -R $genome_file -V $gvcf_combined -O $vcf_dir/allSample.vcf.gz --sample-ploidy $vcf_ploidy --tmp-dir $Tmp");
        ### 添加硬过滤，根据GATK官网推荐https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
        system ("$Java  -jar -Xmx15g $GATK4_Loc VariantFiltration -R $genome_file -V $vcf_dir/allSample.vcf.gz -O $vcf_dir/allSample.filter.vcf.gz --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0' --filter-name 'my_filters'  --tmp-dir $Tmp");

        system("$bcftools index -t $vcf_dir/allSample.filter.vcf.gz"); #建索引
        system("$bcftools norm -m -both -f $genome_file -O z -o $vcf_dir/allSample.vcf.split.gz $vcf_dir/allSample.filter.vcf.gz"); #多态位点拆分
        
        print "\n";
    }
    else
    {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}