package package::gatk_VCF;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start gatk_VCF ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::gatk_VCF($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in gatk_VCF run para, please check\n" if($isOK == 0);

    my $report_dir = $hashConfig->{'Report'}; 
    my $vcf_dir = "$report_dir/vcf";
    package::utils::make_dir($vcf_dir);
    #####
    # 数据状态检测
    #####
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control');
    my %hashCondition; 
    foreach my $sample(@samples){
        my $gVCF     = "$report_dir/gvcf/$sample.g.vcf.gz";
        my $gVCF_idx = "$report_dir/gvcf/$sample.g.vcf.gz.tbi";      
        # 原始数据没问题
        if(package::utils::is_file_ok($gVCF) and package::utils::is_file_ok($gVCF_idx)){
            $hashCondition{"Good2Run"}{$sample} = $gVCF;
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = $gVCF;
    }
    #####
    # 写日志
    #####
    package::utils::write_log("$report_dir/run.log", "gatk_VCF", \%hashCondition);
    die "Sample lost gvcf :" . (join ",", keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));
    
    #####
    # 执行VCF
    #####
    my $VCF_FINAL     = "$vcf_dir/sample.all.final.vcf.gz";
    my $VCF_FINAL_IDX = "$vcf_dir/sample.all.final.vcf.gz.tbi"; 
    if(package::utils::is_file_ok($VCF_FINAL) and package::utils::is_file_ok($VCF_FINAL_IDX) and package::utils::is_force_sample($hashConfig) == 0)
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }
    my $Species    = $hashConfig->{"Species"};
    my $RealignBed = $hashConfig->{"RealignBed"};
    my $vcf_ploidy = (exists($hashConfig->{"VCF_PLOIDY"}) and $hashConfig->{"VCF_PLOIDY"}=~/^\d+$/)         ? $hashConfig->{"VCF_PLOIDY"} : 2;
    my $vcf_split  = (exists($hashConfig->{"VCF_SPLIT"})  and $hashConfig->{"VCF_SPLIT"}=~/^TRUE$|^FALSE$/) ? $hashConfig->{"VCF_SPLIT"}  : 'TRUE';
    my $correct    = (exists($hashConfig->{"Correct"})  and $hashConfig->{"Correct"}=~/^VQSR$|^HARD$/)      ? $hashConfig->{"Correct"}    : 'HARD'; # 校正方式

    my $tmp_dir   = $hashPara->{'Soft'}{'Tmp'};
    my $Java      = $hashPara->{"Soft"}{"Java"};
    my $GATK4_Loc = $hashPara->{"Soft"}{"GATK4_Loc"};
    my $bcftools  = $hashPara->{"Soft"}{"bcftools"};
    my $Genome    = $hashPara->{$Species}{"Genome"};       

    my @gVCFs       = map{ "--variant $hashCondition{'Good2Run'}{$_}"} @samples;
    my $gVCF_string = join " ", @gVCFs;

    #####
    # GVCF 分型，由于分析较慢，故按照染色体拆分处理
    #####
    my %hashBED       = read_bed_chr($RealignBed);
    my @chrs          = sort {$hashBED{$a}{'Chr_Sort_Order'} <=> $hashBED{$b}{'Chr_Sort_Order'}} keys %hashBED;
    my $threshold     = exists $hashPara->{"Process"}{"GATK_VCF"} ? $hashPara->{"Process"}{"GATK_VCF"} : 10;
	   $threshold     = $hashConfig->{"Process_GATK_VCF"} if(exists $hashConfig->{"Process_GATK_VCF"});
    my $chr_split_dir = "$vcf_dir/chr_split";
    package::utils::make_dir($chr_split_dir);

    my @split_vcfs; # vcf列表
    my $pm = Parallel::ForkManager->new($threshold);
       $pm->run_on_start(sub{my ($pid, $chr) = @_; package::utils::process_bar_array($chr, \@chrs)});# 进度条
    foreach my $chr(@chrs)
    {
        my $chr_split_bed        = "$chr_split_dir/CHR.$chr.bed"; # bed拆分文件
        my $gvcf_combined_split  = "$chr_split_dir/CHR.$chr.1.sample.all.g.vcf.gz";
        my $vcf_raw_split        = "$chr_split_dir/CHR.$chr.2.sample.all.raw.vcf.gz";
        push @split_vcfs, $vcf_raw_split;

        $pm->start($chr) and next;

        output_chr_split_bed($hashBED{$chr}{'Chr_Bed_Region'}, $chr_split_bed);
        system ("$Java  -jar  $GATK4_Loc CombineGVCFs  -L $chr_split_bed -R $Genome $gVCF_string        -O $gvcf_combined_split                       --tmp-dir $tmp_dir");
        system ("$Java  -jar  $GATK4_Loc GenotypeGVCFs -L $chr_split_bed -R $Genome -V $gvcf_combined_split   -O $vcf_raw_split --sample-ploidy $vcf_ploidy --tmp-dir $tmp_dir");
        
        $pm->finish;    
    }
    $pm->wait_all_children;

    #####
    # VCF 合并
    #####    
    my $vcf_raw           = "$vcf_dir/2.sample.all.raw.vcf.gz";
    my $split_vcf_list    = "";
    map{ $split_vcf_list .= "-I $_ " }@split_vcfs;
    system ("$Java  -jar -Xmx15g $GATK4_Loc MergeVcfs $split_vcf_list -O $vcf_raw --TMP_DIR $tmp_dir");

    #####
    # VCF 后续过滤、预处理
    #####
    my $filter_vcf = ""; # 过滤后的vcf
    # 1 VQSR 过滤
    if($correct eq 'VQSR')
    {
        my $Rbin      = $hashPara->{"Soft"}{"R"};
        my $RLib      = $hashPara->{"Soft"}{"RLib"};
        my $RDir      = package::utils::get_dirname($Rbin);
        $ENV{'PATH'}   = "$RDir:$ENV{'PATH'}"; # 添加R环境，VQSR 要用到
        $ENV{'R_LIBS'} = $RLib;
        my $HAPMAP       = (exists $hashPara->{$Species}{"VQSR_HAPMAP"})       ? "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hashPara->{$Species}{'VQSR_HAPMAP'}"       :  "";
        my $OMNI_1000G   = (exists $hashPara->{$Species}{"VQSR_1000G_OMNI"})   ? "--resource:omni,known=false,training=true,truth=false,prior=12.0 $hashPara->{$Species}{'VQSR_1000G_OMNI'}"    :  "";
        my $SNP_HC_1000G = (exists $hashPara->{$Species}{"VQSR_1000G_SNP_HC"}) ? "--resource:1000G,known=false,training=true,truth=false,prior=10.0 $hashPara->{$Species}{'VQSR_1000G_SNP_HC'}" :  "";
        my $DBSNP        = (exists $hashPara->{$Species}{"DBsnp"})             ? "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $hashPara->{$Species}{'DBsnp'}"              :  "";
        my $InDelGold    = (exists $hashPara->{$Species}{"InDelGold"})         ? "--resource:mills,known=false,training=true,truth=true,prior=12.0 $hashPara->{$Species}{'InDelGold'}"          :  "";

        die "NO DATABASE IN VQSR MODEL\n" if($HAPMAP eq "" or $OMNI_1000G eq "" or $SNP_HC_1000G eq "" or $DBSNP eq "" or $InDelGold eq "");
        
        # SNP VQSR    
        my $VQSR_SNP_Recal    = "$vcf_dir/3.sample.all.raw.VQSR.snp.recal";
        my $VQSR_SNP_tranches = "$vcf_dir/3.sample.all.raw.VQSR.snp.tranches";
        my $VQSR_SNP_plotR    = "$vcf_dir/3.sample.all.raw.VQSR.snp.plots.R";
        my $vcf_snp_vqsr      = "$vcf_dir/3.sample.all.raw.VQSR.snp.vcf.gz";
        system("$Java -jar -Xmx15g $GATK4_Loc VariantRecalibrator -R $Genome -V $vcf_raw -mode SNP $HAPMAP $OMNI_1000G $SNP_HC_1000G $DBSNP  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --output $VQSR_SNP_Recal --tranches-file $VQSR_SNP_tranches --rscript-file $VQSR_SNP_plotR");
        system("$Java -jar -Xmx15g $GATK4_Loc ApplyVQSR           -R $Genome -V $vcf_raw -mode SNP -O $vcf_snp_vqsr --truth-sensitivity-filter-level 99.0 --tranches-file $VQSR_SNP_tranches --recal-file $VQSR_SNP_Recal");

        # INDEL VQSR    
        my $VQSR_INDEL_Recal    = "$vcf_dir/4.sample.all.raw.VQSR.indel.recal";
        my $VQSR_INDEL_tranches = "$vcf_dir/4.sample.all.raw.VQSR.indel.tranches";
        my $VQSR_INDEL_plotR    = "$vcf_dir/4.sample.all.raw.VQSR.indel.plots.R";
        my $vcf_snp_indel_vqsr      = "$vcf_dir/4.sample.all.raw.VQSR.snp.indel.vcf.gz";
        system("$Java -jar -Xmx15g $GATK4_Loc VariantRecalibrator -R $Genome -V $vcf_snp_vqsr -mode INDEL $InDelGold $DBSNP  -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum --max-gaussians 4  --output $VQSR_INDEL_Recal --tranches-file $VQSR_INDEL_tranches --rscript-file $VQSR_INDEL_plotR");
        system("$Java -jar -Xmx15g $GATK4_Loc ApplyVQSR           -R $Genome -V $vcf_snp_vqsr -mode INDEL -O $vcf_snp_indel_vqsr --truth-sensitivity-filter-level 99.5 --tranches-file $VQSR_INDEL_tranches --recal-file $VQSR_INDEL_Recal");

        $filter_vcf = $vcf_snp_indel_vqsr;
    }
    else # HARD 硬过滤
    {
        # 过滤hard filter ，根据GATK官网推荐https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
        # 虽然叫过滤，实际上只是在FILTER列上添加标签
    
        # SNP （4个外显子，SNP,INDEL,MIXED 三步加起来2分钟左右）
        my $vcf_raw_snp        = "$vcf_dir/3.sample.all.raw.snp.vcf.gz";
        my $vcf_raw_snp_filter = "$vcf_dir/3.sample.all.raw.snp.filter.vcf.gz";
        system ("$Java  -jar -Xmx15g $GATK4_Loc SelectVariants    -R $Genome -V $vcf_raw     -O $vcf_raw_snp        --select-type-to-include SNP   --tmp-dir $tmp_dir");
        system ("$Java  -jar -Xmx15g $GATK4_Loc VariantFiltration -R $Genome -V $vcf_raw_snp -O $vcf_raw_snp_filter --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0' --filter-name 'my_snp_filter'  --tmp-dir $tmp_dir");
      
        # INDEL
        my $vcf_raw_indel        = "$vcf_dir/4.sample.all.raw.indel.vcf.gz";
        my $vcf_raw_indel_filter = "$vcf_dir/4.sample.all.raw.indel.filter.vcf.gz";
        system ("$Java  -jar -Xmx15g $GATK4_Loc SelectVariants    -R $Genome -V $vcf_raw       -O $vcf_raw_indel        --select-type-to-include INDEL --tmp-dir $tmp_dir");
        system ("$Java  -jar -Xmx15g $GATK4_Loc VariantFiltration -R $Genome -V $vcf_raw_indel -O $vcf_raw_indel_filter --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 3.0' --filter-name 'my_indel_filter'  --tmp-dir $tmp_dir");

        # MIXED 采用indel条件   类似于1       9117371 .       A       AAAAG,G 375.57 这种类型的突变记录无法被SNP/INDEL过滤出来,属于MIXED
        my $vcf_raw_mixed        = "$vcf_dir/5.sample.all.raw.mixed.vcf.gz";
        my $vcf_raw_mixed_filter = "$vcf_dir/5.sample.all.raw.mixed.filter.vcf.gz";
        system ("$Java  -jar -Xmx15g $GATK4_Loc SelectVariants    -R $Genome -V $vcf_raw       -O $vcf_raw_mixed        --select-type-to-include MIXED --tmp-dir $tmp_dir");
        system ("$Java  -jar -Xmx15g $GATK4_Loc VariantFiltration -R $Genome -V $vcf_raw_mixed -O $vcf_raw_mixed_filter --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 3.0' --filter-name 'my_indel_mixed'  --tmp-dir $tmp_dir");

        # merge,结果重新合并成一个vcf文件 （4个外显子，1分钟左右）
        my $vcf_raw_snp_filter_para   = package::utils::is_file_ok($vcf_raw_snp_filter)   ? "-I $vcf_raw_snp_filter"   : "";
        my $vcf_raw_indel_filter_para = package::utils::is_file_ok($vcf_raw_indel_filter) ? "-I $vcf_raw_indel_filter" : "";
        my $vcf_raw_mixed_filter_para = package::utils::is_file_ok($vcf_raw_mixed_filter) ? "-I $vcf_raw_mixed_filter" : "";
        my $vcf_merged = "$vcf_dir/6.sample.all.merged.vcf.gz";
        system ("$Java  -jar -Xmx15g $GATK4_Loc MergeVcfs $vcf_raw_snp_filter_para $vcf_raw_indel_filter_para $vcf_raw_mixed_filter_para -O $vcf_merged --TMP_DIR $tmp_dir");

        $filter_vcf = $vcf_merged;
    }


    # 使用bcftools 对vcf的多态位点拆分，左校正
    if($vcf_split eq 'TRUE')
    {
        system("$bcftools norm -m -both -f $Genome -O z -o $VCF_FINAL $filter_vcf");
        system("$bcftools index --tbi $VCF_FINAL");       
    }
    else
    {
        system("cp $filter_vcf $VCF_FINAL");
        system("cp $filter_vcf.tbi $VCF_FINAL_IDX");
    }
}

# 输出染色体bed区域
sub output_chr_split_bed{
    my $hashChrBED    = shift @_;
    my $chr_split_bed = shift @_; 
    open BED, ">$chr_split_bed";
    foreach my $count(sort {$a <=> $b} keys %$hashChrBED)
    {
        print BED "$hashChrBED->{$count}";
    }
    close BED;
}
# 读取bed文件染色体
sub read_bed_chr{
    my $bed_file = shift @_;
    my %hashBED;
    open BED, $bed_file;
    my $count = 0;
    while(<BED>)
    {   
        $count++;
        my ($chr, $tmp) = split /\t/, $_;
        $hashBED{$chr}{'Chr_Sort_Order'}         = $count;
        $hashBED{$chr}{'Chr_Bed_Region'}{$count} = $_;
    }
    close BED;
    return %hashBED;
}

1