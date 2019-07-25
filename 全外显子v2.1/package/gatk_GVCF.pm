package package::gatk_GVCF;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start gatk_GVCF ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::gatk_GVCF($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in gatk_GVCF run para, please check\n" if($isOK == 0);

    my $report_dir = $hashConfig->{'Report'}; 
    my $output_dir = $hashConfig->{'Output'}; 
    my $gvcf_dir = "$report_dir/gvcf";
    package::utils::make_dir($gvcf_dir);

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control');
    my %hashCondition; 
    foreach my $sample(@samples)
    {
        my $gVCF      = "$gvcf_dir/$sample.g.vcf.gz";
        my $gVCF_idx  = "$gvcf_dir/$sample.g.vcf.gz.tbi";
        my $bam_final = "$output_dir/$sample/$sample\_final.bam";
        my $bai_final = "$output_dir/$sample/$sample\_final.bai";        
        # 已完成该步骤
        if(package::utils::is_file_ok($gVCF) and package::utils::is_file_ok($gVCF_idx) and not exists $hashConfig->{'Force_Run'}{$sample})
        {
            $hashCondition{"Finish"}{$sample} = "$gVCF $gVCF_idx";
            next;            
        }
        # 原始数据没问题
        if(package::utils::is_file_ok($bam_final) and package::utils::is_file_ok($bai_final)){
            $hashCondition{"Good2Run"}{$sample} = "$bam_final $bai_final";
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$bam_final $bai_final";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "gatk_GVCF", \%hashCondition); 
    die "Sample lost BAM :" . (join ",", sort keys %{$hashCondition{'Error'}}) . "\n" if(exists($hashCondition{'Error'}));
    
    # 执行 GVCF
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0)
    {
        my $threshold = exists $hashPara->{"Process"}{"GATK_GVCF"} ? $hashPara->{"Process"}{"GATK_GVCF"} : 10;
		   $threshold = $hashConfig->{"Process_GATK_GVCF"} if(exists $hashConfig->{"Process_GATK_GVCF"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::utils::process_bar_array($sample, \@sample_runs)});# 进度条
        foreach my $sample(@sample_runs)
        {
            $pm->start($sample) and next;
            run_gatk_gvcf($hashPara, $hashConfig, $sample);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else
    {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}

sub run_gatk_gvcf{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $sample     = shift @_;
    my $report_dir = $hashConfig->{'Report'}; 
    my $output_dir = "$hashConfig->{'Output'}/$sample"; 
    my $gvcf_dir   = "$report_dir/gvcf";
    my $Species    = $hashConfig->{"Species"};
    my $RealignBed = $hashConfig->{"RealignBed"}; # -L 能显著加快该步骤的运行速度,否则，gatk会根据基因组按照一定窗口进行遍历，很慢(主要原因在于，mapping中，realignment步骤保留了非bed区域的数据)

    my $tmp_dir    = $hashPara->{'Soft'}{'Tmp'};
    my $Java       = $hashPara->{"Soft"}{"Java"};
    my $GATK4_Loc  = $hashPara->{"Soft"}{"GATK4_Loc"};
    my $Genome     = $hashPara->{$Species}{"Genome"};
    my $exclude_96 = $hashPara->{$Species}{"exclude_96_bed"}; # 排除掉96位点所在区域

    my $XL_para    = package::utils::is_file_ok($exclude_96) ? "-XL $exclude_96" : "";
    my $vcf_ploidy = (exists($hashConfig->{"VCF_PLOIDY"}) and $hashConfig->{"VCF_PLOIDY"}=~/^\d+$/) ? $hashConfig->{"VCF_PLOIDY"} : 2;
    my $bam_final  = "$output_dir/$sample\_final.bam";
    my $gVCF       = "$gvcf_dir/$sample.g.vcf.gz";
    # 外显子13G,2.5小时,
    system ("$Java -jar -Xmx6g $GATK4_Loc HaplotypeCaller -I $bam_final -O $gVCF -R $Genome -ERC GVCF --min-base-quality-score 20  --sample-ploidy $vcf_ploidy --tmp-dir $tmp_dir -L $RealignBed $XL_para --max-reads-per-alignment-start 0 --dont-use-soft-clipped-bases true");
}
1