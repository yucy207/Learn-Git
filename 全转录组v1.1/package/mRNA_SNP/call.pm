package package::mRNA_SNP::call;
use strict;
use warnings;
$|=1;

sub run{

    my $metadata = shift;
    my $base     = shift;
    my @samples  = split /,/, $metadata->{'samples'};
    my $organ    = qq{$metadata->{organ}};

    my $GATK4_Loc  = qq{$base->{GATK4_Loc}};
    my $Java       = qq{$base->{Java}};
    my $bcftools   = qq{$base->{bcftools}};

    my $RealignBed = qq{$base->{$organ}{RealignBed}};
    my $Genome     = qq{$base->{$organ}{genome_fasta}};
    
    my $result     = qq{$metadata->{project}/mRNA/snp}; 
    my $report     = qq{$metadata->{report}/03_mRNA_Analysis/09_mRNA_SNP_Analysis};
    my $tmp_dir    = qq{/home/tmp};
    my $call_dir   = qq{$result/result/vcf};
      
    system qq{mkdir -p $report}         if not -d $report;
    system qq{mkdir -p $result/run}     if not -d qq{$result/run};
    system qq{mkdir -p $result/result}  if not -d qq{$result/result};
    system qq{mkdir -p $result/log}     if not -d qq{$result/log};
    system qq{mkdir -p $call_dir}       if not -d $call_dir;
   
    if (-e qq{$call_dir/call.finish}){

        print qq{call已经分析完成!\n};
        exit;
    }
   
    my @gVCFs       = map{"--variant $result/result/$_/$_.g.vcf.gz"} @samples;
    my $gVCF_string = join " ", @gVCFs;

    # 1、执行VCF
    my $gvcf_combined = qq{$call_dir/1.sample.all.g.vcf.gz};
    my $vcf_raw       = qq{$call_dir/2.sample.all.raw.vcf.gz};
    my $VCF_cmd1 = qq{$Java -jar $GATK4_Loc CombineGVCFs  -L $RealignBed -R $Genome $gVCF_string -O $gvcf_combined --tmp-dir $tmp_dir};
    my $VCF_cmd2 = qq{$Java -jar $GATK4_Loc GenotypeGVCFs -L $RealignBed -R $Genome -V $gvcf_combined -O $vcf_raw --sample-ploidy 2 --tmp-dir $tmp_dir};

    # 2、HARD硬过滤
        # 过滤hard filter ，根据GATK官网推荐https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
        # 虽然叫过滤，实际上只是在FILTER列上添加标签
    my $filter_vcf = ""; # 过滤后的vcf
    ## SNP 
    my $vcf_raw_snp        = qq{$call_dir/3.sample.all.raw.snp.vcf.gz};
    my $vcf_raw_snp_filter = qq{$call_dir/3.sample.all.raw.snp.filter.vcf.gz};
    my $SNP_cmd1 = qq{$Java -jar -Xmx15g $GATK4_Loc SelectVariants    -R $Genome -V $vcf_raw     -O $vcf_raw_snp        --select-type-to-include SNP   --tmp-dir $tmp_dir};
    my $SNP_cmd2 = qq{$Java -jar -Xmx15g $GATK4_Loc VariantFiltration -R $Genome -V $vcf_raw_snp -O $vcf_raw_snp_filter --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0' --filter-name 'my_snp_filter'  --tmp-dir $tmp_dir};
   
    ## INDEL
    my $vcf_raw_indel        = qq{$call_dir/4.sample.all.raw.indel.vcf.gz};
    my $vcf_raw_indel_filter = qq{$call_dir/4.sample.all.raw.indel.filter.vcf.gz};
    my $INDEL_cmd1 = qq{$Java -jar -Xmx15g $GATK4_Loc SelectVariants    -R $Genome -V $vcf_raw       -O $vcf_raw_indel        --select-type-to-include INDEL --tmp-dir $tmp_dir};
    my $INDEL_cmd2 = qq{$Java -jar -Xmx15g $GATK4_Loc VariantFiltration -R $Genome -V $vcf_raw_indel -O $vcf_raw_indel_filter --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 3.0' --filter-name 'my_indel_filter'  --tmp-dir $tmp_dir};

    ## MIXED 采用indel条件   类似于1       9117371 .       A       AAAAG,G 375.57 这种类型的突变记录无法被SNP/INDEL过滤出来,属于MIXED
    my $vcf_raw_mixed        = qq{$call_dir/5.sample.all.raw.mixed.vcf.gz};
    my $vcf_raw_mixed_filter = qq{$call_dir/5.sample.all.raw.mixed.filter.vcf.gz};
    my $MIXED_cmd1 = qq{$Java -jar -Xmx15g $GATK4_Loc SelectVariants    -R $Genome -V $vcf_raw       -O $vcf_raw_mixed        --select-type-to-include MIXED --tmp-dir $tmp_dir};
    my $MIXED_cmd2 = qq{$Java -jar -Xmx15g $GATK4_Loc VariantFiltration -R $Genome -V $vcf_raw_mixed -O $vcf_raw_mixed_filter --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 3.0' --filter-name 'my_indel_mixed'  --tmp-dir $tmp_dir};

    ## merge,结果重新合并成一个vcf文件
    my $vcf_raw_snp_filter_para   = package::mRNA_SNP::utils::is_file_ok($vcf_raw_snp_filter)   ? "-I $vcf_raw_snp_filter"   : "";
    my $vcf_raw_indel_filter_para = package::mRNA_SNP::utils::is_file_ok($vcf_raw_indel_filter) ? "-I $vcf_raw_indel_filter" : "";
    my $vcf_raw_mixed_filter_para = package::mRNA_SNP::utils::is_file_ok($vcf_raw_mixed_filter) ? "-I $vcf_raw_mixed_filter" : "";
    my $vcf_merged = qq{$call_dir/6.sample.all.merged.vcf.gz};
    my $merged_cmd = qq{$Java -jar -Xmx15g $GATK4_Loc MergeVcfs $vcf_raw_snp_filter_para $vcf_raw_indel_filter_para $vcf_raw_mixed_filter_para -O $vcf_merged --TMP_DIR $tmp_dir};

    $filter_vcf = $vcf_merged;

    # 3、使用bcftools 对vcf的多态位点拆分，左校正
    my $VCF_FINAL     = qq{$call_dir/sample.all.final.vcf.gz};
    my $VCF_FINAL_IDX = qq{$call_dir/sample.all.final.vcf.gz.tbi}; 
    my $bcf_cmd1 = qq{$bcftools norm -m -both -f $Genome -O z -o $VCF_FINAL $filter_vcf};
    my $bcf_cmd2 = qq{$bcftools index --tbi $VCF_FINAL}; 

    # 4、finish     
    my $finish = qq{touch $call_dir/call.finish};

    open SAVE, qq{>$result/run/call.sh} or die "Can't open $result/run/call.sh\n";
    print SAVE qq{$VCF_cmd1\n$VCF_cmd2\n};
    print SAVE qq{$SNP_cmd1\n$SNP_cmd2\n};
    print SAVE qq{$INDEL_cmd1\n$INDEL_cmd2\n};
    print SAVE qq{$MIXED_cmd1\n$MIXED_cmd2\n};
    print SAVE qq{$merged_cmd\n};
    print SAVE qq{$bcf_cmd1\n$bcf_cmd2\n};
    print SAVE qq{$finish\n};
    close SAVE;

    system qq{bash $result/run/call.sh &> $result/log/call.log};

    print qq{call分析完成!\n};

}

1