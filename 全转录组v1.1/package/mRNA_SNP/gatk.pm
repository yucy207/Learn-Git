package package::mRNA_SNP::gatk;
use Parallel::ForkManager;
use strict;
use warnings;
$|=1;

sub run
{

    my $metadata = shift;
    my $base     = shift;
    my @samples  = split /,/, $metadata->{'samples'};
    my $organ    = qq{$metadata->{organ}};

    my $sambamba  = qq{$base->{sambamba}};
    my $cpulimit  = qq{$base->{cpulimit}};
    my $Java      = qq{$base->{Java}};
    my $Picard    = qq{$base->{Picard}};
    my $GATK3     = qq{$base->{GATK3}};
    my $GATK4_Loc = qq{$base->{GATK4_Loc}};

    my $RealignBed     = qq{$base->{$organ}{RealignBed}};
    my $TargetBed      = qq{$base->{$organ}{TargetBed}};
    my $Genome         = qq{$base->{$organ}{genome_fasta}};
    my $knownDBsnp     = (exists $base->{$organ}{'DBsnp'})     ? "-known $base->{$organ}{'DBsnp'}"     :  "";
    my $knownInDel     = (exists $base->{$organ}{'InDel'})     ? "-known $base->{$organ}{'InDel'}"     :  "";
    my $knownInDelGold = (exists $base->{$organ}{'InDelGold'}) ? "-known $base->{$organ}{'InDelGold'}" :  "";
    my $knownSitesDBsnp     = (exists $base->{$organ}{'DBsnp'})     ? "--known-sites $base->{$organ}{'DBsnp'}"     :  "";
    my $knownSitesInDel     = (exists $base->{$organ}{'InDel'})     ? "--known-sites $base->{$organ}{'InDel'}"     :  "";
    my $knownSitesInDelGold = (exists $base->{$organ}{'InDelGold'}) ? "--known-sites $base->{$organ}{'InDelGold'}" :  "";
    
    my $map      = qq{$metadata->{project}/mapping/result};
    my $result   = qq{$metadata->{project}/mRNA/snp}; 
    my $report   = qq{$metadata->{report}/03_mRNA_Analysis/09_mRNA_SNP_Analysis};

    system qq{mkdir -p $report}         if not -d qq{$report};
    system qq{mkdir -p $result/run}     if not -d qq{$result/run};
    system qq{mkdir -p $result/result}  if not -d qq{$result/result};
    system qq{mkdir -p $result/log}     if not -d qq{$result/log};

    my @res_samples = res_check(qq{$result/result}, \@samples);

    if (scalar @res_samples == 0) {
        print qq{Final bam已经生成!\n};
        return 0;
    }

    my $max_threads = 10;
    my $pm = Parallel::ForkManager->new($max_threads);

    foreach my $x (@res_samples) {

        my $pid = $pm->start and next;

        my $tmp_dir  = qq{/home/tmp};
        my $gatk_dir = qq{$result/result/$x};
        system qq{mkdir -p $gatk_dir} if not -d $gatk_dir;
        
        my $bam = qq{$map/$x/accepted_hits.bam};
        
        # 排序
        my $bam_sort = qq{$map/$x/$x\_sort.bam};
        my $sort_cmd = qq{$sambamba sort $bam -m 6G -t 4 -o $bam_sort --tmpdir $tmp_dir};

        # 质量统计
        my $cover_metrics   = qq{$gatk_dir/$x.metrics}; 
        my $cover_tag       = qq{$gatk_dir/$x.tag};
        my $quality_metrics = qq{$gatk_dir/$x.quality.metrics};
        my $statistics_cmd1 = qq{$cpulimit -l 1000 $Java -jar -Xmx3g $Picard CollectHsMetrics            INPUT=$bam_sort OUTPUT=$cover_metrics R=$Genome BI=$TargetBed TI=$TargetBed PER_TARGET_COVERAGE=$cover_tag MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir};
        my $statistics_cmd2 = qq{$cpulimit -l 1000 $Java -jar -Xmx3g $Picard CollectQualityYieldMetrics  INPUT=$bam_sort OUTPUT=$quality_metrics VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir};
        
        # 去重复
        my $bam_sort_dup         = qq{$map/$x/$x\_sort_dup.bam};
        my $bam_sort_dup_metrics = qq{$gatk_dir/$x\_sort_dup.metrics};
        my $dup_cmd              = qq{$cpulimit -l 1000 $Java -jar -Xmx16g $Picard MarkDuplicates INPUT=$bam_sort OUTPUT=$bam_sort_dup METRICS_FILE=$bam_sort_dup_metrics CREATE_INDEX=true REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir};
        
        # Split'N'Trim
        my $bam_trim     = qq{$map/$x/$x\_sort_dup_trim.bam};
        my $bam_trim_idx = qq{$map/$x/$x\_sort_dup_trim.bai};
        my $trim_cmd1 = qq{$cpulimit -l 1000 $Java -jar $GATK3 -T SplitNCigarReads -R $Genome -I $bam_sort_dup -o $bam_trim -U ALLOW_N_CIGAR_READS -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60};
        my $trim_cmd2 = qq{$cpulimit -l 1000 $Java -jar -Xmx8g $Picard BuildBamIndex INPUT=$bam_trim OUTPUT=$bam_trim_idx VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir} if(not -e $bam_trim_idx);

        # 根据数据库，对已知插入缺失位点校正，同时，只保留目标区域的数据，减少数据量(如果后面分析突变时，只使用GATK的HaplotypeCaller方法，这一步就没有必要做)(实际上，非目标区域的reads也保留了一部分，原因未知)
        my $bam_realign           = qq{$map/$x/$x\_sort_dup_trim_realign.bam};
        my $bam_realign_intervals = qq{$map/$x/$x\_sort_dup_trim_realign.intervals};
        my $realign_cmd1 = qq{$cpulimit -l 1000 $Java -Djava.io.tmpdir=$tmp_dir -jar -Xmx5g $GATK3 -l INFO -T RealignerTargetCreator -R $Genome -I $bam_trim -L $RealignBed -o $bam_realign_intervals $knownInDel $knownInDelGold --validation_strictness LENIENT};
        my $realign_cmd2 = qq{$cpulimit -l 1000 $Java -Djava.io.tmpdir=$tmp_dir -jar -Xmx5g $GATK3 -l INFO -T IndelRealigner         -R $Genome -I $bam_trim -L $RealignBed -o $bam_realign           $knownInDel $knownInDelGold  -targetIntervals $bam_realign_intervals --validation_strictness LENIENT};

        # 碱基质量校正
        my $bam_realign_recalibrator_table = qq{$map/$x/$x\_sort_dup_trim_realign_recalibrator.table};
        my $bam_realign_recalibrator       = qq{$map/$x/$x\_sort_dup_trim_realign_recalibrator.bam};
        my $bam_realign_recalibrator_idx   = qq{$map/$x/$x\_sort_dup_trim_realign_recalibrator.bai};
        my $BQSR_cmd1 = qq{$cpulimit -l 1000 $Java -jar -Xmx5g $GATK4_Loc BaseRecalibrator -R $Genome -I $bam_realign  -O $bam_realign_recalibrator_table $knownSitesDBsnp $knownSitesInDelGold $knownSitesInDel --tmp-dir $tmp_dir}; # 缺失VCF时，无法进行该步骤
        my $BQSR_cmd2 = qq{$cpulimit -l 1000 $Java -jar -Xmx8g $GATK4_Loc ApplyBQSR        -R $Genome -I $bam_realign  -O $bam_realign_recalibrator --bqsr-recal-file $bam_realign_recalibrator_table --tmp-dir $tmp_dir};

        # 最终比对bam
        my $bam_final_input     = ($knownSitesDBsnp=~/\w/) ? $bam_realign_recalibrator : $bam_realign;
        my $bam_final_input_idx = $bam_final_input;
        $bam_final_input_idx =~s/\.bam$/\.bai/;
        
        # 重命名
        my $bam_final     = qq{$map/$x/$x\_final.bam};
        my $bam_final_idx = qq{$map/$x/$x\_final.bai};
        my $mv1 = qq{mv $bam_final_input     $bam_final};
        my $mv2 = qq{mv $bam_final_input_idx $bam_final_idx};

        # finish文件
        my $finish  = qq{touch $gatk_dir/$x.gatk.finish};

        open SAVE, qq{>$result/run/$x.gatk.sh} or die "Can't open $result/run/$x.gatk.sh\n";
        print SAVE qq{$sort_cmd\n};
        print SAVE qq{$statistics_cmd1\n$statistics_cmd2\n};
        print SAVE qq{$dup_cmd\n};
        print SAVE qq{$trim_cmd1\n$trim_cmd2\n};
        print SAVE qq{$realign_cmd1\n$realign_cmd2\n};
        if($knownSitesDBsnp=~/\w/){
            print SAVE qq{$BQSR_cmd1\n$BQSR_cmd2\n};
        }
        print SAVE qq{$mv1\n$mv2\n};
        print SAVE qq{$finish\n};
        close SAVE;

        system qq{bash $result/run/$x.gatk.sh &> $result/log/$x.gatk.log};
        $pm->finish;

    }
 
    $pm->wait_all_children;

    print qq{Final bam生成!\n};

}

sub res_check
{
    my $output = shift;
    my $sample = shift;

    my @temp_sample = ();
    foreach my $x (@{$sample}) {
        next if -e qq{$output/$x/$x.gatk.finish};
        push @temp_sample, $x;
    }

    return @temp_sample;
}

1