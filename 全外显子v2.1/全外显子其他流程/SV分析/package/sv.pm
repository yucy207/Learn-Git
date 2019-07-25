package package::sv;
use strict;
use warnings;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start sv analysis ".package::utils::get_time()." ##########\n";
    my $report_dir  = $hashConfig->{'Report'}; 
    my $output_dir  = $hashConfig->{'Output'}; 
    my $SVdbtype    = (exists $hashConfig->{'GenomeBed'}) ? "SVdbWGS" : "SVdbWES"; #判断是WES 还是WGS
	my $svReportDir = "$report_dir/sv";
	my $tmp_dir     = "$svReportDir/TMP";
	mkdir $svReportDir if(not -e $svReportDir);
	mkdir $tmp_dir if(not -e $tmp_dir);

    # 数据状态检测
    my @samples  = package::utils::get_sample($hashConfig, 'case', 'control');
    my %hashCondition; 
    
    foreach my $sample (@samples) { 
        my ($FinalBam, $SplitBam, $DiscoBam) = ("", "", "");
        if ($SVdbtype eq "SVdbWGS")
        {
            $FinalBam = "$output_dir/$sample/$sample\_final.bam";
            $SplitBam = "$output_dir/$sample/$sample.splitters.bam";
            $DiscoBam = "$output_dir/$sample/$sample.discordants.bam";   
        }
        if ($SVdbtype eq "SVdbWES")
        {
            $FinalBam = "$output_dir/$sample/svbam/$sample\_sv.bam";
            $SplitBam = "$output_dir/$sample/svbam/$sample.splitters.bam";
            $DiscoBam = "$output_dir/$sample/svbam/$sample.discordants.bam";   
        }

		my $sampel_tmp_dir       = "$tmp_dir/$sample";
		my $sample_readdepth_txt = "$sampel_tmp_dir/$sample.sv.$sample\_final.bam.readdepth.txt";
		my $sample_readdepth_bed = "$sampel_tmp_dir/$sample.sv.$sample\_final.bam.readdepth.bed";
		my $sample_sv_vcf        = "$sampel_tmp_dir/$sample.sv.vcf";
		
        # WES 需要生成svbam
        prepare ($hashConfig, $hashPara, $sample) if ($SVdbtype eq "SVdbWES" and (package::utils::is_file_ok($SplitBam, $DiscoBam) == 0));
             
        # 已完成该步骤
        if(package::utils::is_file_ok($sample_readdepth_txt, $sample_readdepth_bed, $sample_sv_vcf)) {
            $hashCondition{"Finish"}{$sample} = "$sample_readdepth_txt $sample_readdepth_bed $sample_sv_vcf";
            next;            
        }
        if(package::utils::is_file_ok($sample_sv_vcf) and $SVdbtype eq "SVdbWES") {
            $hashCondition{"Finish"}{$sample} = "$sample_sv_vcf";
            next;            
        }
        
        # 原始数据没问题
        if(package::utils::is_file_ok($FinalBam, $SplitBam, $DiscoBam)){
            $hashCondition{"Good2Run"}{$sample} = "$FinalBam $SplitBam $DiscoBam";
            next;
        }
        
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$FinalBam $SplitBam $DiscoBam";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "SV", \%hashCondition); 
    
    # 执行SV
    my @sample_runs = sort keys %{$hashCondition{'Good2Run'}};
    if(@sample_runs > 0) {
        my $threshold = 5;
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::utils::process_bar_array($sample, \@sample_runs)});# 进度条
        foreach my $sample(@sample_runs) {
            $pm->start($sample) and next;
            sub_run($hashPara, $hashConfig, $sample, $SVdbtype);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}

sub sub_run{
	my $hashPara        = shift @_;
	my $hashConfig      = shift @_;
	my $sample          = shift @_;
    my $SVdbtype        = shift @_;
	my $speedseq        = $hashPara   -> {'Soft'}{'speedseq'};
	my $Species         = $hashConfig -> {'Species'};
	my $report_dir      = $hashConfig -> {'Report'};
	my $output_dir      = $hashConfig -> {'Output'};
	my $Genome          = $hashPara   -> {$Species}{'Genome'};
	my $excludeinfo     = (exists($hashPara->{$Species}{'Excludebed'})) ?  "-x $hashPara->{$Species}{'Excludebed'}" : "";

	my $svReportDir     = "$report_dir/sv";
	my $tmp_dir         = "$svReportDir/TMP";
	my $sampel_tmp_dir  = "$tmp_dir/$sample";
	mkdir $sampel_tmp_dir if(not -e $sampel_tmp_dir);
    my ($FinalBam, $SplitBam, $DiscoBam) = ("", "", "");
    if ($SVdbtype eq "SVdbWGS")
    {
        $FinalBam = "$output_dir/$sample/$sample\_final.bam";
        $SplitBam = "$output_dir/$sample/$sample.splitters.bam";
        $DiscoBam = "$output_dir/$sample/$sample.discordants.bam";   
    }
    if ($SVdbtype eq "SVdbWES")
    {
        $FinalBam = "$output_dir/$sample/svbam/$sample\_sv.bam";
        $SplitBam = "$output_dir/$sample/svbam/$sample.splitters.bam";
        $DiscoBam = "$output_dir/$sample/svbam/$sample.discordants.bam";   
    }
	
	chdir $sampel_tmp_dir;# 路径跳转至输出目录，因为speedseq会在当前目录下生成一些中间结果
    
    # 外显子或全基因组低深度时，CNV分析会出问题，此时最好把-d参数删掉，深度过高时也要去掉,-k 保留临时目录，因为有时候speedseq会出错，无法生成.gz文件	
    system("$speedseq/speedseq sv -P -t 20 -o  $sampel_tmp_dir/$sample  -B $FinalBam -S $SplitBam -D $DiscoBam -R $Genome  -T  $sampel_tmp_dir -k $excludeinfo") if ($SVdbtype eq "SVdbWES");	
    system("$speedseq/speedseq sv -P -t 20 -d -o  $sampel_tmp_dir/$sample  -B $FinalBam -S $SplitBam -D $DiscoBam -R $Genome  -T  $sampel_tmp_dir -k $excludeinfo") if ($SVdbtype eq "SVdbWGS");
}

sub prepare{
	my $hashConfig = shift @_;
	my $hashPara   = shift @_;
	my $sample     = shift @_;
    
	my $outputdir  = $hashConfig -> {'Output'};
	my $samtools   = $hashPara   -> {'Soft'}{'SamTools'};
	my $samblaster = $hashPara   -> {'Soft'}{'samblaster'};
	my $sambamba   = $hashPara   -> {'Soft'}{'sambamba'};
	my $rawbam     = "$outputdir/$sample/$sample.bam";
	my $dir        = "$outputdir/$sample/svbam";
	my $tmpdir     = "$dir/temp";
	mkdir $dir if(not -e $dir);
	if(-e "$dir/$sample"."_sv.bam"){ return ;}
	mkdir $tmpdir if(not -e $tmpdir);	
	system("$samtools view -h $rawbam|$samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 --splitterFile $dir/spl_pipe --discordantFile $dir/disc_pipe|$sambamba view -S -f bam -l 0 /dev/stdin|$sambamba sort -t 4 --tmpdir=$tmpdir/full -o $dir/$sample"."_sv.bam /dev/stdin");##将原始bam拆分成三个文件：合格的reads.bam、split.sam、disc.sam
	system("gawk '{ if (\$0~\"^@\") { print; next } else { \$10=\"*\"; \$11=\"*\"; print } }' OFS=\"\\t\" $dir/spl_pipe |$sambamba view -S -f bam -l 0 /dev/stdin |$sambamba sort -t 4 -m 1G --tmpdir=$tmpdir/spl -o $dir/$sample.splitters.bam /dev/stdin");##将split.sam、disc.sam分别转化成对应的bam文件，且去除头部信息和冗余的序列、质量信息
	system("gawk '{ if (\$0~\"^@\") { print; next } else { \$10=\"*\"; \$11=\"*\"; print } }' OFS=\"\\t\" $dir/disc_pipe |$sambamba view -S -f bam /dev/stdin |$sambamba sort -t 4 -m 1G --tmpdir=$tmpdir/disc -o $dir/$sample.discordants.bam /dev/stdin");
	system("rm -rf $tmpdir");    
}

1