package package::excavator2;
use strict;
use warnings;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start excavator2 ".package::main::get_time()." ##########\n";
    my $isOK = package::check::excavator2($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in excavator2 run para, please check\n" if($isOK == 0);

	my $hash_group     = $hashConfig->{'Group'};
	my $hash_gender    = $hashConfig->{'Sample'};
	my $report_dir     = $hashConfig->{'Report'}; 
    my $output_dir     = $hashConfig->{'Output'};
	my $species        = $hashConfig->{'Species'};
	my $region_orifile = $hashConfig->{'RealignBed'};
	my $bedtools       = $hashPara->{'Soft'}{'bedtools'};
	my $excavator_dir  = $hashPara->{'Soft'}{'EXCAVATOR2'};
	my $RScript        = $hashPara->{'Soft'}{'RScript'};
	my $genome         = $hashPara->{$species}{'Genome'};
	my $assembly       = $hashPara->{$species}{'AnnovarBuild'};
	my $processors     = (exists $hashPara->{'Process'}{'excavator2'}) ? $hashPara->{'Process'}{'excavator2'} : 6;
	   $processors     = $hashConfig->{"Process_excavator2"} if(exists $hashConfig->{"Process_excavator2"});
	my $window         = (exists $hashConfig->{'Windows'} and $hashConfig->{'Windows'} =~ /\d/) ? $hashConfig->{'Windows'} : 50000;  #设置窗口大小
    my @config_samples = package::main::get_all_sample($hashConfig, 'case', 'control');
	
	my $result_dir     = "$report_dir/excavator2";
	my $region_file    = "$result_dir/merge_".(split /\//, $region_orifile)[-1];
	my $target         = (split /\./, (split /\//, $region_file)[-1])[0]."_w$window";
    package::main::make_dir($result_dir);
	
	# # bed合并重叠区域
	system "$bedtools sort -i $region_orifile > $region_file\_sort"; # # 合并重叠区域
    system "$bedtools merge -i $region_file\_sort -d 1 > $region_file";
	
    # # excavator2运行
    TargetPerla($result_dir, $window, $region_file, $assembly, $genome, $excavator_dir, $target);	
	EXCAVATORDataPrepare($report_dir, $result_dir, $window, $assembly, $excavator_dir, $target, $output_dir, $processors, \@config_samples);		
	EXCAVATORDataAnalysis($hash_group, $report_dir, $result_dir, $window, $assembly, $excavator_dir, $target, $processors);
	
	# # 性别校正
	gender_correct($hash_group, $hash_gender, $report_dir, $result_dir);
	
	# # CNV区域绘图
	excavator2_region_plot($hashConfig, $hashPara, $hash_group, $RScript, $report_dir, $result_dir);
}

sub excavator2_region_plot {
    my $hashConfig       = shift @_;
	my $hashPara         = shift @_;
	my $hash_group       = shift @_;
	my $RScript          = shift @_;
	my $report_dir       = shift @_;
	my $result_dir       = shift @_;
	my $analysis_dir     = "$result_dir/analysis";
	my $bed_for_plot_dir = "$result_dir/bed_for_plot";
	my $document_dir     = "$report_dir/document";
	my $cnv_dir          = "$document_dir/3_CNV";
	my $all_plot_dir     = "$cnv_dir/excavator2_plot";
	my $out_dir          = "$all_plot_dir/EXCAVATOR2_Original_Plot";
	package::main::make_dir($bed_for_plot_dir);
	package::main::make_dir($document_dir);
	package::main::make_dir($cnv_dir);
	package::main::make_dir($all_plot_dir);
	package::main::make_dir($out_dir);
    print "\nProcessing excavator2_region_plot ...\n";
	my %hashCondition = check_plot_file($hash_group, $report_dir, $analysis_dir);
	
	my @case_suffixes = keys %{$hashCondition{"Good2Run"}};
	if(@case_suffixes > 0) { # # 运行
        my $threshold = exists $hashPara->{"Process"}{"excavator2_plot"} ? $hashPara->{"Process"}{"excavator2_plot"} : 10;
		   $threshold = $hashConfig->{"Process_excavator2_plot"} if(exists $hashConfig->{"Process_excavator2_plot"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $case_suffix) = @_; package::main::process_bar_array($case_suffix, \@case_suffixes)});# 进度条
        foreach my $case_suffix(@case_suffixes) {
            $pm->start($case_suffix) and next;
			sub_excavator2_region_plot($hash_group, $RScript, $report_dir, $result_dir, $case_suffix);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}

sub sub_excavator2_region_plot {
    my $hash_group       = shift @_;
	my $RScript          = shift @_;
	my $report_dir       = shift @_;
	my $result_dir       = shift @_;
	my $case_suffix      = shift @_;
	my $analysis_dir     = "$result_dir/analysis";
	my $bed_for_plot_dir = "$result_dir/bed_for_plot";
	my $document_dir     = "$report_dir/document";
	my $cnv_dir          = "$document_dir/3_CNV";
	my $all_plot_dir     = "$cnv_dir/excavator2_plot";
	my ($case, $suffix_tmp) = split /\.VS\./, $case_suffix;
	my $suffix           = "VS.$suffix_tmp";
	
	my @controls         = package::main::get_sample($hash_group->{$suffix}, 'control');
	print "$case_suffix ...";
	my $case_suffix_dir      = "$bed_for_plot_dir/$case_suffix";
	my $out_dir              = "$all_plot_dir/EXCAVATOR2_Original_Plot";
	my $out_case_suffix_dir  = "$out_dir/$case_suffix";
	my $excavator2_plot_file = "$case_suffix_dir/excavator2_PlotResultsRegion.pdf";
    package::main::make_dir($case_suffix_dir);					
    package::main::make_dir($out_case_suffix_dir);			
	my $analysis_hslm        = "$analysis_dir/$suffix/Results/$case/HSLMResults_$case.txt";
	my $analysis_plots       = "$analysis_dir/$suffix/Plots/$case/*";
	my $analysis_file        = "$analysis_dir/$suffix/Results/$case/FastCallResults_$case.txt";
	system "$RScript ./package/excavator2_plot.R $analysis_hslm $analysis_file $case_suffix_dir";
	system "cp $excavator2_plot_file $out_case_suffix_dir"; # # 移动至document文件夹
	system "cp $analysis_plots $out_case_suffix_dir";       # # 移动至document文件夹
	print "OK\n";
}

sub gender_correct {
    my $hash_group   = shift @_;
	my $hash_gender  = shift @_;
	my $report_dir   = shift @_;
	my $result_dir   = shift @_;
	my $analysis_dir = "$result_dir/analysis";
	my $correct_dir  = "$result_dir/result_final";
	package::main::make_dir($correct_dir);
    print "\nProcessing gender_correct ...\n";
	my %hashCondition = check_gender_correct_file($hash_group, $report_dir, $correct_dir, $analysis_dir);
	
	foreach my $suffix(keys %hashCondition) {
	    my @controls        = package::main::get_sample($hash_group->{$suffix}, 'control');
		my $controls_gender = get_control_gender($hash_gender, \@controls);
		
		foreach my $case(keys %{$hashCondition{$suffix}{"Good2Run"}}) {
		    my $case_suffix_dir  = "$correct_dir/$case.$suffix";	
            package::main::make_dir($case_suffix_dir);			
			my $analysis_file    = $hashCondition{$suffix}{"Good2Run"}{$case};
			my $analysis_correct = "$case_suffix_dir/$case.$suffix\_result.txt";		
			print "Reading $analysis_file ...";
			open CORRECT, ">$analysis_correct";
			open ANA, $analysis_file;
			<ANA>;
			while(<ANA>) {
			    $_ =~ s/[\r\n]//g;
				my ($chr, $start, $end, $log2, $cnf, $cnv, $call, $call_prob) = split /\t/, $_;
				# # 扩增/缺失 类型判定
				my $type = "normal";
				if($chr =~ /[Xx]/) { # # X染色体拷贝数校正
				    my $benchmark;
					($benchmark, $cnf, $cnv) = sub_gender_correct($hash_gender, $case, $controls_gender, $cnf, $cnv);
					next if($cnv == $benchmark);
					$type = ($cnv > $benchmark) ? "gain" : "loss";
				}
				else {
				    next if($call == 0);
					$type = ($call > 0) ? "gain" : "loss";
				}
				# # 数据存储
				print CORRECT "$chr\t$start\t$end\t$cnv\t$type\t$log2\t$call_prob\n";
			}
			close CORRECT;
			close ANA;
			print "OK\n";
		}
	}
}

sub sub_gender_correct {
    my $hash_gender     = shift @_;
	my $case            = shift @_;
	my $controls_gender = shift @_;
	my $cnf             = shift @_;
	my $cnv             = shift @_;
	my $case_gender     = $hash_gender->{$case}{'sex'};
	my $benchmark = 2; # # 判断插入缺失的标准
	if($case_gender eq 'male' and $controls_gender eq 'female') { 
	    $benchmark = 1;
	}
	elsif($case_gender eq 'female' and $controls_gender eq 'male') {
	    $cnf /= 2;
		$cnv  = sprintf "%.f", $cnf; # # 四舍五入
	}
	elsif($case_gender eq 'male' and $controls_gender eq 'male') {
	    $cnf /= 2;
		$cnv  = sprintf "%.f", $cnf; # # 四舍五入
		$benchmark = 1;
	}
	return ($benchmark, $cnf, $cnv);
}

sub get_control_gender {
	my $hash_gender = shift @_;
	my $controls    = shift @_;
	my ($male_num, $female_num) = (0, 0);
	foreach my $control(@{$controls}){ # # control性别以多数为准，相等设为女性
		$male_num ++ if($hash_gender->{$control}{'sex'} eq 'male');
		$female_num ++ if($hash_gender->{$control}{'sex'} eq 'female');
	}
	my $controls_gender = $male_num > $female_num ? 'male' : 'female';
}

sub check_plot_file {
    my $hash_group   = shift @_;
	my $report_dir   = shift @_;
    my $analysis_dir = shift @_;
	my %hashCondition; 
    foreach my $suffix(keys %{$hash_group}) {	
		my @cases    = package::main::get_sample($hash_group->{$suffix}, 'case');
		foreach my $case(@cases) {
		    my $analysis_file    = "$analysis_dir/$suffix/Results/$case/FastCallResults_$case.txt";
			# # 原始数据没问题
            if(package::main::is_file_ok($analysis_file)) {
                $hashCondition{"Good2Run"}{"$case.$suffix"} = "$analysis_file";
                next;
            }
            # # 原始数据丢失
            $hashCondition{"Error"}{"$case.$suffix"} = "$analysis_file";
		}
		# # 写日志
        package::main::write_log("$report_dir/run.log", "excavator2_plot", \%hashCondition); 
    }
    return %hashCondition;
}

sub check_gender_correct_file {
    my $hash_group   = shift @_;
	my $report_dir   = shift @_;
	my $correct_dir  = shift @_;
    my $analysis_dir = shift @_;
	my %hashCondition; 
    foreach my $suffix(keys %{$hash_group}) {	
		my @cases    = package::main::get_sample($hash_group->{$suffix}, 'case');
		foreach my $case(@cases) {
		    my $case_suffix_dir  = "$correct_dir/$case.$suffix";
		    my $analysis_file    = "$analysis_dir/$suffix/Results/$case/FastCallResults_$case.txt";
			my $analysis_correct = "$case_suffix_dir/$case.$suffix\_result.txt";
            # 已完成该步骤
            if(package::main::is_file_ok($analysis_correct)) {
                $hashCondition{$suffix}{"Finish"}{$case} = "$analysis_correct";
                next;            
            }
			# # 原始数据没问题
            if(package::main::is_file_ok($analysis_file)) {
                $hashCondition{$suffix}{"Good2Run"}{$case} = "$analysis_file";
                next;
            }
            # # 原始数据丢失
            $hashCondition{$suffix}{"Error"}{$case} = "$analysis_file";
		}
		# # 写日志
        package::main::write_log("$report_dir/run.log", "gender_correct", $hashCondition{$suffix}); 
    }
    return %hashCondition;
}

sub EXCAVATORDataAnalysis {
    my $hash_group     = shift @_;
	my $report_dir     = shift @_;
	my $result_dir     = shift @_;
	my $window         = shift @_;
	my $assembly       = shift @_;
	my $excavator_dir  = shift @_;
	my $target         = shift @_;
	my $processors     = shift @_;
	my $mode           = "pooling"; # 设置模式：pooling/paired
	print "\nProcessing EXCAVATORDataAnalysis.pl ...\n";
	my $input_dir      = "$result_dir/input";
	my $prepare_dir    = "$result_dir/prepare";
	my $analysis_dir   = "$result_dir/analysis";
	package::main::make_dir($input_dir);
	package::main::make_dir($prepare_dir);
	package::main::make_dir($analysis_dir);
	my %hashCondition = check_analysis_file($hash_group, $report_dir, $analysis_dir, $prepare_dir);
	foreach my $suffix(sort {$a cmp $b} keys %hashCondition) {
		my @cases        = keys %{$hashCondition{$suffix}{"Good2Run"}};
		next if(!@cases);
        my @controls     = package::main::get_sample($hash_group->{$suffix}, 'control');
        my $group_dir    = "$analysis_dir/$suffix";
		package::main::make_dir($group_dir);
		my $fileAnalysis = "$input_dir/ExperimentalFileAnalysis_$suffix.w$window.txt";
		write_file_analysis($fileAnalysis, $prepare_dir, $window, \@cases, \@controls);
        print " perl EXCAVATORDataAnalysis.pl $fileAnalysis --processors $processors --target $target --assembly $assembly --output $group_dir --mode $mode\n";
        system "perl $excavator_dir/EXCAVATORDataAnalysis.pl $fileAnalysis --processors $processors --target $target --assembly $assembly --output $group_dir --mode $mode";
	}
}

sub write_file_analysis {
    my $fileAnalysis = shift @_;
	my $prepare_dir  = shift @_;
	my $window       = shift @_;
	my $cases        = shift @_;
	my $controls     = shift @_;
	open OUT,">$fileAnalysis";
	foreach my $i(0..@{$cases} - 1) {
	    my $count = $i + 1;
	    print OUT "T$count $prepare_dir/$cases->[$i] $cases->[$i]\n";
	}
	foreach my $i(0..@{$controls} -1 ) {
	    my $count = $i + 1;
	    print OUT "C$count $prepare_dir/$controls->[$i] $controls->[$i]\n";
	}
	close OUT;
}

sub check_analysis_file {
    my $hash_group   = shift @_;
	my $report_dir   = shift @_;
    my $analysis_dir = shift @_;
	my $prepare_dir  = shift @_;
	my %hashCondition; 
    foreach my $suffix(keys %{$hash_group})
	{	
        my $is_control_ok = 1;
		my @controls = package::main::get_sample($hash_group->{$suffix}, 'control');
		my @cases    = package::main::get_sample($hash_group->{$suffix}, 'case');
		foreach my $control(@controls) {
		    my $control_prepare = "$prepare_dir/$control/RCNorm/$control.NRC.RData";
			if(not package::main::is_file_ok($control_prepare))
			{
			    $is_control_ok = 0;
				last;
			}
		}
		foreach my $case(@cases)
		{
		    my $prepare_file = "$prepare_dir/$case/RCNorm/$case.NRC.RData";
			my $final_file   = "$analysis_dir/$suffix/Results/$case/FastCallResults_$case.txt";
            # # 已完成该步骤
            if(package::main::is_file_ok($final_file))
            {
                $hashCondition{$suffix}{"Finish"}{$case} = "$final_file";
                next;            
            }
            # # 原始数据没问题
            if($is_control_ok == 1 and package::main::is_file_ok($prepare_file)){
                $hashCondition{$suffix}{"Good2Run"}{$case} = "$prepare_file";
                next;
            }
            # # 原始数据丢失
            $hashCondition{$suffix}{"Error"}{$case} = "$prepare_file control_prepare:$is_control_ok";
		}
		# # 写日志
        package::main::write_log("$report_dir/run.log", "EXCAVATORDataAnalysis", $hashCondition{$suffix}); 
    }
    return %hashCondition;
}

sub EXCAVATORDataPrepare {
	my $report_dir     = shift @_;
	my $result_dir     = shift @_;
	my $window         = shift @_;
	my $assembly       = shift @_;
	my $excavator_dir  = shift @_;
	my $target         = shift @_;
	my $output_dir     = shift @_;
	my $processors     = shift @_;
	my $config_samples = shift @_;
	print "\nProcessing EXCAVATORDataPrepare.pl ...\n";
	my $input_dir      = "$result_dir/input";
	my $prepare_dir    = "$result_dir/prepare";
	package::main::make_dir($input_dir);
	package::main::make_dir($prepare_dir);
	
	# # 文件检查
	my @samples      = check_prepare_file($report_dir, $prepare_dir, $output_dir, $config_samples);
	return if(!@samples);
	
	# # 生成配置文件
	my $file_prepare = "$input_dir/ExperimentalFilePrepare.w$window.txt";
	write_file_prepare($file_prepare, $prepare_dir, $output_dir, \@samples);
	
	# # 运行
	print "perl EXCAVATORDataPrepare.pl $file_prepare --processors $processors --target $target --assembly $assembly\n";
    system "perl $excavator_dir/EXCAVATORDataPrepare.pl $file_prepare --processors $processors --target $target --assembly $assembly";
}

sub write_file_prepare {
    my $file_prepare  = shift @_;
	my $prepare_dir   = shift @_;
	my $output_dir    = shift @_;
	my $samples       = shift @_;
	open OUT, ">$file_prepare";
	foreach my $sample(@{$samples})
	{
	    my $final_bam = "$output_dir/$sample/$sample\_final.bam";
		print OUT "$final_bam $prepare_dir/$sample $sample\n";
	}
	close OUT;
}

sub check_prepare_file {
    my $report_dir     = shift @_;
	my $prepare_dir    = shift @_;
	my $output_dir     = shift @_;
	my $config_samples = shift @_;
	my %hashCondition; 
    foreach my $sample(@{$config_samples})
    {
   	 	my $final_bam  = "$output_dir/$sample/$sample\_final.bam"; 
    	my $final_file = "$prepare_dir/$sample/RCNorm/$sample.NRC.RData";
        # # 已完成该步骤
        if(package::main::is_file_ok($final_file))
        {
            $hashCondition{"Finish"}{$sample} = "$final_file";
            next;            
        }
        # # 原始数据没问题
        if(package::main::is_file_ok($final_bam)){
            $hashCondition{"Good2Run"}{$sample} = "$final_bam";
            next;
        }
        # # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$final_bam";
    }
    # # 写日志
    package::main::write_log("$report_dir/run.log", "EXCAVATORDataPrepare", \%hashCondition); 
	my @samples = keys %{$hashCondition{"Good2Run"}};
	return @samples;
}

sub TargetPerla {
	my $result_dir    = shift @_;
	my $window        = shift @_;
	my $region_file   = shift @_;
	my $assembly      = shift @_;
	my $genome        = shift @_;
	my $excavator_dir = shift @_;
	my $target        = shift @_;
	my $target_data   = "$excavator_dir/TargetPerla.p/data/targets/$assembly/$target";
	my $prepare_data  = "$excavator_dir/EXCAVATORDataPrepare.p/data/targets/$assembly/$target";
	my $analysis_data = "$excavator_dir/EXCAVATORDataAnalysis.p/data/targets/$assembly/$target";
	print "\nProcessing TargetPerla.pl ...\n";
	if(not package::main::is_dir_ok($target_data)) {
        my $input_dir     = "$result_dir/input";
		package::main::make_dir($input_dir);
		my $source_target = "$input_dir/SourceTarget.txt";
		
		# # 生成配置文件
		write_source_target($source_target, $excavator_dir, $assembly, $genome);
        
		# # 运行
		print "perl TargetPerla.pl $source_target $region_file $target $window $assembly\n";
        system "perl $excavator_dir/TargetPerla.pl $source_target $region_file $target $window $assembly";
    }
}

sub write_source_target {
    my $source_target = shift @_;
	my $excavator_dir = shift @_;
	my $assembly      = shift @_;
	my $genome        = shift @_;	
	my $bw_dir        = "";
	if($assembly eq 'hg19')
	{
	    $bw_dir="$excavator_dir/TargetPerla.p/data/ucsc.hg19.bw";
	}
	elsif($assembly eq 'hg38')
	{
	    $bw_dir="$excavator_dir/TargetPerla.p/data/GCA_000001405.15_GRCh38.bw";
	}
	else
	{
	    return "please input extra .bw file!\n";
	}
	open OUT, ">$source_target";
	print OUT "$bw_dir $genome\n";
	close OUT;
}

1
