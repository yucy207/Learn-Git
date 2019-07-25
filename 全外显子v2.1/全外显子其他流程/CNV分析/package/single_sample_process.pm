package package::single_sample_process;
use strict;
use warnings;

sub run {
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start single_sample_process ".package::main::get_time()." ##########\n";
    my $isOK = package::check::single_sample_process($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in single_sample_process run para, please check\n" if($isOK == 0);
		
	my $hash_group     = $hashConfig->{'Group'};
	my $report_dir     = $hashConfig->{'Report'};
	my $result_dir     = "$report_dir/excavator2";
	my $correct_dir    = "$result_dir/result_final";
	
	# # 数据状态检测
	my %hashCondition; 
    foreach my $suffix(keys %{$hash_group}) {
		my @cases    = package::main::get_sample($hash_group->{$suffix}, 'case');
		foreach my $case(@cases) {
			my $case_suffix_dir    = "$correct_dir/$case.$suffix";
	        my $analysis_correct   = "$case_suffix_dir/$case.$suffix\_result.txt";
			my $bed_result2target  = "$case_suffix_dir/$case.$suffix\_result2target.bed";
	        my $belongs_statistics = "$case_suffix_dir/$case.$suffix\_belongs.statistics";
            # 已完成该步骤
            if(package::main::is_file_ok($belongs_statistics)) {
                $hashCondition{"Finish"}{"$case.$suffix"} = "$belongs_statistics";
                next;            
            }
			# # 原始数据没问题
            if(package::main::is_file_ok($analysis_correct, $bed_result2target)) {
                $hashCondition{"Good2Run"}{"$case.$suffix"} = "$analysis_correct, $bed_result2target";
                next;
            }
            # # 原始数据丢失
            $hashCondition{"Error"}{"$case.$suffix"} = "$analysis_correct, $bed_result2target";
		}
		# # 写日志
        package::main::write_log("$report_dir/run.log", "single_sample_process", \%hashCondition); 
    }
	
	# # 运行
	my @case_suffixes = keys %{$hashCondition{"Good2Run"}};
	if(@case_suffixes > 0) { # # 运行
        my $threshold = exists $hashPara->{"Process"}{"single_sample_process"} ? $hashPara->{"Process"}{"single_sample_process"} : 5;
		   $threshold = $hashConfig->{"Process_single_sample_process"} if(exists $hashConfig->{"Process_single_sample_process"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $case_suffix) = @_; package::main::process_bar_array($case_suffix, \@case_suffixes)});# 进度条
        foreach my $case_suffix(@case_suffixes) {
            $pm->start($case_suffix) and next;
			sub_run($hashPara, $hashConfig, $case_suffix);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
}

sub sub_run {
    my $hashPara     = shift @_;
	my $hashConfig   = shift @_;
	my $case_suffix  = shift @_;
	my $report_dir   = $hashConfig->{'Report'};
	my $species      = $hashConfig->{'Species'};
	my $DB_version   = $hashConfig->{'DB_version'};
	my $AnnovarDIR   = $hashPara->{"Soft"}{"AnnovarDIR"};
	my $AnnovarBuild = $hashPara->{$species}{"AnnovarBuild"};
	my $result_dir   = "$report_dir/excavator2";
	my $correct_dir  = "$result_dir/result_final";
	my %hash_cnvDB   = package::main::read_cnvDB($hashPara, $DB_version);
	
	# # 单样本区域注释
	my $case_suffix_dir  = "$correct_dir/$case_suffix";
	my $analysis_correct = "$case_suffix_dir/$case_suffix\_result.txt";
	my $library_dir      = "$case_suffix_dir/library";
	my $library          = "$library_dir/library";
	package::main::make_dir($library_dir);
	generate_library($analysis_correct, $library);
	
	my %hashAnnotationList = package::annotation_db::get_annotation(); # # 获取注释列表
	foreach my $anno_model(keys %hashAnnotationList) {
        next if($hashAnnotationList{$anno_model}{'From'} ne 'Annovar');
        my $AnnovarDB = $hashAnnotationList{$anno_model}{'DBDir'}{$species};
        system("$AnnovarDIR/$hashAnnotationList{$anno_model}{'Para'} --buildver $AnnovarBuild $library $AnnovarDB");
    }
	
	# # 单样本信息统计
	my %region_target_num = get_belongs_sample($correct_dir, \%hash_cnvDB, $case_suffix);
	write_belongs_statistics($correct_dir, \%region_target_num, $case_suffix);
}

sub write_belongs_statistics {
    my $correct_dir       = shift @_;
	my $region_target_num = shift @_;
	my $case_suffix       = shift @_;
	my $case_suffix_dir    = "$correct_dir/$case_suffix";
	my $analysis_correct   = "$case_suffix_dir/$case_suffix\_result.txt";
	my $belongs_statistics = "$case_suffix_dir/$case_suffix\_belongs.statistics";
	open CORRECT, $analysis_correct;
	open STAT, ">$belongs_statistics";
	while(<CORRECT>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($chr, $start, $end, $tmp) = split /\t/, $_, 4;
		my $target_num          = $region_target_num->{$chr}{$start}{$end}{'target_num'};
		my $cnvDB_Gain_Freq_sum = $region_target_num->{$chr}{$start}{$end}{'cnvDB_Gain_Freq'};
		my $cnvDB_Loss_Freq_sum = $region_target_num->{$chr}{$start}{$end}{'cnvDB_Loss_Freq'};
		my $cnvDB_Gain_Freq     = $cnvDB_Gain_Freq_sum / $target_num;
		my $cnvDB_Loss_Freq     = $cnvDB_Loss_Freq_sum / $target_num;	
		print STAT "$chr\t$start\t$end\t$tmp\t$target_num\t$cnvDB_Gain_Freq\t$cnvDB_Loss_Freq\n";
	}
	close CORRECT;
	close STAT;
}


sub get_belongs_sample {
    my $correct_dir = shift @_;
	my $hash_cnvDB  = shift @_;
	my $case_suffix = shift @_;
	my %region_target_num;
	my $case_suffix_dir   = "$correct_dir/$case_suffix";
	my $bed_result2target = "$case_suffix_dir/$case_suffix\_result2target.bed";
	open BED, $bed_result2target;
	while(<BED>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($target_chr, $target_start, $target_end, $region_chr, $region_start, $region_end, $cnv, $type, $log2, $call_prob) = split /\t/, $_;
		next if($region_chr !~ /\w/); # # 该target不属于任何区域
		# # 获取CNV数据库注释
		my $cnvDB_Gain_Freq   = exists $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Freq'} ? $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Freq'} : 0;
		my $cnvDB_Loss_Freq   = exists $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Freq'} ? $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Freq'} : 0;
		my $cnvDB_Gain_Sample = exists $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Sample'} ? $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Sample'} : "";
        my $cnvDB_Loss_Sample = exists $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Sample'} ? $hash_cnvDB->{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Sample'} : "";
		# # 各结果区域包含target数目统计
		$region_target_num{$region_chr}{$region_start}{$region_end}{'target_num'} ++;
		$region_target_num{$region_chr}{$region_start}{$region_end}{'cnvDB_Gain_Freq'} += $cnvDB_Gain_Freq;
		$region_target_num{$region_chr}{$region_start}{$region_end}{'cnvDB_Loss_Freq'} += $cnvDB_Loss_Freq;
	}
	close BED;
	return %region_target_num;
}

sub generate_library {
    my $analysis_correct = shift @_;
	my $library          = shift @_;
	open CORRECT, $analysis_correct;
	open LIB, ">$library";
	while(<CORRECT>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($chr, $start, $end, $tmp) = split /\t/, $_;
		print LIB "$chr\t$start\t$end\t0\t-\n";
	}
	close CORRECT;
	close LIB;
}

1