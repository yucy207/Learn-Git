package package::target_filter_combine;
use strict;
use warnings;

sub run {
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start target_filter_combine ".package::main::get_time()." ##########\n";
    my $isOK = package::check::target_filter_combine($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in target_filter_combine run para, please check\n" if($isOK == 0);
	
	my $hash_group     = $hashConfig->{'Group'};
	my $report_dir     = $hashConfig->{'Report'};
	my $region_orifile = $hashConfig->{'RealignBed'};
	my $species        = $hashConfig->{'Species'};
	my $AnnovarDIR     = $hashPara->{"Soft"}{"AnnovarDIR"};
	my $AnnovarBuild   = $hashPara->{$species}{"AnnovarBuild"};
	my $result_dir     = "$report_dir/excavator2";
	my $correct_dir    = "$result_dir/result_final";	
	my $all_target_dir = "$correct_dir/ALL_TARGETS";
	my $target_info    = "$all_target_dir/all_target_infos.txt";
	my $combine_dir    = "$correct_dir/COMBINE_REGION";
	my $region_combine = "$combine_dir/combined_region_infos.txt";
	my $region_file    = "$result_dir/merge_".(split /\//, $region_orifile)[-1];
	package::main::make_dir($combine_dir);
	
	# # 数据状态检测
	my %hashCondition; 
    # 已完成该步骤
    if(package::main::is_file_ok($region_combine)) {
        $hashCondition{"Finish"}{"ALL_TARGETS"} = "$region_combine";          
    }
	# 原始数据没问题
    elsif(package::main::is_file_ok($target_info, $region_file)) {
        $hashCondition{"Good2Run"}{"ALL_TARGETS"} = "$target_info, $region_file";
    }
	# 原始数据丢失
	else {  
        $hashCondition{"Error"}{"ALL_TARGETS"} = "$target_info";
	}
	
	# # 写日志
    package::main::write_log("$report_dir/run.log", "target_filter_combine", \%hashCondition); 
	
	if(not exists $hashCondition{"Good2Run"}) {
	    print "[Note] Process None, for reProcess Delete Result\n";
		return;
	}
	
	# # target筛选
	my $max_mutation_freq  = 0.1; # # gain_freq 或 loss_freq 大于该值的target被过滤掉
	my $target_info_filter = "$combine_dir/filter_target_infos.txt";
	my ($hash_sample_infos, $hash_all_infos) = filter_target($target_info, $target_info_filter, $max_mutation_freq);
	
	# # target合并
	my @case_suffixes  = package::main::get_all_case_suffix($hashConfig, $correct_dir);
	my %combine_region = comebine_target($hash_sample_infos, $region_file, \@case_suffixes);
	get_combine_infos(\%combine_region, $hash_all_infos);
	my $library_dir    = "$combine_dir/library";
	my $library        = "$library_dir/library";
	package::main::make_dir($library_dir);
	write_combine_region($all_target_dir, \%combine_region, $region_combine, $library);
	
	# # 合并后区域注释	
	my %hashAnnotationList = package::annotation_db::get_annotation(); # # 获取注释列表
	foreach my $anno_model(keys %hashAnnotationList) {
        next if($hashAnnotationList{$anno_model}{'From'} ne 'Annovar');
        my $AnnovarDB = $hashAnnotationList{$anno_model}{'DBDir'}{$species};
        system("$AnnovarDIR/$hashAnnotationList{$anno_model}{'Para'} --buildver $AnnovarBuild $library $AnnovarDB");
    }
}

sub write_combine_region {
    my $all_target_dir = shift @_;
	my $combine_region = shift @_;
	my $region_combine = shift @_;
	my $library        = shift @_;
	open COM, ">$region_combine";
	open LIB, ">$library";
	foreach my $chr(sort {$a cmp $b} keys %{$combine_region}) {
	    foreach my $start(sort {$a <=> $b} keys %{$combine_region->{$chr}}) {
		    foreach my $end(sort {$a <=> $b} keys %{$combine_region->{$chr}{$start}}) {
				my $Gain_Freq         = $combine_region->{$chr}{$start}{$end}{'Infos'}{'Gain_Freq'};
				my $Loss_Freq         = $combine_region->{$chr}{$start}{$end}{'Infos'}{'Loss_Freq'};
				my $gain_sample       = $combine_region->{$chr}{$start}{$end}{'Infos'}{'Gain_Sample'};
				my $loss_sample       = $combine_region->{$chr}{$start}{$end}{'Infos'}{'Loss_Sample'};
				my $cnvDB_Gain_Freq   = $combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Gain_Freq'};
				my $cnvDB_Loss_Freq   = $combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Loss_Freq'};
				my $cnvDB_Gain_Sample = $combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Gain_Sample'};
				my $cnvDB_Loss_Sample = $combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Loss_Sample'};
				my $sample_infos      = $combine_region->{$chr}{$start}{$end}{'Infos'}{'sample_infos'};
				my $target_num        = $combine_region->{$chr}{$start}{$end}{'Infos'}{'target_num'};
				print COM "$chr\t$start\t$end\t$Gain_Freq\t$Loss_Freq\t$gain_sample\t$loss_sample\t$cnvDB_Gain_Freq\t$cnvDB_Loss_Freq\t$cnvDB_Gain_Sample\t$cnvDB_Loss_Sample\t$sample_infos\t$target_num\n";
				print LIB "$chr\t$start\t$end\t0\t-\n";
			}
		}
	}
	close COM;
	close LIB;
}

sub get_combine_infos {
    my $combine_region = shift @_;
	my $hash_all_infos = shift @_;
	foreach my $chr(keys %{$combine_region}) {
	    foreach my $start(keys %{$combine_region->{$chr}}) {
		    foreach my $end(keys %{$combine_region->{$chr}{$start}}) {
			    my ($Gain_Freq_sum, $Loss_Freq_sum, $cnvDB_Gain_Freq_sum, $cnvDB_Loss_Freq_sum) = (0, 0, 0, 0);
				my (%cnvDB_Gain_Samples, %cnvDB_Loss_Samples);
				# # 汇总样本信息
				foreach my $case_suffix(keys %{$combine_region->{$chr}{$start}{$end}}) {
				    next if($case_suffix eq 'Infos');
					my $type = $combine_region->{$chr}{$start}{$end}{$case_suffix}{'type'};
			        my $cnv  = $combine_region->{$chr}{$start}{$end}{$case_suffix}{'cnv'};
					$combine_region->{$chr}{$start}{$end}{'Infos'}{'sample_infos'} .= "$case_suffix:$cnv:$type,";
				}			
				foreach my $target_pos(split /,/, $combine_region->{$chr}{$start}{$end}{'Infos'}{'target_pos'}) {
				    my ($target_chr, $target_start, $target_end) = split /:/, $target_pos;
					# # 合并前target信息读取
					my $Gain_Freq         = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'Gain_Freq'};        
	                my $Loss_Freq         = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'Loss_Freq'};        
	                my $gain_sample       = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'Gain_Sample'};      
	                my $loss_sample       = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'Loss_Sample'};      
	                my $cnvDB_Gain_Freq   = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Freq'};  
	                my $cnvDB_Loss_Freq   = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Freq'};  
	                my $cnvDB_Gain_Sample = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Sample'};
	                my $cnvDB_Loss_Sample = $hash_all_infos->{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Sample'};
					# # 各突变频率求和，用于求均值
					$Gain_Freq_sum       += $Gain_Freq;
					$Loss_Freq_sum       += $Loss_Freq;
					$cnvDB_Gain_Freq_sum += $cnvDB_Gain_Freq;
					$cnvDB_Loss_Freq_sum += $cnvDB_Loss_Freq;
					# # CNV数据库突变样本汇总、去重（项目样本合并的区域里样本一定相同，不必处理）
					map{$cnvDB_Gain_Samples{$_} ++} (split /,/, $cnvDB_Gain_Sample);
					map{$cnvDB_Loss_Samples{$_} ++} (split /,/, $cnvDB_Loss_Sample);
					$combine_region->{$chr}{$start}{$end}{'Infos'}{'Gain_Sample'}   = $gain_sample;
				    $combine_region->{$chr}{$start}{$end}{'Infos'}{'Loss_Sample'}   = $loss_sample;
				}
				my @cnvDB_Gain_Samples = sort keys %cnvDB_Gain_Samples;
				my @cnvDB_Loss_Samples = sort keys %cnvDB_Loss_Samples;
				my $target_num = $combine_region->{$chr}{$start}{$end}{'Infos'}{'target_num'};
				$combine_region->{$chr}{$start}{$end}{'Infos'}{'Gain_Freq'}         = $Gain_Freq_sum / $target_num;
				$combine_region->{$chr}{$start}{$end}{'Infos'}{'Loss_Freq'}         = $Loss_Freq_sum / $target_num;
				$combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Gain_Freq'}   = $cnvDB_Gain_Freq_sum / $target_num;
				$combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Loss_Freq'}   = $cnvDB_Loss_Freq_sum / $target_num;
				$combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Gain_Sample'} = join ",", @cnvDB_Gain_Samples;
				$combine_region->{$chr}{$start}{$end}{'Infos'}{'cnvDB_Loss_Sample'} = join ",", @cnvDB_Loss_Samples;
			}
		}
	}
}

sub comebine_target {
    my $hash_sample_infos = shift @_;
	my $region_file       = shift @_;
	my $case_suffixes     = shift @_;
	my %last_info;
	my %combine_region;
	my ($last_chr, $last_start, $last_end, $target_num) = (0, 0, 0, 0);
	open REGION, $region_file;
	while(<REGION>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($now_chr, $now_start, $now_end, $tmp) = split /\t/, $_;
		if(exists $hash_sample_infos->{$now_chr}{$now_start}{$now_end}) {
		    if($last_start == 0) { # # 初始化 或 上一个区域不符合要求
			    $target_num = 1;
				($last_chr, $last_start, $last_end) = ($now_chr, $now_start, $now_end);
				initialise_sampe_infos($now_chr, $now_start, $now_end, \%last_info, $hash_sample_infos, $case_suffixes);
			}
			else { # # 上一区域符合要求
				my $if_combine = decide_if_combine($last_chr, $now_chr, $now_start, $now_end, \%last_info, $hash_sample_infos, $case_suffixes);
				if($if_combine == 1) { # # 该target可以合并
				    $last_start = ($last_start == 0) ? $now_start : $last_start; # # 用于染色体变化的情况
					$last_end   = $now_end;
					$target_num ++;
					map{$last_info{$_}{'cnv_sum'}     += $hash_sample_infos->{$now_chr}{$now_start}{$now_end}{$_}{'cnv'} if(exists $hash_sample_infos->{$now_chr}{$now_start}{$now_end}{$_})} @$case_suffixes;
				    $last_info{'Infos'}{'target_pos'} .= "$now_chr:$now_start:$now_end,";
				}
				else { # # 该target不可以合并
				    save_last_region($last_chr, $last_start, $last_end, $target_num, \%last_info, \%combine_region, $case_suffixes);
					$target_num = 1;
					($last_chr, $last_start, $last_end) = ($now_chr, $now_start, $now_end);
					initialise_sampe_infos($now_chr, $now_start, $now_end, \%last_info, $hash_sample_infos, $case_suffixes);
				}
			}
		}
		elsif((not exists $hash_sample_infos->{$now_chr}{$now_start}{$now_end}) and ($last_start != 0)) { # # 该target(没有突变 或 没有通过过滤) 且 上一个target也不符合条件
		    save_last_region($last_chr, $last_start, $last_end, $target_num, \%last_info, \%combine_region, $case_suffixes);
			$last_start = 0;
		}
	}
	save_last_region($last_chr, $last_start, $last_end, $target_num, \%last_info, \%combine_region, $case_suffixes) if($last_start != 0);
	close REGION;
	return %combine_region;
}

sub decide_if_combine {
    my $last_chr          = shift @_;
    my $now_chr           = shift @_;
	my $now_start         = shift @_;
	my $now_end           = shift @_;
	my $last_info         = shift @_;
	my $hash_sample_infos = shift @_;
	my $case_suffixes     = shift @_;
	my $if_combine        = 1;
	$if_combine           = 0 if($now_chr ne $last_chr); # # 当前染色体与上一target染色体不一致，即不合并
	foreach my $case_suffix(@$case_suffixes) { # # 任一样本组的CNV类型与上一target的CNV类型不一致，即不合并
	    my $type    = 'normal';
		$type       = $hash_sample_infos->{$now_chr}{$now_start}{$now_end}{$case_suffix}{'type'} if(exists $hash_sample_infos->{$now_chr}{$now_start}{$now_end}{$case_suffix});
		$if_combine = 0 if($type ne $last_info->{$case_suffix}{'type'});
	}
	return $if_combine;
}

sub save_last_region {
    my $last_chr       = shift @_;
	my $last_start     = shift @_;
	my $last_end       = shift @_;
	my $target_num     = shift @_;
	my $last_info      = shift @_;
	my $combine_region = shift @_;
	my $case_suffixes  = shift @_;
	
	foreach my $case_suffix(@$case_suffixes) {
	    if($last_info->{$case_suffix}{'type'} ne 'normal') {
		    my $type       = $last_info->{$case_suffix}{'type'};
			my $cnv_sum    = $last_info->{$case_suffix}{'cnv_sum'};
			my $mean_cnv   = $cnv_sum / $target_num;
			$combine_region->{$last_chr}{$last_start}{$last_end}{$case_suffix}{'type'}  = $type;
			$combine_region->{$last_chr}{$last_start}{$last_end}{$case_suffix}{'cnv'}   = $mean_cnv;
		}
	}
	my $target_pos = $last_info->{'Infos'}{'target_pos'};
	$combine_region->{$last_chr}{$last_start}{$last_end}{'Infos'}{'target_num'} = $target_num;
	$combine_region->{$last_chr}{$last_start}{$last_end}{'Infos'}{'target_pos'} = $target_pos;
	$last_info->{'Infos'}{'target_pos'} = ""; # # 清空上一区域的target位置
}

sub initialise_sampe_infos {
    my $now_chr           = shift @_;
	my $now_start         = shift @_;
	my $now_end           = shift @_;
	my $last_info         = shift @_;
	my $hash_sample_infos = shift @_;
	my $case_suffixes     = shift @_;
	foreach my $case_suffix(@$case_suffixes) {
	    my $type = 'normal';
		if(exists $hash_sample_infos->{$now_chr}{$now_start}{$now_end}{$case_suffix})	{
		    $type = $hash_sample_infos->{$now_chr}{$now_start}{$now_end}{$case_suffix}{'type'};
			$last_info->{$case_suffix}{'cnv_sum'} = $hash_sample_infos->{$now_chr}{$now_start}{$now_end}{$case_suffix}{'cnv'};
		}
		$last_info->{$case_suffix}{'type'} = $type;
	}
	$last_info->{'Infos'}{'target_pos'} .= "$now_chr:$now_start:$now_end,";
}

sub filter_target {
    my $target_info = shift @_;
	my $target_info_filter = shift @_;
	my $max_mutation_freq  = shift @_;
	my %hash_sample_infos;
	my %hash_all_infos;
	open INFO, $target_info;
	open FILTER, ">$target_info_filter";
	while(my $line = <INFO>) {
	    next if($line !~ /\w/);
		$line =~ s/[\r\n]//g;
		my (
		    $chr, $start, $end,  
			$Gain_Freq, $Loss_Freq, $gain_sample, $loss_sample,
			$cnvDB_Gain_Freq, $cnvDB_Loss_Freq, $cnvDB_Gain_Sample, $cnvDB_Loss_Sample,
			$sample_infos
			) = split /\t/, $line;
		# # CNV数据路频率过滤
		next if(($cnvDB_Gain_Freq > $max_mutation_freq) or ($cnvDB_Loss_Freq > $max_mutation_freq));
		# # 记录通过过滤的target信息
		$hash_all_infos{$chr}{$start}{$end}{'Gain_Freq'}         = $Gain_Freq;
		$hash_all_infos{$chr}{$start}{$end}{'Loss_Freq'}         = $Loss_Freq;
		$hash_all_infos{$chr}{$start}{$end}{'Gain_Sample'}       = $gain_sample;
		$hash_all_infos{$chr}{$start}{$end}{'Loss_Sample'}       = $loss_sample;
		$hash_all_infos{$chr}{$start}{$end}{'cnvDB_Gain_Freq'}   = $cnvDB_Gain_Freq;
		$hash_all_infos{$chr}{$start}{$end}{'cnvDB_Loss_Freq'}   = $cnvDB_Loss_Freq;
		$hash_all_infos{$chr}{$start}{$end}{'cnvDB_Gain_Sample'} = $cnvDB_Gain_Sample;
		$hash_all_infos{$chr}{$start}{$end}{'cnvDB_Loss_Sample'} = $cnvDB_Loss_Sample;
		print FILTER "$line\n";
		# # 获取通过过滤的样本信息，用于后续区域合并
        get_sample_infos($chr, $start, $end, $sample_infos, \%hash_sample_infos);		
	}
	close INFO;
	close FILTER;
	return(\%hash_sample_infos, \%hash_all_infos);
}

sub get_sample_infos {
    my $chr                = shift @_;
	my $start              = shift @_;
	my $end                = shift @_;
	my $sample_infos       = shift @_;
	my $hash_sample_infos  = shift @_;
	my @sample_infos = split /,/, $sample_infos;
	foreach my $sample_info(@sample_infos) {
	    my ($case_suffix, $cnv, $type) = split /:/, $sample_info;
		$hash_sample_infos->{$chr}{$start}{$end}{$case_suffix}{'cnv'}  = $cnv;
		$hash_sample_infos->{$chr}{$start}{$end}{$case_suffix}{'type'} = $type;		
	}
}


1