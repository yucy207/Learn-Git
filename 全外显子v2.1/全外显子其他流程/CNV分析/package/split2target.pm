package package::split2target;
use strict;
use warnings;

sub run {
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start split2target ".package::main::get_time()." ##########\n";
    my $isOK = package::check::split2target($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in split2target run para, please check\n" if($isOK == 0);
	
	my $hash_group     = $hashConfig->{'Group'};
	my $report_dir     = $hashConfig->{'Report'};
	my $region_orifile = $hashConfig->{'RealignBed'};
	my $DB_version     = $hashConfig->{'DB_version'};
	my $bedtools       = $hashPara->{"Soft"}{"bedtools"};  
	my $result_dir     = "$report_dir/excavator2";
	my $analysis_dir   = "$result_dir/analysis";
	my $correct_dir    = "$result_dir/result_final";
	my $region_file    = "$result_dir/merge_".(split /\//, $region_orifile)[-1];
	my %hash_cnvDB     = package::main::read_cnvDB($hashPara, $DB_version);
	
	# # 数据状态检测
	my %hashCondition; 
    foreach my $suffix(keys %{$hash_group}) {	
		my @cases    = package::main::get_sample($hash_group->{$suffix}, 'case');
		foreach my $case(@cases) {
			my $case_suffix_dir   = "$correct_dir/$case.$suffix";
			my $analysis_correct  = "$case_suffix_dir/$case.$suffix\_result.txt";
			my $bed_result2target = "$case_suffix_dir/$case.$suffix\_result2target.bed";
            # 已完成该步骤
            if(package::main::is_file_ok($bed_result2target)) {
                $hashCondition{"Finish"}{"$case.$suffix"} = "$bed_result2target";
                next;            
            }
			# # 原始数据没问题
            if(package::main::is_file_ok($analysis_correct, $region_file)) {
                $hashCondition{"Good2Run"}{"$case.$suffix"} = "$analysis_correct, $region_file";
                next;
            }
            # # 原始数据丢失
            $hashCondition{"Error"}{"$case.$suffix"} = "$analysis_correct";
		}
		# # 写日志
        package::main::write_log("$report_dir/run.log", "split2target", \%hashCondition); 
    }
	
	# # 确定target归属
	my @case_suffixes = keys %{$hashCondition{"Good2Run"}};
	if(@case_suffixes > 0) { # # 运行
        my $threshold = exists $hashPara->{"Process"}{"split2target"} ? $hashPara->{"Process"}{"split2target"} : 10;
		   $threshold = $hashConfig->{"Process_split2target"} if(exists $hashConfig->{"Process_split2target"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $case_suffix) = @_; package::main::process_bar_array($case_suffix, \@case_suffixes)});# 进度条
        foreach my $case_suffix(@case_suffixes) {
            $pm->start($case_suffix) and next;
			my $case_suffix_dir   = "$correct_dir/$case_suffix";
	        my $analysis_correct  = "$case_suffix_dir/$case_suffix\_result.txt";
	        my $bed_result2target = "$case_suffix_dir/$case_suffix\_result2target.bed";
	        system("$bedtools intersect -wa -loj -a $region_file -b $analysis_correct > $bed_result2target");
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
	
	# # 未过滤的所有target注释
	my @case_suffixes_belongs = (keys %{$hashCondition{"Finish"}}, keys %{$hashCondition{"Good2Run"}});
	my $case_suffixes_num     = @case_suffixes_belongs;
	my %all_target_infos      = get_belongs_target($correct_dir, \%hash_cnvDB, \@case_suffixes_belongs);
	annotate_target($hashPara, $hashConfig, $correct_dir, \%all_target_infos, $case_suffixes_num);
}

sub annotate_target {
    my $hashPara          = shift @_;
	my $hashConfig        = shift @_;
    my $correct_dir       = shift @_;
	my $all_target_infos  = shift @_;
	my $case_suffixes_num = shift @_;	
	my $species           = $hashConfig->{'Species'};
	my $AnnovarDIR        = $hashPara->{"Soft"}{"AnnovarDIR"};
	my $AnnovarBuild      = $hashPara->{$species}{"AnnovarBuild"};
	my $all_target_dir    = "$correct_dir/ALL_TARGETS";
	package::main::make_dir($all_target_dir);
	
	# # 生成library
	my $library_dir = "$all_target_dir/library";
	my $library     = "$library_dir/library";
	my $target_info = "$all_target_dir/all_target_infos.txt";
	package::main::make_dir($library_dir);
	open LIB, ">$library";
	open INFO, ">$target_info";
	foreach my $chr(1 .. 22, 'X') {
	    foreach my $start(sort {$a <=> $b} keys %{$all_target_infos->{$chr}}) {
		    foreach my $end(sort {$a <=> $b} keys %{$all_target_infos->{$chr}{$start}}) {
			    my $sample_infos      = $all_target_infos->{$chr}{$start}{$end}{'Sample_Infos'};
				my $gain_sample       = exists $all_target_infos->{$chr}{$start}{$end}{'Gain_Sample'} ? $all_target_infos->{$chr}{$start}{$end}{'Gain_Sample'} : "";
				my $loss_sample       = exists $all_target_infos->{$chr}{$start}{$end}{'Loss_Sample'} ? $all_target_infos->{$chr}{$start}{$end}{'Loss_Sample'} : "";
				my $gain_sample_num   = exists $all_target_infos->{$chr}{$start}{$end}{'Gain_Sample_num'} ? $all_target_infos->{$chr}{$start}{$end}{'Gain_Sample_num'} : 0;
				my $loss_sample_num   = exists $all_target_infos->{$chr}{$start}{$end}{'Loss_Sample_num'} ? $all_target_infos->{$chr}{$start}{$end}{'Loss_Sample_num'} : 0;
				my $cnvDB_Gain_Freq   = $all_target_infos->{$chr}{$start}{$end}{'cnvDB_Gain_Freq'};
				my $cnvDB_Loss_Freq   = $all_target_infos->{$chr}{$start}{$end}{'cnvDB_Loss_Freq'};
				my $cnvDB_Gain_Sample = $all_target_infos->{$chr}{$start}{$end}{'cnvDB_Gain_Sample'};
				my $cnvDB_Loss_Sample = $all_target_infos->{$chr}{$start}{$end}{'cnvDB_Loss_Sample'};
				my $Gain_Freq         = $gain_sample_num / $case_suffixes_num;
				my $Loss_Freq         = $loss_sample_num / $case_suffixes_num;
				print LIB "$chr\t$start\t$end\t0\t-\n";
				print INFO "$chr\t$start\t$end\t$Gain_Freq\t$Loss_Freq\t$gain_sample\t$loss_sample\t$cnvDB_Gain_Freq\t$cnvDB_Loss_Freq\t$cnvDB_Gain_Sample\t$cnvDB_Loss_Sample\t$sample_infos\n";
		    }
		}
	}
	close LIB;
	close INFO;
	
	# # 注释	
	my %hashAnnotationList = package::annotation_db::get_annotation(); # # 获取注释列表
	foreach my $anno_model(keys %hashAnnotationList) {
        next if($hashAnnotationList{$anno_model}{'From'} ne 'Annovar');
        my $AnnovarDB = $hashAnnotationList{$anno_model}{'DBDir'}{$species};
        system "$AnnovarDIR/$hashAnnotationList{$anno_model}{'Para'} --buildver $AnnovarBuild $library $AnnovarDB";
    }
}

sub get_belongs_target {
    my $correct_dir           = shift @_;
	my $hash_cnvDB            = shift @_;
	my $case_suffixes_belongs = shift @_;
	my %all_target_infos;
	foreach my $case_suffix_belongs(@$case_suffixes_belongs) {
	    my $case_suffix_dir   = "$correct_dir/$case_suffix_belongs";
		my $bed_result2target = "$case_suffix_dir/$case_suffix_belongs\_result2target.bed";
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
			if(not exists $all_target_infos{$target_chr}{$target_start}{$target_end}) {	
				$all_target_infos{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Freq'}   = $cnvDB_Gain_Freq;
				$all_target_infos{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Freq'}   = $cnvDB_Loss_Freq;
				$all_target_infos{$target_chr}{$target_start}{$target_end}{'cnvDB_Gain_Sample'} = $cnvDB_Gain_Sample;
				$all_target_infos{$target_chr}{$target_start}{$target_end}{'cnvDB_Loss_Sample'} = $cnvDB_Loss_Sample;				
			}
			# # 其他target CNV信息获取
			$all_target_infos{$target_chr}{$target_start}{$target_end}{'Sample_Infos'}      .= "$case_suffix_belongs:$cnv:$type,"; 
			if($type eq 'gain') {
			    $all_target_infos{$target_chr}{$target_start}{$target_end}{'Gain_Sample'}   .= "$case_suffix_belongs,";
				$all_target_infos{$target_chr}{$target_start}{$target_end}{'Gain_Sample_num'} ++;
			}
			elsif($type eq 'loss') {
			    $all_target_infos{$target_chr}{$target_start}{$target_end}{'Loss_Sample'}   .= "$case_suffix_belongs,";
				$all_target_infos{$target_chr}{$target_start}{$target_end}{'Loss_Sample_num'} ++;
			}
		}
		close BED;
	}
	return %all_target_infos;
}

1
