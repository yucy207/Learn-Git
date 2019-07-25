package package::prepare_for_plot_n_snv;
use strict;
use warnings;

sub run {
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start prepare_for_plot_n_snv ".package::main::get_time()." ##########\n";
    my $isOK = package::check::prepare_for_plot_n_snv($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in prepare_for_plot_n_snv run para, please check\n" if($isOK == 0);
	
	my $hash_group     = $hashConfig->{'Group'};
	my $report_dir     = $hashConfig->{'Report'};
	my $output_dir     = $hashConfig->{'Output'};
	my $species        = $hashConfig->{'Species'};
	my $region_orifile = $hashConfig->{'RealignBed'};
	my $DB_version     = $hashConfig->{'DB_version'};
	
	my $AnnovarDIR     = $hashPara->{"Soft"}{"AnnovarDIR"};
	my $genome         = $hashPara->{$species}{'Genome'};
	my $region_dict    = $hashPara->{$species}{'Dict'};
	my $AnnovarBuild   = $hashPara->{$species}{"AnnovarBuild"};
	
	my $result_dir        = "$report_dir/excavator2";
	my $correct_dir       = "$result_dir/result_final";
	my $status_dir        = "$result_dir/status";
	my $plot_bed_dir      = "$result_dir/bed_for_plot";
	my $snp_filter_DB_dir = "$result_dir/snp_filter_DB_dir";
	my %hash_cnvDB        = package::main::read_cnvDB($hashPara, $DB_version);
    package::main::make_dir($result_dir);
	package::main::make_dir($status_dir);
	package::main::make_dir($plot_bed_dir);
	package::main::make_dir($snp_filter_DB_dir);
	
	# # 进行tag统计	
    generate_tag($hashConfig, $hashPara, $report_dir, $output_dir, $result_dir, $status_dir, $hash_group, $region_dict, $region_orifile, $genome);
	
	my @case_suffixes = package::main::get_all_case_suffix($hashConfig, $correct_dir);   
	my %hash_all_sample_region;
	foreach my $case_suffix(@case_suffixes) {
	    my %sample_pos;
		my $region_bed_filter = "$plot_bed_dir/$case_suffix";
        my $original_data_dir = "$region_bed_filter/original_data";
		package::main::make_dir($region_bed_filter);
	    package::main::make_dir($original_data_dir);
		# # 单样本信息统计
	    my %region_target_num = get_belongs_sample($correct_dir, \%hash_cnvDB, $case_suffix);
		# # 区域过滤
	    region_filter($correct_dir, $original_data_dir, \%region_target_num, $case_suffix, \%hash_all_sample_region, \%sample_pos);
		# # 输出绘图原始文件
		write_region_bed($original_data_dir, \%sample_pos);
    }
    # # 输出SNV过滤文件
	write_bed_for_filter($snp_filter_DB_dir, $AnnovarBuild, \%hash_all_sample_region);
	all_snv_annotation($AnnovarDIR, $AnnovarBuild, $report_dir, $snp_filter_DB_dir, $result_dir);
}

sub all_snv_annotation {
	my $AnnovarDIR        = shift @_;
	my $AnnovarBuild      = shift @_;
	my $report_dir        = shift @_;
	my $snp_filter_DB_dir = shift @_;
	my $result_dir        = shift @_;
	my $all_library       = "$report_dir/library/library";
	my $database_snv      = "$report_dir/database.snv";
	my $all_library_new   = "$snp_filter_DB_dir/library";
	die "wait until database.snv done\n" if(not -e $database_snv);
	my %hash_database       = read_database_snv($database_snv);
	system("ln -s $all_library $all_library_new") if(not -e $all_library_new);
	# # 注释
	system("$AnnovarDIR/annotate_variation.pl -regionanno -dbtype cnv_gain_info --buildver $AnnovarBuild $all_library_new $snp_filter_DB_dir");
	system("$AnnovarDIR/annotate_variation.pl -regionanno -dbtype cnv_loss_info --buildver $AnnovarBuild $all_library_new $snp_filter_DB_dir");
	# # 生成annotation文件
	my %hash_anno;
	my $gain_anno_file = "$all_library_new.$AnnovarBuild\_cnv_gain_info";
	my $loss_anno_file = "$all_library_new.$AnnovarBuild\_cnv_loss_info";
	read_gain_loss_annotation(\%hash_anno, $gain_anno_file, 'CNV_gain');
	read_gain_loss_annotation(\%hash_anno, $loss_anno_file, 'CNV_loss');
	my $gain_annotation_file = "$result_dir/CNV_gain_annotation.txt";
	my $loss_annotation_file = "$result_dir/CNV_loss_annotation.txt";
	write_annotation(\%hash_database, \%hash_anno, $gain_annotation_file, 'CNV_gain');
	write_annotation(\%hash_database, \%hash_anno, $loss_annotation_file, 'CNV_loss');
}

sub write_annotation {
    my $hash_database   = shift @_;
	my $hash_anno       = shift @_;
	my $annotation_file = shift @_;
	my $type            = shift @_;
	open FILE, ">$annotation_file";
	foreach my $title(keys %{$hash_database}) {
	    my ($chr, $pos, $ref) = split /\|/, $title;
		my $start = (split /\-/, $pos)[0];
		foreach my $alt(keys %{$hash_database->{$title}}) {
		    if(exists $hash_anno->{$type}{$chr}{$start}{$ref}{$alt}) {
			    my $to_annotate = $hash_anno->{$type}{$chr}{$start}{$ref}{$alt};
				print FILE "$title\t$alt\t$type\t$to_annotate\n";
			}
		}
	}
	close FILE;
}


sub read_gain_loss_annotation {
    my $hash_anno = shift @_;
	my $anno_file = shift @_;
	my $type      = shift @_;
	open ANNO, $anno_file;
	while(<ANNO>) {
	   next if($_ !~ /\w/);
	   $_ =~ s/[\r\n]//g;
	   my ($mark, $infos, $chr, $start, $end, $ref, $alt, $tmp)= split /\t/, $_;
	   my $to_annotate = (split /=/, $infos)[1];
	   $hash_anno->{$type}{$chr}{$start}{$ref}{$alt} = $to_annotate;
	}
}

sub read_database_snv {
    my $database_snv = shift @_;
	my %hash_database;
	open DATA, $database_snv;
	while(<DATA>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my @infos = split /\t/, $_;
		my $title = $infos[0];
		my ($chr, $pos, $ref) = split /\|/, $title;
		$pos = (split /\-/, $pos)[0];
		my $alt   = $infos[5];
		my $freq  = $infos[3];
		$hash_database{$title}{$alt} ++;
	}
	close DATA;
	return %hash_database;
}

sub write_bed_for_filter {
     my $snp_filter_DB_dir      = shift @_;
	 my $AnnovarBuild           = shift @_;
	 my $hash_all_sample_region = shift @_;
	 foreach my $type(keys %{$hash_all_sample_region}) {
	     my $count          = 1;
		 my $bed_for_filter = "$snp_filter_DB_dir/$AnnovarBuild\_cnv_$type\_info.txt";
		 open BED, ">$bed_for_filter";
		 foreach my $chr(1..22,'X') {
		     foreach my $start(sort {$a <=> $b} keys %{$hash_all_sample_region->{$type}{$chr}}) {
			     foreach my $end(sort {$a <=> $b} keys %{$hash_all_sample_region->{$type}{$chr}{$start}}) {
				     my $info = $hash_all_sample_region->{$type}{$chr}{$start}{$end};
					 print BED "$count\t$chr\t$start\t$end\t$info\n";
					 $count ++;
				 }
			 }
		 }
		 close BED;
	 }
}

sub write_region_bed {
    my $original_data_dir = shift @_;
	my $sample_pos        = shift @_;
	foreach my $chr(keys %{$sample_pos}) {
	    my $chr_bed_filter = "$original_data_dir/$chr";
		open BED, ">$chr_bed_filter";
        foreach my $start(sort {$a <=> $b} keys %{$sample_pos->{$chr}}) {
            foreach my $end(sort {$a <=> $b} keys %{$sample_pos->{$chr}{$start}}) {
                print BED $sample_pos->{$chr}{$start}{$end};
            }
        }
        close BED;
	}
}

sub region_filter {
    my $correct_dir            = shift @_;
	my $original_data_dir      = shift @_;
	my $region_target_num      = shift @_;
	my $case_suffix            = shift @_;
	my $hash_all_sample_region = shift @_;
	my $sample_pos             = shift @_;
	my $case_suffix_dir        = "$correct_dir/$case_suffix";
	my $analysis_correct       = "$case_suffix_dir/$case_suffix\_result.txt";
	my $correct_filter         = "$case_suffix_dir/$case_suffix\_result_filter.txt";
	open CORRECT, $analysis_correct;
	open FILTER, ">$correct_filter";
	while(<CORRECT>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($chr, $start, $end, $cnv, $type, $log2, $call_prob) = split /\t/, $_;
		my $target_num = $region_target_num->{$chr}{$start}{$end}{'target_num'};
        next if($target_num < 3); # # 过滤：可信区域要覆盖3个及以上target区域（用于绘图和突变注释）
		print FILTER "$chr\t$start\t$end\t$cnv\t$type\n";
		$hash_all_sample_region->{$type}{$chr}{$start}{$end} .= "$case_suffix|CNV_$type:$chr-$start-$end,";
		$sample_pos->{$chr}{$start}{$end}                     = "$chr\t$start\t$end\t$cnv\t$type\n";
	}
	close CORRECT;
	close FILTER;
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
		# # 各结果区域包含target数目统计
		$region_target_num{$region_chr}{$region_start}{$region_end}{'target_num'} ++;
	}
	close BED;
	return %region_target_num;
}


sub generate_tag {  
    my $hashConfig     = shift @_;
	my $hashPara       = shift @_;
    my $report_dir     = shift @_;
    my $output_dir     = shift @_;
	my $result_dir     = shift @_;
    my $status_dir     = shift @_;
    my $hash_group     = shift @_;
	my $region_dict    = shift @_;
	my $region_orifile = shift @_;
    my $genome         = shift @_;
	my $bedtools       = $hashPara->{'Soft'}{'bedtools'};
	my $java           = $hashPara->{"Soft"}{'Java'};
    my $picard         = $hashPara->{"Soft"}{'PicardDIR'};
    my $Tmp            = $hashPara->{"Soft"}{'Tmp'};
	my $region_file    = "$result_dir/merge_".(split /\//, $region_orifile)[-1];
	my $region_file_tp = "$result_dir/tmp_".(split /\//, $region_file)[-1];
	my $target_bed     = "$result_dir/target_".(split /\//, $region_file_tp)[-1];
    my %all_cases;
	print "\nProcessing generate_tag ...\n";
	
	# # 生成target文件
	system "$bedtools merge -i $region_file\_sort -c 4,5 -o distinct,distinct -d 1 > $region_file_tp"; # # picard的target文件需要5列，文件较大，用后删除
	system "cat $region_dict $region_file_tp > $target_bed";
	
    foreach my $suffix(keys %{$hash_group}) {	
		map{$all_cases{$_}++} package::main::get_sample($hash_group->{$suffix}, 'case');
        map{$all_cases{$_}++} package::main::get_sample($hash_group->{$suffix}, 'control');		
    }
    my @config_samples = keys %all_cases;
    my @samples        = check_tag_file($report_dir, $status_dir, $output_dir, \@config_samples);
    
    if(@samples > 0) { # # 运行
        my $threshold = exists $hashPara->{"Process"}{"tag"} ? $hashPara->{"Process"}{"tag"} : 20;
		   $threshold = $hashConfig->{"Process_tag"} if(exists $hashConfig->{"Process_tag"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::main::process_bar_array($sample, \@samples)});# 进度条
        foreach my $sample(@samples) {
            $pm->start($sample) and next;
			my $bam_final   = "$output_dir/$sample/$sample\_final.bam";
            my $tag         = "$status_dir/$sample\_final.tag";
            my $metrics     = "$status_dir/$sample\_final.metrics";
	        system "cpulimit -l 1000 $java -jar -Xmx3g $picard CollectHsMetrics INPUT=$bam_final OUTPUT=$metrics R=$genome BI=$target_bed TI=$target_bed PER_TARGET_COVERAGE=$tag MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 VALIDATION_STRINGENCY=LENIENT TMP_DIR=$Tmp";
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else {
        print "[Note] Tag Process None, for reProcess Delete Result\n";
    }
	# # 删除区域文件（文件较大，且后续用不到）
	system "rm $region_file_tp" if(-e $region_file_tp);
	system "rm $target_bed" if(-e $target_bed);
}

sub check_tag_file {
    my $report_dir     = shift @_;
    my $status_dir     = shift @_;
	my $output_dir     = shift @_;
	my $config_samples = shift @_;
	my %hashCondition; 
    foreach my $sample(@{$config_samples})
    {
   	 	my $final_bam  = "$output_dir/$sample/$sample\_final.bam"; 
    	my $final_file = "$status_dir/$sample\_final.tag";
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
    package::main::write_log("$report_dir/run.log", "EXCAVATORTag", \%hashCondition); 
	my @samples = keys %{$hashCondition{"Good2Run"}};
	return @samples;
}

1
