package package::main;
use strict;
use warnings;

# 获取默认参数
sub get_para{
    my %hashPara;

    # 并行线程数配置
    $hashPara{'Process'}{'excavator2'}            = 6;
    $hashPara{'Process'}{'excavator2_plot'}       = 10;
	$hashPara{'Process'}{'split2target'}          = 10;
	$hashPara{'Process'}{'single_sample_process'} = 5;
	
	# 物种配置
	$hashPara{'Human_hg19'}{'Genome'}       = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.fa";
	$hashPara{'Human_hg19'}{'Dict'}         = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.dict";
    $hashPara{'Human_hg19'}{'AnnovarBuild'} = "hg19";
	
	$hashPara{'Human_hg38'}{'Genome'}       = "/home/genesky/database/ucsc/hg38/genome/hg38.fa";
	$hashPara{'Human_hg38'}{'Dict'}         = "/home/genesky/database/ucsc/hg38/genome/hg38.dict";
	$hashPara{'Human_hg38'}{'AnnovarBuild'} = "hg38";
    
	# CNV数据库配置
	$hashPara{'V6r2'}{'CNV_DB'}     = "/home/genesky/database/self_build_database/cnv/Exon300_CNVdb.txt";
	
	# 软件配置
	$hashPara{"Soft"}{'Tmp'}        = "/home/tmp";
    $hashPara{"Soft"}{'PicardDIR'}  = "/home/genesky/software/picard/2.9.2/picard.jar";
    $hashPara{"Soft"}{'Java'}       = "/home/genesky/software/java/1.8.0_181/bin/java";
	$hashPara{'Soft'}{'EXCAVATOR2'} = "/home/genesky/software/excavator2/1.1.2";
	$hashPara{'Soft'}{'AnnovarDIR'} = "/home/genesky/software/annovar/2018Apr16";
	$hashPara{'Soft'}{'bedtools'}   = "/home/genesky/software/bedtools/2.28.0/bin/bedtools";
	$hashPara{'Soft'}{'R'}          = "/home/genesky/software/r/3.5.1/bin/R";
	$hashPara{'Soft'}{'R_dir'}      = "/home/genesky/software/r/3.5.1/bin";
	$hashPara{'Soft'}{'RScript'}    = "/home/genesky/software/r/3.5.1/bin/Rscript";
	$hashPara{'Soft'}{'RLib'}       = "/home/genesky/software/r/3.5.1/lib64/R/library";

    return %hashPara;
}

# 读取配置文件
sub read_config{
	my $config_file   = shift @_;
	my %hashConfig    = ();
	my %sample_record = ();
	open FILE, $config_file or die "[ERR] Loss Config File, $config_file\n";
	while(my $line = <FILE>) {
		$line =~ s/[\s\t\r\n]//g;
		my @split_line = split /\=/, $line;
		next if(@split_line !=2 and @split_line != 4);
		my $group_num = 1;
		if($split_line[0] eq 'Group') { # 有些control样本相同的case被分开配置，在此处合并
			my @split_sample = split /;/,$split_line[1];
			my @cases        = split /,/, $split_sample[0];
			my @controls     = split /,/, $split_sample[1];
			my $control_list = join ",", @controls;
			$sample_record{$group_num}{$control_list}{'case'}    .= (join ",", @cases).",";
		    $sample_record{$group_num}{$control_list}{'control'} = (join ",", @controls).",";
			$group_num ++;
		}
		elsif($split_line[0] eq 'Sample') {
		    $hashConfig{'Sample'}{$split_line[1]}{$split_line[2]} = $split_line[3];
		}else{
			$hashConfig{$split_line[0]} = $split_line[1];
		}
	}
	close FILE;
	save_sample(\%hashConfig, \%sample_record);
	return %hashConfig;
}

sub read_cnvDB {
    my $hashPara   = shift @_;
	my $DB_version = shift @_;
	my %hash_cnvDB;
	my $cnvDB      = $hashPara->{$DB_version}{'CNV_DB'};
	open DB, $cnvDB;
	while(<DB>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($chr, $start, $end, $gain_freq, $loss_freq, $gain_sample, $loss_sample) = split /\t/, $_;
		$hash_cnvDB{$chr}{$start}{$end}{'cnvDB_Gain_Freq'}   = $gain_freq;
		$hash_cnvDB{$chr}{$start}{$end}{'cnvDB_Loss_Freq'}   = $loss_freq;
		$hash_cnvDB{$chr}{$start}{$end}{'cnvDB_Gain_Sample'} = $gain_sample;
		$hash_cnvDB{$chr}{$start}{$end}{'cnvDB_Loss_Sample'} = $loss_sample;		
	}
	close DB;
	return %hash_cnvDB;
}

# 按control样本数目标注后缀
sub save_sample {
    my $hashConfig    = shift @_;
	my $sample_record = shift @_;
	my $suffix        = "";
	my $group_count   = 1;
	foreach my $group_num(sort {$a <=> $b} keys %{$sample_record}) { # # $group_num 用于保证分组编号顺序与配置文件顺序一致
	    foreach my $control_list(keys %{$sample_record->{$group_num}}) {
	    	my @controls  = get_sample($sample_record->{$group_num}{$control_list}, 'control');
	    	my $control1 = $controls[0];
	    	$suffix = "VS.control$group_count" if(@controls > 1);
	    	$suffix = "VS.$control1" if(@controls <= 1);
	    	$hashConfig->{'Group'}{$suffix}{'case'}    .= $sample_record->{$group_num}{$control_list}{'case'};
	        $hashConfig->{'Group'}{$suffix}{'control'} .= $sample_record->{$group_num}{$control_list}{'control'};
	    	$group_count ++;
	    }
	}
}

# 从屏幕获得输入
sub get_input{
    my $input = <STDIN>;
    $input =~ s/[\r\n]//g;
    return $input;
}

# 获取时间
sub get_time{
    my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
    $year+=1900;
    $mon+=1;
    return "$year-$mon-$day $hour:$min:$sec";
}

# 创建目录
sub make_dir{
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 获取当前项目分析的样本
sub get_all_sample{
    my $hashConfig = shift @_;
    my @need_types = @_;
    my @all_samples;
    
	foreach my $suffix(sort {$a cmp $b} keys %{$hashConfig->{"Group"}})
	{
        foreach my $need_type(@need_types)
		{
            next if(!exists($hashConfig->{"Group"}{$suffix}{$need_type}));
            foreach my $sample(split /,/, $hashConfig->{"Group"}{$suffix}{$need_type})
			{
                push @all_samples, $sample if($sample=~/\w/);
            }
        }
    }
	my %forunique = ();
    my @samples   = grep {!$forunique{$_}++} @all_samples; #除去重复样本
    return @samples;
}

# 获取当前组的样本
sub get_sample{
    my $hashConfig = shift @_;
    my @need_types = @_;
    my @samples;

    foreach my $need_type(@need_types)
    {
        next if(!exists($hashConfig->{$need_type}));
        foreach my $sample(split /,/, $hashConfig->{$need_type})
        {
            push @samples, $sample if($sample=~/\w/);
        }
    }

    return @samples;
}

# 获取当前项目分析的样本组
sub get_all_case_suffix{
    my $hashConfig  = shift @_;
	my $correct_dir = shift @_;
    my @case_suffixes;
    
	foreach my $suffix(sort {$a cmp $b} keys %{$hashConfig->{"Group"}})
	{
        next if(!exists($hashConfig->{"Group"}{$suffix}{'case'}));
        foreach my $case(split /,/, $hashConfig->{"Group"}{$suffix}{'case'})
	    {   
	        my $case_suffix_dir   = "$correct_dir/$case.$suffix";
	        my $analysis_correct = "$case_suffix_dir/$case.$suffix\_result.txt";
			next if(not is_file_ok($analysis_correct)); # # 跳过没有软件结果的样本组
            push @case_suffixes, "$case.$suffix" if($case=~/\w/);
        }
    }
    return @case_suffixes;
}

# 写运行日志
sub write_log{
    my $log_file = shift @_;
    my $title    = shift @_;
    my $hashCondition = shift @_;

    my $tab="##########";
    open LOG,">>$log_file";
    print LOG "$tab Start $title ".get_time()." $tab\n" if($title=~/\w/);
    # Finish
    my @finish_samples = exists $hashCondition->{'Finish'} ? keys %{$hashCondition->{'Finish'}} : ();
    print LOG 'Finish   : '. (join ",", @finish_samples) . "\n"; 
    # Good2Run
    my @good2run_samples = exists $hashCondition->{'Good2Run'} ? keys %{$hashCondition->{'Good2Run'}} : (); 
    print LOG 'Good2Run : '. (join ",", @good2run_samples) . "\n"; 
    # Error
    my @error_samples = exists $hashCondition->{'Error'} ? keys %{$hashCondition->{'Error'}} : ();
    print LOG 'Error    : '. (join ",", @error_samples) . "\n"; 

    print LOG "\n";
    close LOG;
}

# 进度条， 根据输入向量计算
sub process_bar_array{
    my $process_name   = shift @_; # 需要分析的对象
    my $process_arrays = shift @_; # 总列表
    my $process_all_count = @$process_arrays;# 需要分析的总列表数量
    my ($pos) = grep $$process_arrays[$_] eq $process_name, 0..$process_all_count-1;
    my $process_count = $pos+ 1; # 当前对象的位置
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count; # 进度
    my $windows = 100; # 窗口宽度
    my $finished = $windows * int($process_perc) / 100; # 完成
    my $unfinished = $windows - $finished; # 未完成
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\n" if($process_count == $process_all_count);   
}

# 设定表格样式
sub formatRun{
    my ($workbook)=@_;
    my %format=();
    $format{'title'} = $workbook->add_format();
    $format{'title'} ->set_align('center');
    $format{'title'} ->set_align('vcenter');
    $format{'title'} ->set_size(12);
    $format{'title'} ->set_font("Times New Roman");
    $format{'title'} ->set_border();
    $format{'title'} ->set_bg_color("yellow");
    $format{'title'} ->set_color("black");
    
    $format{'normal'} = $workbook->add_format();
    $format{'normal'} ->set_align('center');
    $format{'normal'} ->set_align('vcenter');
    $format{'normal'} ->set_size(12);
    $format{'normal'} ->set_font("Times New Roman");
    $format{'normal'} ->set_border();
    
    $format{'small'} = $workbook->add_format();
    $format{'small'} ->set_align('vcenter');
    $format{'small'} ->set_size(10);
    $format{'small'} ->set_font("Times New Roman");
    $format{'small'} ->set_border();
    
    $format{'seq'} = $workbook->add_format();
    $format{'seq'} ->set_align('vcenter');
    $format{'seq'} ->set_size(11);
    $format{'seq'} ->set_font("Courier New");
    $format{'seq'} ->set_border();
    
    $format{'left'} = $workbook->add_format();
    $format{'left'} ->set_align('vcenter');
    $format{'left'} ->set_size(12);
    $format{'left'} ->set_font("Times New Roman");
    $format{'left'} ->set_border();
    
    $format{'orange'} = $workbook->add_format();
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();

    $format{'skyblue'} = $workbook->add_format();
    $format{'skyblue'} ->set_align('vcenter');
    $format{'skyblue'} ->set_size(12);
    $format{'skyblue'} ->set_font("Times New Roman");
    $format{'skyblue'} ->set_bg_color("#538ed5");
    $format{'skyblue'} ->set_border();

    $format{'bold'} = $workbook->add_format( bold => 1 );
    $format{'blue'} = $workbook->add_format( color => "#538ed5" );
    $format{'redbold'} = $workbook->add_format( color => "#ff0000", bold => 1, );
    $format{'italic'} = $workbook->add_format( italic => 1 );
    $format{'boldblue'} = $workbook->add_format( bold => 1, color => "#538ed5" );
    $format{'bolditalic'} = $workbook->add_format( bold => 1, italic => 1 );
    $format{'blueitalic'} = $workbook->add_format( color => "#538ed5", italic => 1 );
    $format{'boldblueitalic'} = $workbook->add_format( bold => 1, color => "#538ed5", italic => 1 );
    
    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    
    $format{'readme1'} = $workbook->add_format();
    $format{'readme1'}->set_align('center');
    $format{'readme1'}->set_align('vcenter');
    $format{'readme1'}->set_bold();
    $format{'readme1'}->set_size(14);
    $format{'readme1'}->set_font("Times New Roman");
    $format{'readme1'}->set_border();

    $format{'readme2'} = $workbook->add_format();
    $format{'readme2'}->set_align('vcenter');
    $format{'readme2'}->set_bold();
    $format{'readme2'}->set_size(14);
    $format{'readme2'}->set_font("Times New Roman");

    $format{'readme2tmp'} = $workbook->add_format();
    $format{'readme2tmp'}->set_right();

    $format{'readme3'} = $workbook->add_format();
    $format{'readme3'}->set_align('center');
    $format{'readme3'}->set_align('vcenter');
    $format{'readme3'}->set_bold();
    $format{'readme3'}->set_size(11);
    $format{'readme3'}->set_font("Times New Roman");
    $format{'readme3'}->set_border();

    $format{'readme4'} = $workbook->add_format();
    $format{'readme4'}->set_align('vcenter');
    $format{'readme4'}->set_bold();
    $format{'readme4'}->set_size(11);
    $format{'readme4'}->set_font("Times New Roman");
    $format{'readme4'}->set_border();

    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    return %format;
}

1