package package::output;
use Encode;
use strict;
use warnings;

sub run {
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start output ".package::main::get_time()." ##########\n";
    my $isOK = package::check::output($hashPara, $hashConfig); # 运行条件检验
    die "[ERR] Error exists in output run para, please check\n" if($isOK == 0);
	
	my $hash_group     = $hashConfig->{'Group'};
	my $report_dir     = $hashConfig->{'Report'};
	my $species        = $hashConfig->{'Species'};
	my $AnnovarBuild   = $hashPara->{$species}{"AnnovarBuild"};
	my $document_dir   = "$report_dir/document";
	my $cnv_dir        = "$document_dir/3_CNV";
	my $result_dir     = "$report_dir/excavator2";
	my $correct_dir    = "$result_dir/result_final";
	package::main::make_dir($document_dir);
	package::main::make_dir($cnv_dir);
    my @case_suffixes  = package::main::get_all_case_suffix($hashConfig, $correct_dir);	
		
	my $workbook       = Excel::Writer::XLSX->new("$cnv_dir/CNV_EXCAVATOR2.xlsx");
	my %format         = package::main::formatRun($workbook);
	
	# # 过滤的汇总表
	my $combine_dir    = "$correct_dir/COMBINE_REGION";
	my ($anno_filter_region, $hash_filter_region_infos) = read_filter_region_infos($combine_dir, $AnnovarBuild);
	write_filter_sheet($workbook, \%format, $anno_filter_region, $hash_filter_region_infos, \@case_suffixes);
	
	# # 未过滤的汇总表
	my $all_target_dir = "$correct_dir/ALL_TARGETS";
	my ($anno_all_target, $hash_all_target_infos) = read_all_target_infos($all_target_dir, $AnnovarBuild);
	write_cnv_sheet($workbook, \%format, $anno_all_target, $hash_all_target_infos, \@case_suffixes);
	
	# # 样本表
	my ($anno_sample, $hash_sample_infos) = read_sample_infos($correct_dir, $AnnovarBuild, \@case_suffixes);
	write_sample_sheet($workbook, \%format, $anno_sample, $hash_sample_infos, \@case_suffixes);
    
	# # 其他
	write_readme($workbook, \%format, "readme.txt");
	write_control_group($cnv_dir, $hash_group); # # control样本为多个时用group代替，生成文件对group中control样本进行说明

}

sub write_control_group {
    my $cnv_dir      = shift @_;
	my $hash_group   = shift @_;
	my $if_need      = 0;
	open FILE, ">$cnv_dir/control_groups.txt";
	foreach my $suffix(sort {$a cmp $b} keys %{$hash_group}) {
	    if($suffix =~ /control\d/) {
		    $if_need = 1;
			my $control_list = $hash_group->{$suffix}{'control'};
			$suffix  =~ s/^VS\.//g;
			print FILE "$suffix = $control_list\n";
		}
	}
	close FILE;
	system("rm $cnv_dir/control_groups.txt") if($if_need == 0);
}

sub write_readme {
    my $workbook = shift @_;
	my $format   = shift @_;
	my $file     = shift @_;
	open IN, "$file";
	my $readme   = $workbook->add_worksheet("ReadMe");
	my $row      = 0;
	$readme->set_column(0,0,15);
	$readme->set_column(1,1,100);
	while(<IN>) {
	    $_ =~ s/[\r\n]//g;
		my @data = split /\t/,$_;
		if($row == 0) {
		    $readme->write($row, 0, decode("gb2312", $data[0]), $format->{"title"});
	        $readme->write($row, 1, decode("gb2312", $data[1]), $format->{"title"});
		}
		else {
		    $readme->write($row, 0, decode("gb2312", $data[0]), $format->{"normal"});
	        $readme->write($row, 1, decode("gb2312", $data[1]), $format->{"normal"});
        }		
		$row ++;
	}
	close IN;
}

sub write_sample_sheet {
    my $workbook          = shift @_;
	my $format            = shift @_;
	my $anno_sample       = shift @_;
	my $hash_sample_infos = shift @_;
	my $case_suffixes     = shift @_;
	print "##### generating Sample-sheet #####\n";
	my @db      = ('cytoBand', 'dgvFreq', 'iscaPathGainCum', 'iscaPathLossCum', 'iscaLikelyPathogenic', 'iscaPathogenic', 'iscaCuratedPathogenic', 'CNVD', 'DECIPHER');
    my @cnv_db  = ('CNV_DB_Gain_Freq', 'CNV_DB_Loss_Freq');
	my @titles1 = ('Chromosome', 'Start', 'End', 'Size', 'Target_Num', 'Gene', 'Gene Region', 'Exon Cover', @db, @cnv_db, 'log2', 'Posterior Probability');
	foreach my $case_suffix(@$case_suffixes) {
		my @titles = (@titles1, $case_suffix);
		my $CNV    = $workbook->add_worksheet($case_suffix);
        $CNV->set_row(0, 60);
        my $row    = 0;
        foreach my $col(0 .. @titles - 1) {
        	$CNV->write($row, $col, $titles[$col], $format->{"title"});
        }
        $row ++;
        foreach my $chr(sort{$a cmp $b}keys %{$hash_sample_infos->{$case_suffix}}) {
        	foreach my $start(sort {$a <=> $b} keys %{$hash_sample_infos->{$case_suffix}{$chr}}) {
        		foreach my $end(keys %{$hash_sample_infos->{$case_suffix}{$chr}{$start}}) {
        			foreach my $col(0 .. @titles - 1) {
						my $v = " ";
        				$v = $anno_sample->{$case_suffix}{$chr}{$start}{$end}{$titles[$col]} if(exists $anno_sample->{$case_suffix}{$chr}{$start}{$end}{$titles[$col]});
        				$v = $hash_sample_infos->{$case_suffix}{$chr}{$start}{$end}{$titles[$col]} if(exists $hash_sample_infos->{$case_suffix}{$chr}{$start}{$end}{$titles[$col]});
        				my $color = exists $hash_sample_infos->{$case_suffix}{$chr}{$start}{$end}{$titles[$col]."_color"} ? $hash_sample_infos->{$case_suffix}{$chr}{$start}{$end}{$titles[$col]."_color"} : "normal";
        				$CNV->write($row,$col,$v,$format->{$color});
        			}
        			$row++;
        		}
        	}
        }
	}
}

sub write_cnv_sheet {
    my $workbook              = shift @_;
	my $format                = shift @_;
	my $anno_all_target       = shift @_;
	my $hash_all_target_infos = shift @_;
	my $case_suffixes         = shift @_;
	my @db     = ('cytoBand', 'dgvFreq', 'iscaPathGainCum', 'iscaPathLossCum', 'iscaLikelyPathogenic', 'iscaPathogenic', 'iscaCuratedPathogenic', 'CNVD', 'DECIPHER');
	my @cnv_db = ('CNV_DB_Gain_Freq', 'CNV_DB_Loss_Freq');
    my @titles = ('Chromosome', 'Start', 'End', 'Size', 'Gene', 'Gene Region', 'Exon Cover', 'Gain_Freq', 'Loss_Freq', 'Gain_Sample', 'Loss_Sample', @cnv_db, @db, @{$case_suffixes});
    print "##### generating CNV_Original sheet #####\n";
    my $CNV    = $workbook->add_worksheet("CNV Original");
    $CNV->set_row(0, 60);
    my $row    = 0;
    foreach my $col(0 .. @titles - 1) {
    	$CNV->write($row, $col, $titles[$col], $format->{"title"});
    }
    $row ++;
    foreach my $chr(sort{$a cmp $b}keys %{$hash_all_target_infos}) {
    	foreach my $start(sort {$a <=> $b} keys %{$hash_all_target_infos->{$chr}}) {
    		foreach my $end(keys %{$hash_all_target_infos->{$chr}{$start}}) {
    			foreach my $col(0 .. @titles - 1) {
    				my $v = " ";
    				$v = $anno_all_target->{$chr}{$start}{$end}{$titles[$col]} if(exists $anno_all_target->{$chr}{$start}{$end}{$titles[$col]});
    				$v = $hash_all_target_infos->{$chr}{$start}{$end}{$titles[$col]} if(exists $hash_all_target_infos->{$chr}{$start}{$end}{$titles[$col]});
    				my $color = exists $hash_all_target_infos->{$chr}{$start}{$end}{$titles[$col]."_color"} ? $hash_all_target_infos->{$chr}{$start}{$end}{$titles[$col]."_color"} : "normal";
    				$CNV->write($row,$col,$v,$format->{$color});
    			}
    			$row++;
    		}
    	}
    }
}

sub write_filter_sheet {
    my $workbook                 = shift @_;
	my $format                   = shift @_;
	my $anno_filter_region       = shift @_;
	my $hash_filter_region_infos = shift @_;
	my $case_suffixes            = shift @_;
	my @db     = ('cytoBand', 'dgvFreq', 'iscaPathGainCum', 'iscaPathLossCum', 'iscaLikelyPathogenic', 'iscaPathogenic', 'iscaCuratedPathogenic', 'CNVD', 'DECIPHER');
	my @cnv_db = ('CNV_DB_Gain_Freq', 'CNV_DB_Loss_Freq');
    my @titles = ('Chromosome', 'Start', 'End', 'Size', 'Target_Num', 'Gene', 'Gene Region', 'Exon Cover', 'Gain_Freq', 'Loss_Freq', 'Gain_Sample', 'Loss_Sample', @cnv_db, @db, @{$case_suffixes});
    print "##### generating CNV_Filter sheet #####\n";
    my $CNV    = $workbook->add_worksheet("CNV Filter");
    $CNV->set_row(0, 60);
    my $row    = 0;
    foreach my $col(0 .. @titles - 1) {
    	$CNV->write($row, $col, $titles[$col], $format->{"title"});
    }
    $row ++;
    foreach my $chr(sort{$a cmp $b}keys %{$hash_filter_region_infos}) {
    	foreach my $start(sort {$a <=> $b} keys %{$hash_filter_region_infos->{$chr}}) {
    		foreach my $end(keys %{$hash_filter_region_infos->{$chr}{$start}}) {
				foreach my $col(0 .. @titles - 1) {
    				my $v = " ";
    				$v = $anno_filter_region->{$chr}{$start}{$end}{$titles[$col]} if(exists $anno_filter_region->{$chr}{$start}{$end}{$titles[$col]});
    				$v = $hash_filter_region_infos->{$chr}{$start}{$end}{$titles[$col]} if(exists $hash_filter_region_infos->{$chr}{$start}{$end}{$titles[$col]});
    				my $color = exists $hash_filter_region_infos->{$chr}{$start}{$end}{$titles[$col]."_color"} ? $hash_filter_region_infos->{$chr}{$start}{$end}{$titles[$col]."_color"} : "normal";
    				$CNV->write($row,$col,$v,$format->{$color});
    			}
    			$row++;
    		}
    	}
    }
}

sub read_sample_infos {
    my $correct_dir     = shift @_;
	my $AnnovarBuild    = shift @_;
	my $case_suffixes   = shift @_;
	my %anno_sample;
	my %hash_sample_infos;
	foreach my $case_suffix(@$case_suffixes) {
	    my $case_suffix_dir = "$correct_dir/$case_suffix";
		# # 注释信息读取
	    my $library         = "$case_suffix_dir/library/library";
		my %anno_sample_tmp;
	    annotation($library, \%anno_sample_tmp, $AnnovarBuild);
		$anno_sample{$case_suffix} = \%anno_sample_tmp;
	    # # 其他信息读取
	    my $belongs_statistics = "$case_suffix_dir/$case_suffix\_belongs.statistics";
	    open INFO, $belongs_statistics;
	    while(<INFO>) {
	        next if($_ !~ /\w/);
	        $_ =~ s/[\r\n]//g;
	    	my ($chr, $start, $end, $cnv, $type, $log2, $call_prob, $target_num, $cnvDB_Gain_Freq, $cnvDB_Loss_Freq) = split /\t/, $_;
	        next if($cnv == 2); # # 排除掉拷贝数为2的区域
			my $color        = 'normal';
		       $color        = 'orange' if($type eq 'gain');
		       $color        = 'skyblue' if($type eq 'loss');
			$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'Chromosome'}            = $chr;
	    	$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'Start'}                 = $start;
	    	$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'End'}                   = $end;
	    	$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'Size'}                  = $end - $start + 1;
            $hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'log2'}                  = $log2;   
			$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'Target_Num'}            = $target_num;	        		
	    	$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'CNV_DB_Gain_Freq'}      = $cnvDB_Gain_Freq;
	    	$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'CNV_DB_Loss_Freq'}      = $cnvDB_Loss_Freq;
			$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{'Posterior Probability'} = $call_prob;	
            $hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{$case_suffix}            = $cnv;
			$hash_sample_infos{$case_suffix}{$chr}{$start}{$end}{$case_suffix."\_color"}  = $color;
	    }
	    close INFO;
	}
	return (\%anno_sample, \%hash_sample_infos);
}

sub read_all_target_infos {
    my $all_target_dir  = shift @_;
	my $AnnovarBuild    = shift @_;
	# # 注释信息读取
	my $library         = "$all_target_dir/library/library";
	my %anno_all_target;
	annotation($library, \%anno_all_target, $AnnovarBuild);
	# # 其他信息读取
	my $all_target_infos = "$all_target_dir/all_target_infos.txt";
	my %hash_all_target_infos;
	open INFO, $all_target_infos;
	while(<INFO>) {
	    next if($_ !~ /\w/);
	    $_ =~ s/[\r\n]//g;
		my ($chr, $start, $end, $Gain_Freq, $Loss_Freq, $gain_sample, $loss_sample, $cnvDB_Gain_Freq, $cnvDB_Loss_Freq, $cnvDB_Gain_Sample, $cnvDB_Loss_Sample, $sample_infos) = split /\t/, $_;
		$hash_all_target_infos{$chr}{$start}{$end}{'Chromosome'}        = $chr;
		$hash_all_target_infos{$chr}{$start}{$end}{'Start'}             = $start;
		$hash_all_target_infos{$chr}{$start}{$end}{'End'}               = $end;
		$hash_all_target_infos{$chr}{$start}{$end}{'Size'}              = $end - $start + 1;
		$hash_all_target_infos{$chr}{$start}{$end}{'Gain_Freq'}         = $Gain_Freq;
		$hash_all_target_infos{$chr}{$start}{$end}{'Loss_Freq'}         = $Loss_Freq;
		$hash_all_target_infos{$chr}{$start}{$end}{'Gain_Sample'}       = $gain_sample;
		$hash_all_target_infos{$chr}{$start}{$end}{'Loss_Sample'}       = $loss_sample;
		$hash_all_target_infos{$chr}{$start}{$end}{'CNV_DB_Gain_Freq'}  = $cnvDB_Gain_Freq;
		$hash_all_target_infos{$chr}{$start}{$end}{'CNV_DB_Loss_Freq'}  = $cnvDB_Loss_Freq;
		$hash_all_target_infos{$chr}{$start}{$end}{'cnvDB_Gain_Sample'} = $cnvDB_Gain_Sample;
		$hash_all_target_infos{$chr}{$start}{$end}{'cnvDB_Loss_Sample'} = $cnvDB_Loss_Sample;
		my $if_delete_region = get_sample_infos($sample_infos, $hash_all_target_infos{$chr}{$start}{$end});
		delete $hash_all_target_infos{$chr}{$start}{$end} if($if_delete_region == 1); # # 删除区域
	}
	close INFO;
	return (\%anno_all_target, \%hash_all_target_infos);
}

sub read_filter_region_infos {
    my $combine_dir         = shift @_;
	my $AnnovarBuild        = shift @_;
	# # 注释信息读取        
	my $library             = "$combine_dir/library/library";
	my %anno_filter_region;
	annotation($library, \%anno_filter_region, $AnnovarBuild);
	# # 其他信息读取
	my $combined_region_infos = "$combine_dir/combined_region_infos.txt";
	my %hash_filter_region_infos;
	open INFO, $combined_region_infos;
	while(<INFO>) {
	    next if($_ !~ /\w/);
	    $_ =~ s/[\r\n]//g;
		my ($chr, $start, $end, $Gain_Freq, $Loss_Freq, $gain_sample, $loss_sample, $cnvDB_Gain_Freq, $cnvDB_Loss_Freq, $cnvDB_Gain_Sample, $cnvDB_Loss_Sample, $sample_infos, $target_num) = split /\t/, $_;
	    $hash_filter_region_infos{$chr}{$start}{$end}{'Chromosome'}        = $chr;
		$hash_filter_region_infos{$chr}{$start}{$end}{'Start'}             = $start;
		$hash_filter_region_infos{$chr}{$start}{$end}{'End'}               = $end;
		$hash_filter_region_infos{$chr}{$start}{$end}{'Size'}              = $end - $start + 1;
		$hash_filter_region_infos{$chr}{$start}{$end}{'Target_Num'}        = $target_num;
		$hash_filter_region_infos{$chr}{$start}{$end}{'Gain_Freq'}         = $Gain_Freq;
		$hash_filter_region_infos{$chr}{$start}{$end}{'Loss_Freq'}         = $Loss_Freq;
		$hash_filter_region_infos{$chr}{$start}{$end}{'Gain_Sample'}       = $gain_sample;
		$hash_filter_region_infos{$chr}{$start}{$end}{'Loss_Sample'}       = $loss_sample;
		$hash_filter_region_infos{$chr}{$start}{$end}{'CNV_DB_Gain_Freq'}  = $cnvDB_Gain_Freq;
		$hash_filter_region_infos{$chr}{$start}{$end}{'CNV_DB_Loss_Freq'}  = $cnvDB_Loss_Freq;
		$hash_filter_region_infos{$chr}{$start}{$end}{'cnvDB_Gain_Sample'} = $cnvDB_Gain_Sample;
		$hash_filter_region_infos{$chr}{$start}{$end}{'cnvDB_Loss_Sample'} = $cnvDB_Loss_Sample;
		my $if_delete_region = get_sample_infos($sample_infos, $hash_filter_region_infos{$chr}{$start}{$end});
        delete $hash_filter_region_infos{$chr}{$start}{$end} if($if_delete_region == 1); # # 删除区域
	}
	close INFO;
	return (\%anno_filter_region, \%hash_filter_region_infos);
}

sub get_sample_infos {
    my $sample_infos     = shift @_;
	my $hash             = shift @_;
	my @all_sample_infos = split /,/, $sample_infos;
	my $cnv2_sample_num  = 0;
	foreach my $all_sample_info(@all_sample_infos) {
	    my ($sample, $cnv, $type) = split /:/, $all_sample_info;
		if($cnv == 2) {# # 排除掉拷贝数为2的区域
		    $cnv2_sample_num ++;
			next;
		}
		my $color        = 'normal';
		   $color        = 'orange' if($type eq 'gain');
		   $color        = 'skyblue' if($type eq 'loss');
		   $hash->{$sample}          = $cnv;
		   $hash->{"$sample\_color"} = $color;
	}
	my $if_delete_region = ($cnv2_sample_num == @all_sample_infos) ? 1 : 0; # # 若该区域全部样本拷贝数均为0，则删除该区域
	return $if_delete_region;
}

sub annotation {
    my $library      = shift @_;
	my $anno         = shift @_;
	my $AnnovarBuild = shift @_;
	return if(not -e $library);
    open FUNC, "$library.variant_function";
    while(<FUNC>) {
    	$_ =~ s/[\r\n]//g;
    	my ($region, $gene, $chr, $start, $end, $tmp) = split /\t/, $_, 6;
    	$gene =~ s/\(.*\)//g;
    	$anno->{$chr}{$start}{$end}{"Gene Region"} = $region;
    	$anno->{$chr}{$start}{$end}{"Gene"}        = $gene;
    }
    close FUNC;  
	open EXO, "$library.exonic_variant_function";
    while(<EXO>) {
    	$_ =~ s/[\r\n]//g;
		my ($id, $function, $infos, $chr, $start, $end, $tmp) = split /\t/, $_, 7;
    	my $info_list  = $infos;
    	next if($info_list !~ /\w/ or $info_list =~ /UNKNOWN/);
    	my @infos      = split /,/, $info_list;
    	my $exon_cover = "";
    	foreach my $info(@infos) {
    	    my ($gene, $transcript, $exon, $tmp) = split /:/, $info, 4;
    		$exon_cover .= "$gene:$transcript:$exon,";
    	}
    	$anno->{$chr}{$start}{$end}{"Exon Cover"} .=$exon_cover;
    }
    close EXO;
    my @db = ("cytoBand", "dgvFreq", "iscaPathGainCum", "iscaPathLossCum", "iscaLikelyPathogenic", "iscaPathogenic", "iscaCuratedPathogenic", "CNVD", "DECIPHER");
    foreach(@db) {
    	open FILE, "$library.$AnnovarBuild\_$_";
    	while(my $line=<FILE>) {
    		$line =~ s/[\r\n]//g;
			my ($database, $value, $chr, $start, $end, $tmp) =  split /\t/, $line, 6;
    		$value =~ s/^Name\=//;
    		if($_ eq "dgvFreq") {
    			my ($gainAdd, $lossAdd, $allAdd)=(0, 0, 0);
    			foreach my $key(split /,/, $value) {
    				my ($gain, $loss, $all) = split /\;/, $key;
    				$gainAdd += $gain;
    				$lossAdd += $loss;
    				$allAdd  += $all;
    			}
    			$value = 0;
				$value = ($gainAdd + $lossAdd) / $allAdd if($allAdd > 0);
    		}
    		$anno->{$chr}{$start}{$end}{$_} = $value;
    	}
    	close FILE;
    }
}

1