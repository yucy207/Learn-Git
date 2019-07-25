package package::report_methyl;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start Report Methylation Result ".package::utils::get_time()." ##########\n";
    my $isOK = package::check::report_methyl($hashConfig); # 运行条件检验
    die "[ERR] Error exists in Report Methyl run para, please check\n" if($isOK == 0);

    # 路径准备
    my $report_dir       = $hashConfig->{'Report'};
    my $methyl_dir       = "$report_dir/methylation";
    my $txt_dir          = "$report_dir/methylation_result";
    my $document_dir     = "$report_dir/document";
    my $qc_dir           = "$report_dir/document/01_QC";    
    my $methyl_excel_dir = "$report_dir/document/02_Methylation_Quantification";
    my $inter_qc_dir     = "$report_dir/document/Internal_QC";
    package::utils::make_dir($txt_dir, $document_dir, $qc_dir, $methyl_excel_dir, $inter_qc_dir);

    # 软件
    my $Rbin = $hashPara->{'Soft'}{'R'};
    my $RLib = $hashPara->{'Soft'}{'RLib'};

    # 样本
    my @allsamples   = package::utils::get_sample($hashConfig, 'case', 'control', 'Ycontrol'); # 所有样本
    my @samples      = package::utils::get_sample($hashConfig, 'case', 'control');             # 客户样本
    my @pair_samples = package::utils::get_pair_sample($hashConfig);                           # 内参配对样本

    # 检查必需文件
    my $isOK_sample = check_sample_result($hashConfig, \@allsamples);
    die if($isOK_sample == 0);

    # 读取数据
    my %hashPoint      = ();  # 需分析的c位点
    my %hashMethyData  = ();  # 甲基化结果
    my %hashEfficiency = ();  # 盐化效率
    my %hashHaplotype  = ();  # 单倍型结果
    my %hashSampleHapN = ();  # 每个样本中前N种最多的单倍型
    read_methy_base($methyl_dir, \@allsamples, \%hashMethyData, \%hashEfficiency);       # 所有甲基化值
    read_analysis_point("$report_dir/c.point.analysis.txt", \%hashPoint);  # 所有甲基化位点信息
    read_haplotype($methyl_dir, \@samples, \%hashHaplotype);               # 所有单倍型信息
    get_sample_firstn(\%hashHaplotype, \@samples, \%hashSampleHapN, 3);    # 获取每个样本中前N种最多的单倍型

    # 文件
    my $seq_region         = "$hashConfig->{'Report'}/seq.region";
    my $efficiency_file    = "$txt_dir/1.efficiency.txt";    # 1. 盐化处理效率
    my $methyl_site_file   = "$txt_dir/2.methyl.site.txt";   # 2. 单位点甲基化值
    my $methyl_target_file = "$txt_dir/3.methyl.target.txt"; # 3. 片段平均甲基化值
    my $methyl_gene_file   = "$txt_dir/4.methyl.gene.txt";   # 4. 基因平均甲基化值
    my $haplotype_file     = "$txt_dir/5.haplotype.txt";     # 5. 单倍型
    my $methyl_info_file   = "$txt_dir/6.methyl.info.txt";   # 6. 单位点甲基化值详细信息
    my $pair_qc_file       = "$txt_dir/7.pair.qc.txt";       # 7. 内部质控结果

    # 生成所需 txt 文件
    generate_efficiency_file(\@allsamples, \%hashEfficiency, $efficiency_file);
    generate_methyl_table_file($hashConfig, \@allsamples, \%hashMethyData, \%hashPoint, $methyl_site_file, $methyl_target_file, $methyl_gene_file);
    generate_haplotype_file(\@allsamples, \%hashHaplotype, \%hashSampleHapN, $haplotype_file);
    generate_methyl_info_file(\@allsamples, \%hashMethyData, \%hashPoint, $methyl_info_file);
    generate_pair_qc_file($txt_dir, $hashConfig->{'pairSamples'}, $pair_qc_file, $Rbin, $RLib) if(scalar(@pair_samples) > 0);

    my $is_file_ok = package::utils::is_file_ok($seq_region, $efficiency_file, $methyl_site_file, $haplotype_file, $methyl_info_file);
    die "Some file lost in $txt_dir\n" if($is_file_ok == 0);

    plot_haplotype_dot($haplotype_file, $methyl_excel_dir, $Rbin, $RLib) if(package::utils::is_file_ok($haplotype_file.".4plot"));

    # 输出excel表格
    my $excel            = "$methyl_excel_dir/Methylation.xlsx";      # 提供给客户
    my $excel_complete   = "$inter_qc_dir/Methylation_complete.xlsx"; # 生信内部阅读,包含Ycontrol样本
    my $excel_pair_qc    = "$inter_qc_dir/QC_Correlation.xlsx";       # 配对样本相关性检验
    my $excel_primer     = "$qc_dir/Target_Primer_infomation.xlsx";   # 片段、引物信息
    my $excel_efficiency = "$qc_dir/Efficiency.xlsx";                 # 盐化效率

    print "output $excel\n";
    print "output $excel_complete\n" if(scalar(@pair_samples) > 0); # 只有存在内部质控样本时，该文件才输出
    print "output $excel_pair_qc\n"  if(scalar(@pair_samples) > 0);
    print "output $excel_primer\n";
    print "output $excel_efficiency\n";
    print "output process : \n";
    my $workbook            = Excel::Writer::XLSX->new($excel);
    my $workbook_complete   = Excel::Writer::XLSX->new($excel_complete) if(scalar(@pair_samples) > 0);
    my $workbook_pair_qc    = Excel::Writer::XLSX->new($excel_pair_qc)  if(scalar(@pair_samples) > 0);
    my $workbook_primer     = Excel::Writer::XLSX->new($excel_primer);
    my $workbook_efficiency = Excel::Writer::XLSX->new($excel_efficiency);
    my %format              = package::utils::sheet_format($workbook);  
    my %format_complete     = package::utils::sheet_format($workbook_complete) if(scalar(@pair_samples) > 0);
    my %format_pair_qc      = package::utils::sheet_format($workbook_pair_qc)  if(scalar(@pair_samples) > 0); 
    my %format_primer       = package::utils::sheet_format($workbook_primer);
    my %format_efficiency   = package::utils::sheet_format($workbook_efficiency);

    # (1) 片段和引物信息输出
    print "    1 - target primer info\n";
    excel_target_primer($workbook_primer, \%format_primer, 'Target.Info', $seq_region, \%hashPoint);

    # (2) 盐化效率信息输出
    print "    2 - efficiency\n";
    excel_efficiency_methyl_info($workbook_efficiency, \%format_efficiency, 'Efficiency', $efficiency_file, \@samples);
    excel_efficiency_methyl_info($workbook_complete, \%format_complete, 'Efficiency', $efficiency_file, \@allsamples) if(scalar(@pair_samples) > 0);
    
    #（3）单位点甲基化值输出
    print "    3 - site methyl\n";
    excel_methyl_table($workbook, \%format, 'Site', $methyl_site_file, \@samples);
    excel_methyl_table($workbook_complete, \%format_complete, 'site', $methyl_site_file, \@allsamples) if(scalar(@pair_samples) > 0);
    
    #（4）片段甲基化值输出
    print "    4 - target methyl\n";
    excel_methyl_table($workbook, \%format, 'Target', $methyl_target_file, \@samples);
    excel_methyl_table($workbook_complete, \%format_complete, 'Target', $methyl_target_file, \@allsamples) if(scalar(@pair_samples) > 0);

    #（5）基因甲基化值输出
    if(exists $hashConfig->{'Gene'}){
        print "    5 - gene methyl\n";
        excel_methyl_table($workbook, \%format, 'Gene', $methyl_gene_file, \@samples);
        excel_methyl_table($workbook_complete, \%format_complete, 'Gene', $methyl_gene_file, \@allsamples) if(scalar(@pair_samples) > 0);
    }
    
    #（6）单倍型信息输出
    print "    6 - haplotype\n";
    excel_haplotype($workbook, \%format, 'Haplotype', $haplotype_file, \@samples);
    excel_haplotype($workbook_complete, \%format_complete, 'Haplotype', $haplotype_file, \@samples) if(scalar(@pair_samples) > 0);

    #（7）详细甲基化信息输出
    print "    7 - methyl detail\n";
    excel_efficiency_methyl_info($workbook, \%format, 'Methyl.detail', $methyl_info_file, \@samples);
    excel_efficiency_methyl_info($workbook_complete, \%format_complete, 'Methyl.detail', $methyl_info_file, \@allsamples) if(scalar(@pair_samples) > 0);

    # (8) 内部对照相关性检验
    print "    8 - pair QC\n";
    excel_pair_qc_methyl($workbook_pair_qc, \%format_pair_qc, 'Raw.methyl.data', \%hashMethyData, \%hashPoint, \@pair_samples) if(scalar(@pair_samples) > 0);
    excel_pair_qc_reads($workbook_pair_qc, \%format_pair_qc, 'Raw.reads.data', $hashConfig, $hashConfig->{'pairSamples'}) if(scalar(@pair_samples) > 0);
    excel_pair_qc_cor($workbook_pair_qc, \%format_pair_qc, 'Pearson.correlation.result', $pair_qc_file) if(scalar(@pair_samples) > 0);

    # (9) readme
    print "    9 - readme\n";
    readme($workbook,            \%format,            'readme/readme.txt');
    readme($workbook_complete,   \%format_complete,   'readme/readme.txt') if(scalar(@pair_samples) > 0);
    readme($workbook_pair_qc,    \%format_pair_qc,    'readme/readme_correlation.txt') if(scalar(@pair_samples) > 0); 
    readme($workbook_primer,     \%format_primer,     'readme/readme_primer.txt');
    readme($workbook_efficiency, \%format_efficiency, 'readme/readme_efficiency.txt');

    print "    *Finished* \n";

}

sub plot_haplotype_dot{
    my $haplotype_file   = shift @_;
    my $methyl_excel_dir = shift @_;
    my $Rbin             = shift @_;
    my $RLib             = shift @_;
    my $hap4plot_file    = $haplotype_file.".4plot";
    my %hashHapPlot;
    my $max_length = 0; 
    open PLOT,$hap4plot_file;
    while(<PLOT>){
        $_=~s/[\r\n]//g;
        my($target, $haps, $support_reads) = split/\t/,$_;
        my @haps = split/,/,$haps;
        my @support_reads = split/,/,$support_reads;
        $hashHapPlot{$target}{'hap'} = \@haps;
        $hashHapPlot{$target}{'reads'} = \@support_reads;
        my $length = length($haps[0]);
        $max_length = $length if($length > $max_length);
    }
    close PLOT;
    my $R = Statistics::R->new(bin => $Rbin);
    $R->startR;
    $R->send(qq` .libPaths("$RLib") \n`) ;
    # $R->send(qq` pdf("$methyl_excel_dir/haplotype.dot.pdf", width=0.25*$max_length+2, height=3) \n`) ;
    $R->send(qq` pdf("$methyl_excel_dir/haplotype.dot.pdf", width=0.25*$max_length+6, height=3) \n`) ;
    $R->send(qq` par(mar = c(3, 2, 2, 2)) \n`) ;
    foreach my $target(sort keys %hashHapPlot){
        $R->send(qq` plot(c(1:$max_length), xlim=c(0,$max_length+2), ylim=c(0,5), type = "n", xaxt="n", yaxt="n", bty="n", ylab="", xlab="",main="") \n`) ;
        my $y = 3;
        my @haps = @{$hashHapPlot{$target}{'hap'}};
        my @reads = @{$hashHapPlot{$target}{'reads'}};
        foreach my $i(0..$#haps){
            my $length = length($haps[$i]);
            my @words  = split//, $haps[$i];
            # $R->send(qq` abline(h = $y, lwd = 1.5) \n`) ;
            my $xmax = $length + 1;
            my $yrep = $length + 2;
            $R->send(qq` lines(x = c(0:$xmax),y = rep($y,$yrep), lwd = 1.5) \n`) ;            
            foreach my $x(1..@words){
                $R->send(qq` points(x = $x, y = $y, cex = 2.8, bg = "white", pch=21) \n`) if($words[$x-1] eq "t");
                $R->send(qq` points(x = $x, y = $y, cex = 2.8, bg = "red", pch=21) \n`) if($words[$x-1] eq "c");
            }
            $R->send(qq` text(x = $xmax+0.2, y = $y, labels="$reads[$i]", adj=0) \n`);
            $y--; 
        }
        $R->send(qq` text(x = 1+$max_length/2, y = 4, labels="$target", cex=1.8) \n`);
    }
    $R->send(qq` dev.off() \n`) ;
}

sub excel_pair_qc_cor{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;

    my %hashData = read_file($input_file);

    my $sheet = $workbook->add_worksheet($sheet_name);

    # 表头
    my @heads = split/\t/, $hashData{'HEAD'};
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {
        my @datas;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            push @datas, $data;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;
    }

}

sub excel_pair_qc_reads{
    my $workbook     = shift @_;
    my $format       = shift @_;
    my $sheet_name   = shift @_;
    my $hashConfig   = shift @_;
    my $pair_info    = shift @_;

    my $status_dir = "$hashConfig->{'Report'}/status";
    my %hashReads;
    my @uniq_pair_samples; # 防止重复样本重复输出
    foreach my $pairSample(split /;/, $pair_info){
        foreach my $sample(split/,/,$pairSample){
            push @uniq_pair_samples, $sample;
            my $status_file = "$status_dir/$sample.status";
            open STATUS, $status_file;
            while(<STATUS>){
                $_=~s/[\r\n]//;
                next if($_=~/^#/);
                my($target, $reads) = split/\t/, $_;
                $hashReads{$target}{$sample}  = $reads;
                $hashReads{$target}{'Target'} = $target;
            }
            close STATUS;
        }       
    }
    
    my $sheet = $workbook->add_worksheet($sheet_name);
    # 表头
    my %count;
    @uniq_pair_samples = grep { ++$count{ $_ } < 2; } @uniq_pair_samples;
    my @heads = ("Target", @uniq_pair_samples);
    my $row   = 0;
    foreach my $col(0..$#heads){
        $sheet->write($row, $col, $heads[$col], $format->{"title"});
    }
    $row++;
    foreach my $target(sort keys %hashReads){
        foreach my $col(0..$#heads){
            my $value = "";
               $value = $hashReads{$target}{$heads[$col]} if(exists $hashReads{$target}{$heads[$col]});
            $sheet->write($row, $col, $value, $format->{"normal"});
        }
        $row++;
    }
}

sub excel_pair_qc_methyl{
    my $workbook      = shift @_;
    my $format        = shift @_;
    my $sheet_name    = shift @_;
    my $hashMethyData = shift @_;
    my $hashPoint     = shift @_;
    my $pair_samples  = shift @_;
    
    my $sheet = $workbook->add_worksheet($sheet_name);

    # 表头
    my @heads = ("Target","Pos",@$pair_samples);
    my $row   = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $target(sort keys %$hashPoint){
        foreach my $position(sort { $a <=> $b } keys %{$hashPoint->{$target}}){
            $sheet->write($row, 0, $target, $format->{"normal"});
            $sheet->write($row, 1, $position, $format->{"normal"});
            foreach my $col(2..$#heads){
                my $value = "";
                   $value = (exists($hashMethyData->{$target}{$position}{$heads[$col]}{'Allcount'}) and $hashMethyData->{$target}{$position}{$heads[$col]}{'Allcount'}>=10) ? $hashMethyData->{$target}{$position}{$heads[$col]}{'Methycount'}/$hashMethyData->{$target}{$position}{$heads[$col]}{'Allcount'}:"";
                $sheet->write($row, $col, $value, $format->{"normal"});
            }
            $row++;
        }
    }
}

sub excel_haplotype{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;
    my $samples    = shift @_;

    my %hashData = read_file($input_file);

    my $sheet = $workbook->add_worksheet($sheet_name);
       $sheet->set_row(0, 60);
       $sheet->set_column(0,1,15);

    # 表头
    my @heads = ("Target","Haplotype","Depth",@$samples);
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    my %hashFlag = ();
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {
        my $target = $hashData{'Value'}{$row_count}{'Target'};
        $hashFlag{$target}++; # 只标记每个片段第一个单倍型（深度最高）

        my @datas;
        my @colors;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            my $color = "normal";
               $color = "orange" if($data =~ /\w/ and $head eq 'Haplotype' and $hashFlag{$target}==1 );
            push @datas, $data;
            push @colors, $color;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{$colors[$col]});
        }
        $row++;
    }
}

sub excel_methyl_table{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;
    my $samples    = shift @_;

    my %hashData = read_file($input_file);

    my $sheet = $workbook->add_worksheet($sheet_name);
       $sheet->set_row(0, 60);
       $sheet->set_column(0,0,15);
       $sheet->set_column(3,4,13) if($sheet_name eq 'Site');

    # 表头
    my @heads = ("Target","Position","Chr","GenomePosition","Distance2TSS","Type",@$samples);
       @heads = ("Target",@$samples) if($sheet_name eq 'Target');
       @heads = ("Gene",@$samples) if($sheet_name eq 'Gene');
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {
        my @datas;
        my @colors;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            my $color = "normal";
               $color = "orange" if($data =~ /\d/ and $head~~@$samples and $data > 0.5 );
            push @datas, $data;
            push @colors, $color;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{$colors[$col]});
        }
        $row++;
    }
}

sub excel_efficiency_methyl_info{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;
    my $samples    = shift @_;

    my %hashData = read_file($input_file);

    my $sheet = $workbook->add_worksheet($sheet_name);
       $sheet->set_row(0, 60);
       $sheet->set_column(0,3,15) if($sheet_name eq "Efficiency");
       $sheet->set_column(0,0,15) if($sheet_name eq "Methyl.detail");

    # 表头
    my @heads = split/\t/, $hashData{'HEAD'};
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {   
        my $sample = $hashData{'Value'}{$row_count}{'Sample'};
        next if(!($sample~~@$samples));

        my @datas;
        my @colors;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            my $color = "normal";
               $color = "orange" if($data =~ /\d/ and (($head eq 'Efficiency' and $data < 80) or ($head eq 'Percentage' and $data > 0.5)) );
            $data = $data."%" if($head eq 'Efficiency');
            push @datas, $data;
            push @colors, $color;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{$colors[$col]});
        }
        $row++;
    }
}

sub excel_target_primer{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $seq_region = shift @_;
    my $hashPoint  = shift @_;
    my $sheet = $workbook->add_worksheet($sheet_name);
       $sheet->set_row(0, 60);
    my %hashRegion = package::utils::read_seq_region($seq_region);
    my @titles=("Target","Chr","Gene","mRNA","mRNA_Strand","TSS","TES","Start","End","Length","Target_Strand","Distance2TSS","PrimerF","PrimerR","TargetSeq");
    $sheet->set_column(0,$#titles-3,10);  
    $sheet->set_column($#titles-2,$#titles-1,30);  
    $sheet->set_column($#titles,$#titles,200);  
    my $row=0;
    for my $i(0..$#titles){
        $sheet->write($row,$i,$titles[$i],$format->{"title"});
    }
    $row++;
    foreach my $target(sort keys %hashRegion){
        my @CGformats = mark_cg($hashRegion{$target}{'TargetSeq'}, $hashPoint->{$target}, $format);# 高亮产物区域内分析的C
        foreach my $col(0..$#titles){
            my $color = "normal";
            my $value=(exists($hashRegion{$target}{$titles[$col]})) ? $hashRegion{$target}{$titles[$col]}:"";
            $color = "left" if($titles[$col]=~/^Primer/);
            $sheet->write($row,$col,$value,$format->{$color}) if($titles[$col] ne "TargetSeq");
            $sheet->write_rich_string($row,$col,@CGformats,$format->{"left"}) if($titles[$col] eq "TargetSeq");
        }
        $row++;     
    }
}

sub generate_pair_qc_file{
    my $txt_dir       = shift @_;
    my $pair_info     = shift @_;
    my $pair_qc_file  = shift @_;
    my $Rbin          = shift @_;
    my $RLib          = shift @_;
    my @title         = ("pairSample","Site_num","t.value","df","p.value","estimate","L95\-U95");
    my $title = join"\t",@title;
    open OUT, ">$pair_qc_file\n";
    print OUT "$title\n";

    my $R = Statistics::R->new(bin => $Rbin);
       $R->startR;
       $R->send( qq` .libPaths("$RLib") ` );
       $R->send( qq` data = read.table("$txt_dir/2.methyl.site.txt", sep = "\t", header = TRUE, check.names = FALSE ) `);
    foreach my $pairSample(split /;/, $pair_info){
        my($sample1,$sample2)=split/,/,$pairSample;
        $R->send( qq` sample1datas = data\$"$sample1" ` );
        $R->send( qq` sample2datas = data\$"$sample2" ` );
        $R->send( qq` pairdatas = cbind(sample1datas,sample2datas) ` );
        $R->send( qq` pairdatas = pairdatas[complete.cases(pairdatas),] ` );
        $R->send( qq` site_num = nrow(pairdatas) ` );
        
        my($t_value, $df, $p_value, $estimate, $L95, $U95) = map{'NA'}1..6;
        my $L95_U95  = "$L95\-$U95";
        my $site_num = $R->get('site_num');
        next if($site_num<=1);
        $R->send( qq` result = cor.test(pairdatas[,1], pairdatas[,2], method = "pearson", conf.level = 0.95) ` );
        $R->send( qq` t.value = result\$statistic[[1]] ` );
        $R->send( qq` df = result\$parameter[[1]] ` );
        $R->send( qq` p.value = result\$p.value[[1]] ` );
        $R->send( qq` cor = result\$estimate[[1]] ` );
        $R->send( qq` if(site_num>=4) l95 = result\$conf.int[[1]]` );
        $R->send( qq` if(site_num>=4) u95 = result\$conf.int[[2]] ` );
        $t_value  = $R->get('t.value');
        $df       = $R->get('df');
        $p_value  = $R->get('p.value');
        $estimate = $R->get('cor');
        $L95      = $R->get('l95');
        $U95      = $R->get('u95');
        $L95      = sprintf "%0.4f", $L95 if($L95 =~/\d/);
        $U95      = sprintf "%0.4f", $U95 if($U95 =~/\d/);
        $L95_U95  = "$L95\-$U95";
        print OUT "$pairSample\t$site_num\t$t_value\t$df\t$p_value\t$estimate\t$L95_U95\n";
    }
    close OUT;
}

sub generate_methyl_info_file{
    my $allsamples       = shift @_;
    my $hashMethyData    = shift @_;
    my $hashPoint        = shift @_;
    my $methyl_info_file = shift @_;
    my @title = ("Sample","Target","Position","Chr","GenomePosition","Distance2TSS","Type","Mathy. Base(C)","All Base","Percentage");
    my $title = join"\t",@title;
    open OUT, ">$methyl_info_file\n";
    print OUT "$title\n";
    foreach my $sample(@$allsamples){
        foreach my $target(sort keys %$hashPoint){
            foreach my $position(sort { $a <=> $b }keys %{$hashPoint->{$target}}){
                my @datas;
                push @datas, $sample;
                push @datas, $target;
                push @datas, $position;
                push @datas, $hashPoint->{$target}{$position}{'Chr'};
                push @datas, $hashPoint->{$target}{$position}{'GenomePosition'};
                push @datas, $hashPoint->{$target}{$position}{'Distance2TSS'};
                push @datas, $hashPoint->{$target}{$position}{'Type'};
                my $methyc   = (exists ($hashMethyData->{$target}{$position}{$sample}{'Methycount'}) ) ? $hashMethyData->{$target}{$position}{$sample}{'Methycount'} : "";
                my $Allcount = (exists ($hashMethyData->{$target}{$position}{$sample}{'Allcount'}) ) ? $hashMethyData->{$target}{$position}{$sample}{'Allcount'} : "";
                my $pct      = ($methyc =~ /\d/ and $Allcount =~ /\d/ and $Allcount > 0) ? $methyc/$Allcount : "";
                push @datas, $methyc;
                push @datas, $Allcount;
                push @datas, $pct;
                my $datas = join"\t", @datas;
                print OUT "$datas\n";       
            }
        }
    }
    close OUT;
}

sub generate_haplotype_file{
    my $allsamples     = shift @_;
    my $hashHaplotype  = shift @_;
    my $hashSampleHapN = shift @_;
    my $haplotype_file = shift @_;
    my $hap4plot_file  = $haplotype_file.".4plot";   #用于绘图的单倍型文件
    my @title = ("Target","Haplotype","Depth",@$allsamples);
    my $title = join"\t",@title;
    open OUT, ">$haplotype_file\n";
    open PLOT, ">$hap4plot_file\n";
    print OUT "$title\n";
    foreach my $target(sort keys %$hashSampleHapN){
        my @haplotypes = split /,/, $hashSampleHapN->{$target};
        my @support_reads ;
        foreach my $haplotype(@haplotypes){
            my @datas;
            push @datas, $target;
            push @datas, $haplotype;
            push @datas, $hashHaplotype->{'hap_sum'}{$target}{$haplotype};
            push @support_reads, $hashHaplotype->{'hap_sum'}{$target}{$haplotype};
            foreach my $sample(@$allsamples){
                my $reads = (exists ($hashHaplotype->{'hap1'}{$target}{$haplotype}{$sample})) ? $hashHaplotype->{'hap1'}{$target}{$haplotype}{$sample}: 0;
                my $all   = (exists ($hashHaplotype->{'sample_sum'}{$target}{$sample})) ? $hashHaplotype->{'sample_sum'}{$target}{$sample} : "";
                my $pct   = ($all =~ /\d/ and $all > 0) ? $reads/$all : "";
                push @datas, $pct;
            }
            my $datas = join"\t", @datas;
            print OUT "$datas\n";
        }
        
        my $haplotype_s    = "";
        my $support_reads_s = "";
	foreach my $count(0..$#haplotypes)
        {
            $haplotype_s .= "$haplotypes[$count],";
            $support_reads_s .= "$support_reads[$count],";
            last if($count >=2); # 仅保留前3个
        }
        $haplotype_s =~ s/,$//;
        $support_reads_s =~ s/,$//;
        print PLOT "$target\t$haplotype_s\t$support_reads_s\n";
    }
    close OUT;
    close PLOT;
}

sub generate_methyl_table_file{
    my $hashConfig         = shift @_;
    my $allsamples         = shift @_;
    my $hashMethyData      = shift @_;
    my $hashPoint          = shift @_;
    my $methyl_site_file   = shift @_;
    my $methyl_target_file = shift @_;
    my $methyl_gene_file   = shift @_;
    my $MeanSNPStat        = (exists($hashConfig->{'MeanSNPStat'})) ? $hashConfig->{'MeanSNPStat'} : 'keep';
    
    # 单位点甲基化值输出
    my @title_site = ("Target","Position","Chr","GenomePosition","Distance2TSS","Type",@$allsamples);
    my $title_site = join"\t", @title_site;
    open OUT, ">$methyl_site_file\n";
    print OUT "$title_site\n";
    foreach my $target(sort keys %$hashPoint){
        foreach my $position(sort { $a <=> $b } keys %{$hashPoint->{$target}}){
            my @datas;
            push @datas, $target;
            push @datas, $position;
            push @datas, $hashPoint->{$target}{$position}{'Chr'};
            push @datas, $hashPoint->{$target}{$position}{'GenomePosition'};
            push @datas, $hashPoint->{$target}{$position}{'Distance2TSS'};
            push @datas, $hashPoint->{$target}{$position}{'Type'};
            foreach my $sample(@$allsamples){
                my $pct=(exists($hashMethyData->{$target}{$position}{$sample}) and $hashMethyData->{$target}{$position}{$sample}{'Allcount'}>=10) ? $hashMethyData->{$target}{$position}{$sample}{'Methycount'}/$hashMethyData->{$target}{$position}{$sample}{'Allcount'}:"";
                push @datas, $pct;      
            }
            my $datas = join"\t", @datas;
            print OUT "$datas\n";
        }
    }
    close OUT;

    # 片段平均甲基化值输出
    my @title_target = ("Target",@$allsamples);
    my $title_target = join"\t", @title_target;
    open OUT, ">$methyl_target_file\n";
    print OUT "$title_target\n";
    foreach my $target(sort keys %$hashPoint){
        my @datas_target;
        push @datas_target, $target;
        foreach my $sample(@$allsamples){
            my $sum = 0;my $num = 0;
            foreach my $position(sort { $a <=> $b } keys %{$hashPoint->{$target}}){
                next if($MeanSNPStat eq 'remove' and $hashPoint->{$target}{$position}{'Type'}=~/\/rs\d+/); # 排除SNP位点      
                my $pct = (exists($hashMethyData->{$target}{$position}{$sample}) and $hashMethyData->{$target}{$position}{$sample}{'Allcount'}>=10) ? $hashMethyData->{$target}{$position}{$sample}{'Methycount'}/$hashMethyData->{$target}{$position}{$sample}{'Allcount'}:"";
                if($pct =~ /\d+/){$sum+=$pct;$num++;}   
            }
            my $mean = ($num > 0) ? $sum/$num : "";
            push @datas_target, $mean;
        }
        my $datas_target = join"\t", @datas_target;
        print OUT "$datas_target\n";
    }
    close OUT;

    # 基因平均甲基化值输出
    return if(not exists $hashConfig->{'Gene'});
    my @title_gene = ("Gene",@$allsamples);
    my $title_gene = join"\t", @title_gene;
    open OUT, ">$methyl_gene_file\n";
    print OUT "$title_gene\n";
    foreach my $gene(sort keys %{$hashConfig->{'Gene'}}){
        my @datas_gene;
        push @datas_gene, $gene;
        foreach my $sample(@$allsamples){
            my $sum = 0;my $num = 0;
            foreach my $target(sort keys %{$hashConfig->{'Gene'}{$gene}}){
                foreach my $position(sort { $a <=> $b } keys %{$hashPoint->{$target}}){
                    next if($MeanSNPStat eq 'remove' and $hashPoint->{$target}{$position}{'Type'}=~/\/rs\d+/); # 排除SNP位点
                    my $pct = (exists($hashMethyData->{$target}{$position}{$sample}) and $hashMethyData->{$target}{$position}{$sample}{'Allcount'}>=10) ? $hashMethyData->{$target}{$position}{$sample}{'Methycount'}/$hashMethyData->{$target}{$position}{$sample}{'Allcount'}:"";
                    if($pct =~ /\d+/){$sum+=$pct;$num++;}   
                }           
            }
            my $mean = ($num > 0) ? $sum/$num : "";
            push @datas_gene, $mean;
        }
        my $datas_gene = join"\t", @datas_gene;
        print OUT "$datas_gene\n";
    }
    close OUT;
}

sub generate_efficiency_file{
    my $allsamples      = shift @_;
    my $hashEfficiency  = shift @_;
    my $efficiency_file = shift @_;
    my @title = ("Sample","Transferred(C->T)","All Base(C)","Efficiency");
    my $title = join"\t",@title;
    open OUT, ">$efficiency_file\n";
    print OUT "$title\n";
    foreach my $sample(@$allsamples){
        my @datas = ();
        push @datas, $sample;
        push @datas, $hashEfficiency->{$sample}{'C2T'};
        push @datas, $hashEfficiency->{$sample}{'AllPositionCount'};
        push @datas, $hashEfficiency->{$sample}{'PCT'};
        my $datas = join"\t", @datas;
        print OUT "$datas\n";
    }
    close OUT;
}

sub read_analysis_point{
    my $c_point_file = shift @_;
    my $hashPoint    = shift @_;
    my @head;
    my $count=0;
    open C_POINT,$c_point_file;
    while(<C_POINT>){
        $_=~s/[\r\n]//g;$count++;
        my @data = split /\t/, $_;
        if($count == 1){
           @head = @data;
           next;
        }
        my $target   = $data[0];
        my $position = $data[1];
        foreach my $col(0..$#head){
            my $value = exists( $data[$col] ) ? $data[$col] : "";           
            $hashPoint->{$target}{$position}{$head[$col]} = $value;
        }
    }
    close C_POINT;
}

sub read_haplotype{
    my $methyl_dir    = shift @_;
    my $samples       = shift @_;
    my $hashHaplotype = shift @_;
    print "Reading Sample.haplotype ...";
    foreach my $sample(@$samples){
        open HAP,"$methyl_dir/$sample.haplotype";
        while(<HAP>){
            $_=~s/[\r\n]//g;
            next if($_!~/\w/);
            my($target, $haplotype, $reads) = split /\t/, $_;
            $hashHaplotype->{'hap1'}{$target}{$haplotype}{$sample} = $reads;  # 保存每个片段，每个单倍型在每个样本上的reads数量
            $hashHaplotype->{'hap2'}{$target}{$sample}{$haplotype} = $reads;  # 键值顺序变化
            $hashHaplotype->{'sample_sum'}{$target}{$sample}       += $reads; # 每个片段上，样本的所有reads
            $hashHaplotype->{'hap_sum'}{$target}{$haplotype}       += $reads; # 每个片段上单倍型总和
        }
        close HAP;
    }
    print "OK\n";
}

sub read_methy_base{
    my $methyl_dir     = shift @_;
    my $allsamples     = shift @_;
    my $hashMethyData  = shift @_;
    my $hashEfficiency = shift @_;
    print "Reading Sample.base ...";
    # 读甲基化数据
    foreach my $sample (@$allsamples){
        open METH_BASE, "$methyl_dir/$sample.base";
        while(<METH_BASE>){
            $_ =~ s/[\r\n]//g;
            next if($_ !~ /\w/);
            if($_ =~ /^summary\(/){
               my @array = split /\t/,$_;
               shift @array;
               my %hashthis;
               map{my ($name,$value) = split /=/,$_; $hashthis{$name} = $value} @array;
               my $pct = $hashthis{'AllPositionCount'} > 0 ? sprintf("%.2f",($hashthis{'AllPositionCount'} - $hashthis{'MethyC'})*100/$hashthis{'AllPositionCount'}) : 0;
               $hashEfficiency->{$sample}{'C2T'}              = $hashthis{'AllPositionCount'} - $hashthis{'MethyC'};
               $hashEfficiency->{$sample}{'AllPositionCount'} = $hashthis{'AllPositionCount'};
               $hashEfficiency->{$sample}{'PCT'}              = $pct;       
               next;
            }
            my ($target, $position, $type, $methycount, $allcount) = split /\t/,$_;
            $hashMethyData->{$target}{$position}{$sample}{'Methycount'} = $methycount;
            $hashMethyData->{$target}{$position}{$sample}{'Allcount'}   = $allcount;
        }
        close METH_BASE;
    }
    print "OK\n";
}

sub get_sample_firstn{
    my $hashHaplotype  = shift @_;
    my $samples        = shift @_;
    my $hashSampleHapN = shift @_;
    my $n              = shift @_;  
    foreach my $target(keys %{$hashHaplotype->{'hap2'}}){
        my %hashChoose;# 最终筛选出的单倍型
        foreach my $sample(@$samples){
            next if(!exists($hashHaplotype->{'hap2'}{$target}{$sample}));
            my $count = 0;
            foreach my $haplotype(sort{ $hashHaplotype->{'hap2'}{$target}{$sample}{$b} <=> $hashHaplotype->{'hap2'}{$target}{$sample}{$a} } keys %{$hashHaplotype->{'hap2'}{$target}{$sample}}){
                $hashChoose{$haplotype}++;
                $count++;
                last if($count>=$n);# 列出前n种
            }
        }
        #按总reads数排序
        my @chooses;
        foreach my $haplotype(sort { $hashHaplotype->{'hap_sum'}{$target}{$b} <=> $hashHaplotype->{'hap_sum'}{$target}{$a} } keys %{$hashHaplotype->{'hap_sum'}{$target}}){
            push @chooses, $haplotype if(exists $hashChoose{$haplotype});
        }
        $hashSampleHapN->{$target} = join ",",@chooses;
    }
}

sub mark_cg{
    my $target_seq = shift @_;
    my $hashPoint  = shift @_;
    my $format     = shift @_;
    my @alleles    = split //,$target_seq;
    my $position=0;
    my @formats;
    my @tmps;
    my $type="";
    foreach my $allele(@alleles){
        $position++;
        my $nowtype = exists($hashPoint->{$position}) ? "redbold" : "normal";# 当前位点颜色
        if($nowtype eq $type){
            push @tmps, $allele;
            next;
        }
        if(@tmps==0){
            @tmps = ($allele);
            $type = $nowtype;
            next;
        }else{
            my $string = join "", @tmps;
            push @formats, $format->{$type};
            push @formats, $string;
            @tmps = ($allele);
            $type = $nowtype;
        } 
    }
    if(@tmps > 0){
        my $string = join "",@tmps;
        push @formats, $format->{$type};
        push @formats, $string;     
    }
    return @formats;
}

# 读取文件
sub read_file{
    my $file = shift @_;
    my %hashData;

    my $row = 0;
    open INPUT, $file;
    my $line0 = <INPUT>;
       $line0 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line0;

    $row++;
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        foreach my $col(0..$#heads)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashData{'Value'}{$row}{$heads[$col]} = $value;
        }
        $row++;
    }
    close INPUT;

    $hashData{'HEAD'} = join "\t", @heads;
    return %hashData;
}

sub check_sample_result{
    my $hashConfig = shift @_;
    my $allsamples = shift @_;
    print "Check Sample Set...";
    my $warnings = "";
    my $isOK_sample = 1;
    foreach my $sample(@$allsamples){
        my $methyl_dir = "$hashConfig->{'Report'}/methylation";
        $warnings.="$sample," if(package::utils::is_file_ok("$methyl_dir/$sample.base")==0 or package::utils::is_file_ok("$methyl_dir/$sample.haplotype")==0);
    }
    if($warnings eq ""){
        print "OK\n";
    }else{
        print "\n[ERR] Loss Methylation Result, Sample:$warnings\n";
        $isOK_sample=0;
    }
    return $isOK_sample;
}

# 输出readme
sub readme{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $readme_file = shift @_;
    # package路径
    my $package_dir = Cwd::abs_path(package::utils::get_dirname(File::Spec->rel2abs(__FILE__)));

    my $sheet = $workbook->add_worksheet('ReadMe');
       $sheet->set_row(0, 65);
       $sheet->set_column('A:A', 35);
       $sheet->set_column('B:B', 25);
       $sheet->set_column('C:C', 110);

    my $row = 0;
    open FILE,"$package_dir/$readme_file";
    while(<FILE>){
        $_=~ s/\^/\n/g;
        my @datas = split /\t/, $_;
        my $col = 0;
        foreach my $data(@datas)
        {
            my $text = Encode::decode("UTF-8",$data);
            $sheet->write($row, $col, $text, $format->{'readme1'})    if($row == 0 and $col == 0);
            $sheet->write($row, $col, $text, $format->{'readme2'})    if($row == 0 and $col == 1);
            $sheet->write($row, $col, $text, $format->{'readme2tmp'}) if($row == 0 and $col == 2);

            $sheet->write($row, $col, $text, $format->{'readme3'})    if($row >  0 and $col == 0);
            $sheet->write($row, $col, $text, $format->{'readme4'})    if($row >  0 and $col == 1);
            $sheet->write($row, $col, $text, $format->{'readme5'})    if($row >  0 and $col == 2);
            $col++;
        }
        $row++;
    }
    close FILE;   
}

1
