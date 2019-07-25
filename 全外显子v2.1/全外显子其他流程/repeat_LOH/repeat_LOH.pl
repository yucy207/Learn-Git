use package::main;
use package::check;
use package::format;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use package::annotation_db;
use Encode;
use strict;
use warnings;
use Statistics::R;
$|=1;

run();

sub run {
	die "Usage perl $0 Config\n" if(@ARGV < 1);
	my $config_file    = shift @ARGV;
	my %hashPara       = package::main::get_para();
	my %hashConfig     = package::main::read_config($config_file);
	
	package::check::repeat_LOH(\%hashPara, \%hashConfig);
	my @samples        = package::main::get_sample(\%hashConfig, 'case', 'control');
	my $species        = $hashConfig{'Species'};
	my $report_dir     = $hashConfig{'Report'};
    my $analysis_all   = "$report_dir/analysis_all.txt";
	my $all_library    = "$report_dir/library/library";
    my $database_snv   = "$report_dir/database.snv";
	my $repeat_loh_dir = "$report_dir/Repeat_LOH";
	my $document_dir   = "$report_dir/document";
	my $cnv_dir        = "$document_dir/3_CNV";
	my $out_dir        = "$cnv_dir/repeat_LOH_plot";
	package::main::make_dir($repeat_loh_dir);
	package::main::make_dir($document_dir);
	package::main::make_dir($cnv_dir);
	package::main::make_dir($out_dir);
	
	my $step           = 10000;
	my $window         = 50000;
	my %hash_analysis  = read_analysis_all($analysis_all);
    my %hash_database  = read_database_snv($database_snv);
	
	# # 重复区域判定
	print "Analyzing Repeat Region ...\n";
	my $greater_than    = 1.25;
	my $less_than       = 0.75;
	my $min_snp_num     = 10;
	my $right_snp_ratio = 0.8;
	my $mutation_type   = 'HET';
	my $mutation        = 'Repeat';	
	my @requirement     = ($greater_than, $less_than, $min_snp_num, $right_snp_ratio, $mutation_type, $mutation);
	run_analysis(\%hashConfig, \%hashPara, $repeat_loh_dir, $cnv_dir, $all_library, $species, $step, $window, \%hash_analysis, \%hash_database, \@samples, \@requirement);
	print "Repeat Region Analysis done.\n\n";
	
	# # LOH区域判定
	print "Analyzing LOH Region ...\n";
	$greater_than    = 10;
	$less_than       = 0.1;
	$min_snp_num     = 20;
	$right_snp_ratio = 0.9;
	$mutation_type   = 'ALL'; # # 是突变位点即可，包括HET和HOMA
	$mutation        = 'LOH';
    @requirement     = ($greater_than, $less_than, $min_snp_num, $right_snp_ratio, $mutation_type, $mutation);
	run_analysis(\%hashConfig, \%hashPara, $repeat_loh_dir, $cnv_dir, $all_library, $species, $step, $window, \%hash_analysis, \%hash_database, \@samples, \@requirement);	
	print "LOH Region Analysis done.\n\n";
	
	repeat_loh_overlap_plot_bed(\%hashPara, $repeat_loh_dir, \@samples);
    
    # # 绘图
    print "Repeat LOH Region Plot...\n";
    my $bed_for_plot_dir  = "$repeat_loh_dir/bed_for_plot";
    my $addLength         = 500000;# 区域绘图左右两边扩展长度
	
	if(@samples > 0) { # # 注释运行
        my $threshold = exists $hashPara{"Process"}{"repeat_LOH"} ? $hashPara{"Process"}{"repeat_LOH"} : 20;
		   $threshold = $hashConfig{"Process_repeat_LOH"} if(exists $hashConfig{"Process_repeat_LOH"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::main::process_bar_array($sample, \@samples)});# 进度条
        foreach my $i(0 .. @samples-1) {
            $pm->start($samples[$i]) and next;
			plot (\%hashPara, $species, $bed_for_plot_dir, $out_dir, $samples[$i], $addLength);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else {
        print "[Note] Process None, for reProcess Delete Result\n";
    }
    print "Plot done.\n\n";
}

sub plot{
    my $hashPara         = shift @_;
	my $species          = shift @_;
	my $bed_for_plot_dir = shift @_;
	my $out_dir          = shift @_;
    my $sample           = shift @_;
    my $addLength        = shift @_;
	my $out_sample_dir   = "$out_dir/$sample";
	package::main::make_dir($out_sample_dir);
    my $Rbin             = $hashPara -> {'Soft'}{'R'};
    my $RLib             = $hashPara -> {'Soft'}{'RLib'};
	my $dict_file        = $hashPara -> {$species}{'Dict'};
	my %chr_end          = read_dict($dict_file);
    my $R = Statistics::R->new(bin => $Rbin);
    $R -> startR;	
    $R->send(qq`.libPaths("$RLib") \n`);
    $R->send(qq` library(karyoploteR) \n`) ; 
    $R->send(qq` library(Gviz) \n`) ; 
    my $sample_plot_dir    = "$bed_for_plot_dir/$sample";
	my $original_data_dir  = "$sample_plot_dir/original_data";
    my $chrlist            = "$sample_plot_dir/chr.list";
    my %chrs               = readlist ($chrlist);   
    foreach my $chr (sort keys %chrs)
    {
        print "$sample:$chr plot.....";
        # 原始数据存取
        my $repeat_bed  = "$original_data_dir/$chr-Repeat";
		my $loh_bed     = "$original_data_dir/$chr-LOH";
        my $normal      = "$original_data_dir/$chr-normal";
		my $overlap_bed = "$original_data_dir/$chr-Overlap";
        my %allSNP      = ();
        readbedfile ($repeat_bed, \%allSNP);
        readbedfile ($loh_bed, \%allSNP);
        readbedfile ($overlap_bed, \%allSNP) if (-e $overlap_bed);
        readnormal  ($normal, \%allSNP);
        my %overlap     = readoverlop ($overlap_bed);
        # 绘图数生成
        my $plotDataDir = "$sample_plot_dir/plot_data";          
        mkdir $plotDataDir if (not -e $plotDataDir);
        my $regionRepeat  = "$plotDataDir/$chr-regionRepeat";
        my $regionLoh     = "$plotDataDir/$chr-regionLOH";
        my $regionOverlop = "$plotDataDir/$chr-regionOverlop";
        my $pointDatas    = "$plotDataDir/$chr-points";
        my %regionfiles   = generatePlotData(\%allSNP, \%overlap, \%chr_end, $addLength ,$regionRepeat, $regionLoh, $regionOverlop, $pointDatas, $plotDataDir, $chr);
        #
        my $sumCount    = keys %{$regionfiles{'c'}};
        my $pdf         = "$sample_plot_dir/$sample.$chr.pdf";          
        $R -> send(qq` pdf("$pdf",width=15) \n`);
        #总图 #red blue black
        $R -> send(qq` plot.params <- getDefaultPlotParams(plot.type=2) \n`);
        $R -> send(qq` plot.params\$ideogramheight <- 10 \n`);
        $R -> send(qq` plot.params\$data2height <- 600 \n`);
        $R -> send(qq` plot.params\$data1height <- 100 \n`);
        $R -> send(qq` plot.params\$bottommargin <- 0 \n`);
        $R -> send(qq` kp <- plotKaryotype(chromosomes = c("chr$chr"),plot.type = 2, plot.params=plot.params) \n`);
        $R -> send(qq` kpAddBaseNumbers(kp)  \n`);
        $R -> send(qq` repeatFile  = makeGRangesFromDataFrame(read.table("$regionRepeat", head = T)) \n`) if (counts($regionRepeat) >= 2);
        $R -> send(qq` LOHFile     = makeGRangesFromDataFrame(read.table("$regionLoh", head = T))  \n`) if (counts($regionLoh) >= 2);
        $R -> send(qq` OverlopFile = makeGRangesFromDataFrame(read.table("$regionOverlop", head = T))  \n`) if (counts($regionOverlop) >= 2);#判断行数是否大于2
        $R -> send(qq` point_loh = read.table("$pointDatas-loh", head = T)  \n`) if (counts("$pointDatas-loh") >= 2);
        $R -> send(qq` point_repeat = read.table("$pointDatas-repeat", head = T)  \n`) if (counts("$pointDatas-repeat") >= 2);
        $R -> send(qq` point_overlap = read.table("$pointDatas-overlap", head = T)  \n`) if (counts("$pointDatas-overlap") >= 2);
        $R -> send(qq` point_normal = read.table("$pointDatas-normal", head = T)  \n`);
        $R -> send(qq` kpPlotRegions(kp, repeatFile, col="red",data.panel=1, r1=0.5)  \n`) if (counts($regionRepeat) >= 2);
        $R -> send(qq` kpPlotRegions(kp, LOHFile, col="blue",data.panel=1, r1=0.5)  \n`) if (counts($regionLoh) >= 2);
        $R -> send(qq` kpPlotRegions(kp, OverlopFile, col="black",data.panel=1, r1=0.5)  \n`) if (counts($regionOverlop) >= 2);
        $R -> send(qq` kpAxis(kp, ymin=-10, ymax = 10, data.pane=2, r1=0.8)  \n`);
        $R -> send(qq` kpPoints(kp, chr = point_repeat\$chr,  x= point_repeat\$pos,  y= point_repeat\$value, ymin=-10, ymax=10, pch=16, cex=0.3, data.panel=2, col = "red")  \n`) if (counts("$pointDatas-repeat") >= 2);
        $R -> send(qq` kpPoints(kp, chr = point_loh\$chr,     x= point_loh\$pos,     y= point_loh\$value,    ymin=-10, ymax=10, pch=16, cex=0.3, data.panel=2, col = "blue")  \n`) if (counts("$pointDatas-loh") >= 2);
        $R -> send(qq` kpPoints(kp, chr = point_overlap\$chr, x= point_overlap\$pos, y= point_overlap\$value,ymin=-10, ymax=10, pch=16, cex=0.3, data.panel=2, col = "brown")  \n`) if (counts("$pointDatas-overlap") >= 2);
        $R -> send(qq` kpPoints(kp, chr = point_normal\$chr,  x= point_normal\$pos,  y= point_normal\$value, ymin=-10, ymax=10, pch=16, cex=0.2, data.panel=2, col = "grey")  \n`);
        foreach my $i (1..$sumCount)
        {
            my $filename = $regionfiles{'c'}{$i}{'file'};
            my ($start, $end) =split /\-/, $filename;
            my $starnew    = $start - $addLength; 
               $starnew    = 0 if ($starnew < 0);
            my $endnew     = $end + $addLength;
			    $endnew     = $chr_end{$chr} if(exists $chr_end{$chr} and $chr_end{$chr} < $endnew);
            my $type       = $regionfiles{'c'}{$i}{'type'};
            my $color      = "red";
               $color      = "blue" if ($type =~ /loh/ig);
            my $file       = "$sample_plot_dir/plot_data/$chr\_region/$filename\_$type.txt";
            my $filenormal = "$sample_plot_dir/plot_data/$chr\_region/$filename\_normal.txt";
            $R -> send(qq` ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr$chr") \n`);
            $R -> send(qq` axisTrack <- GenomeAxisTrack() \n`);
            $R -> send(qq` data1 = read.table("$file",header = T) \n`);  
            $R -> send(qq` data2 = read.table("$filenormal",header = T) \n`) if (counts($filenormal) >= 2); 
            $R -> send(qq` dtrack1 = DataTrack(data = data1\$log2r,start = data1\$pos,end = data1\$pos, chr = "chr$chr", genome = "hg19",col = "$color", name = "$filename\_$type") \n`);
            $R -> send(qq` dtrack2 = DataTrack(data = data2\$log2r,start = data2\$pos,end = data2\$pos, chr = "chr$chr", genome = "hg19",col = "gray") \n`) if (counts($filenormal) >= 2);
            $R -> send(qq` ot <- OverlayTrack(trackList = list(dtrack1)) \n`);
            $R -> send(qq` ot <- OverlayTrack(trackList = list(dtrack1,dtrack2)) \n`) if (counts($filenormal) >= 2);
            $R -> send(qq` ylims <- NULL \n`);
            $R -> send(qq` ylims <- extendrange(range(c(values(dtrack1), values(dtrack2)))) \n`) if (counts($filenormal) >= 2);
            $R -> send(qq` plotTracks(list(ideoTrack, axisTrack, ot),from = c($starnew) , to = c($endnew), ylim = ylims, showBandId = TRUE) \n`);
        }             
        $R -> send(qq` dev.off() \n`);
        system "cp $pdf $out_sample_dir";
		print "ok\n";
    }            
}

sub read_dict {
    my $dict_file = shift @_;
	my %chr_end;
	open DICT, $dict_file;
	while(<DICT>) {
	    next if($_ !~ /^\@SQ/); # # 跳过非染色体行
		$_ =~ s/[\r\n]//g;
	    $_ =~ /SN:(.*)\tLN:(\d+)\t/;
		my ($chr, $chr_length) = ($1, $2);
		$chr_end{$chr} = $chr_length;
	}
	close DICT;
	return %chr_end;
}

sub generatePlotData{
    my $allSNP       = shift @_;
    my $overlap      = shift @_;
	my $chr_end      = shift @_;
    my $addLength    = shift @_;
    my $regionRepeat = shift @_;
    my $regionLOH    = shift @_;
    my $regionOverlop= shift @_;
    my $pointDatas   = shift @_;
    my $plotDataDir  = shift @_;
    my $chr          = shift @_;   
    my %regionfiles  = ();
    #总图的散点图
    open POINT1, ">$pointDatas-overlap";
    open POINT2, ">$pointDatas-loh";
    open POINT3, ">$pointDatas-repeat";
    open POINT4, ">$pointDatas-normal";
    print POINT1 "chr\tpos\tvalue\n";
    print POINT2 "chr\tpos\tvalue\n";
    print POINT3 "chr\tpos\tvalue\n";
    print POINT4 "chr\tpos\tvalue\n";
    foreach my $snp (sort {$a<=>$b} keys %{$allSNP->{'all'}})
    {
        my $log2r    = "";
           $log2r    = $allSNP->{'all'}{$snp}{'Repeat'} if(exists $allSNP->{'all'}{$snp}{'Repeat'}) ; 
           $log2r    = $allSNP->{'all'}{$snp}{'LOH'} if(exists $allSNP->{'all'}{$snp}{'LOH'}) ; 
           $log2r    = $allSNP->{'all'}{$snp}{'Normal'} if(exists $allSNP->{'all'}{$snp}{'Normal'}) ; 
           $log2r    = 10 if ($log2r >10); 
           $log2r    = -10 if ($log2r < -10); 
        if (exists $overlap->{$snp})
        {        
            print POINT1 "chr$chr\t$snp\t$log2r\n";                    
        }
        elsif (exists $allSNP->{'all'}{$snp}{'LOH'})
        {
            print POINT2 "chr$chr\t$snp\t$log2r\n";
        }
        elsif ($allSNP->{'all'}{$snp}{'Repeat'})
        {  
            print POINT3 "chr$chr\t$snp\t$log2r\n";
        }
        else
        {
            print POINT4 "chr$chr\t$snp\t$log2r\n";
        }        
    } 
    close POINT1;
    close POINT2;
    close POINT3;
    close POINT4;
    
    #
    open REGIONR,   ">$regionRepeat";
    open REGIONLOH, ">$regionLOH";
    open REGIONO,   ">$regionOverlop";    
    print REGIONR "chr\tstart\tend\n";
    print REGIONLOH "chr\tstart\tend\n";
    print REGIONO "chr\tstart\tend\n";  
    
    my $region_dir = "$plotDataDir/$chr\_region";
    mkdir $region_dir if (not -e $region_dir);
    my $count = 0;  
    foreach my $start (sort {$a<=>$b} keys %{$allSNP->{'Region'}})
    {       
        my $minStart = $start -$addLength;
           $minStart = 0 if ($minStart < 0);
           my %tempStart = getRegionSNP($allSNP, $minStart, $start);
        foreach my $end (sort {$a<=>$b} keys %{$allSNP->{'Region'}{$start}})
        {
            my $maxEnd  = $end + $addLength;
			   $maxEnd  = $chr_end->{$chr} if(exists $chr_end->{$chr} and  $chr_end->{$chr} < $maxEnd);
            my %tempEnd  = getRegionSNP($allSNP, $end, $maxEnd);
            foreach my $name (sort keys %{$allSNP->{'Region'}{$start}{$end}})
            {
                 print REGIONR "chr$chr\t$start\t$end\n" if ($name =~ /Repeat/);
                 print REGIONLOH "chr$chr\t$start\t$end\n" if ($name =~ /LOH/);    
                 print REGIONO "chr$chr\t$start\t$end\n" if ($name =~ /Overlap/);    
                 next if ($name =~ /Overlap/);
                 $count ++;
                 $regionfiles{'c'}{$count}{'type'} = $name ;
                 $regionfiles{'c'}{$count}{'file'} = "$start-$end";
                 open OUT1, ">$region_dir/$start-$end\_$name.txt";
                 open OUT2, ">$region_dir/$start-$end\_normal.txt";
                 print OUT1 "pos\tlog2r\n";
                 print OUT2 "pos\tlog2r\n";
                 foreach my $pos (sort {$a <=> $b} keys %tempStart)
                {
                    my $log2r = $tempStart {$pos};
                    print OUT2 "$pos\t$log2r\n";
                }
                foreach my $snp (sort {$a<=>$b} keys %{$allSNP->{'Region'}{$start}{$end}{$name}})
                {
                    my $log2r = exists $allSNP->{'Region'}{$start}{$end}{'Repeat'}{$snp} ? $allSNP->{'Region'}{$start}{$end}{'Repeat'}{$snp} : $allSNP->{'Region'}{$start}{$end}{'LOH'}{$snp}; 
                    print OUT1 "$snp\t$log2r\n";                                  
                }
                foreach my $pos (sort {$a <=> $b} keys %tempEnd)
                {
                    my $log2r = $tempEnd {$pos};
                    print OUT2 "$pos\t$log2r\n";
                }                
                close OUT1;
                close OUT2;
           }
        }            
    }  
    close REGIONR;
    close REGIONLOH;
    close REGIONO;
    return %regionfiles;    
}
sub getRegionSNP{
    my $allSNP = shift @_;
    my $up     = shift @_;
    my $low    = shift @_;
    my %hash   = ();
    foreach my $snp (sort {$a<=>$b} keys %{$allSNP->{'all'}})
    {
        if (exists $allSNP->{'all'}{$snp}{'Normal'} and $snp >= $up and $snp <= $low)
        {
            my $log2r = $allSNP->{'all'}{$snp}{'Normal'};
            $hash{$snp} = $log2r;
        }  
    }
    return %hash;
}
sub readoverlop{
    my $file = shift @_;
    my %hash = ();
    return if (not -e $file);
    open IN, $file;
    while (<IN>)
    { 
           $_ =~ s/[\r\n]//g;
           my ($start, $end)   = (split /\t/, $_)[1,2,];           
           map{ $hash{$_} = 1 }($start..$end);
    }
    close IN;
    return %hash;
}
sub readnormal{
    my $normal = shift @_;
    my $allSNP = shift @_;
    open IN, $normal;
    my $count = 0;
    while (<IN>)
    {
          $count ++;
          next if ($count == 1);
          $_ =~ s/[\r\n]//g;
          my ($chr, $pos, $log2r) = split /\t/, $_;
          $allSNP -> {'all'}{$pos}{'Normal'} = $log2r;
    }
     close IN;
}
sub readbedfile{
    my $bed    = shift @_;
    my $allSNP = shift @_;
    return if (not -e $bed);
    open IN, $bed;
    while (<IN>)
    {
           $_ =~ s/[\r\n]//g;
           if ($bed =~ /Overlap/)
           {
              my ($chr, $start, $end, $temp) =split /\t/, $_, 4;
              $allSNP -> {'Region'}{$start}{$end}{'Overlap'} = 1;           
           }
           else
           {
               my ($chr, $start, $end, $snps, $depths, $log2r, $name) =split /\t/, $_;
               my @snp_all  = split /,/, $snps;
               my @r_all    = split /,/, $log2r;
               foreach my $i (0..$#snp_all)
               {
                   $allSNP -> {'all'}{$snp_all[$i]}{$name}       = $r_all[$i];                            
                   $allSNP -> {'Region'}{$start}{$end}{$name}{$snp_all[$i]}  =$r_all[$i];                            
               }                     
           }

    }
    close IN;
}


sub readlist{
    my $chrlist = shift @_;
    my %hash    = ();
    open IN,$chrlist;
    while (<IN>)
    {
          $_ =~ s/[\r\n]//;
          $hash{$_} = 1;          
    }close IN;
    return %hash;
}

sub repeat_loh_overlap_plot_bed {
    my $hashPara       = shift @_;
	my $repeat_loh_dir = shift @_;
	my $samples        = shift @_;
	my $bed_for_plot_dir   = "$repeat_loh_dir/bed_for_plot";
	foreach my $sample(@$samples) {
	    my $sample_plot_dir    = "$bed_for_plot_dir/$sample";
	    my $original_data_dir  = "$sample_plot_dir/original_data";
		my $bedtools           = $hashPara->{'Soft'}{'bedtools'};
        my $chrlist            = "$sample_plot_dir/chr.list";
        open OUT,">$chrlist";
		foreach my $chr(1..22,'X','Y') {
		    my $repeat_bed  = "$original_data_dir/$chr-Repeat";
			my $loh_bed     = "$original_data_dir/$chr-LOH";
            my $overlap_bed = "$original_data_dir/$chr-Overlap";
            print OUT "$chr\n" if(-e $repeat_bed or -e $loh_bed);
            if (-e $repeat_bed and -e $loh_bed) {
			    `$bedtools intersect -a $repeat_bed -b $loh_bed > $overlap_bed`;
			    `rm $overlap_bed` if(not -s $overlap_bed); # # 文件为空则删除
			}
		}
        close OUT;
	}
}

sub counts{
    my $file = shift @_;
    return 0 if (!(-e $file));
    my $info = `wc -l $file`;
    my $count = (split /\s+/, $info)[0];
    return $count ;
}

sub run_analysis {
    my $hashConfig     = shift @_;
	my $hashPara       = shift @_;
	my $repeat_loh_dir = shift @_;
	my $cnv_dir        = shift @_;
	my $all_library    = shift @_;
	my $species        = shift @_;
	my $step           = shift @_;
	my $window         = shift @_;
	my $hash_analysis  = shift @_;
	my $hash_database  = shift @_;
	my $samples        = shift @_;
	my $requirement    = shift @_;
	my $mutation       = $requirement->[-1];	
	my $mutation_dir   = "$repeat_loh_dir/$mutation";
	package::main::make_dir($mutation_dir);
	# # 区域筛选
	if(@$samples > 0) {
        my $threshold = exists $hashPara->{"Process"}{"repeat_LOH"} ? $hashPara->{"Process"}{"repeat_LOH"} : 20;
		   $threshold = $hashConfig->{"Process_repeat_LOH"} if(exists $hashConfig->{"Process_repeat_LOH"});
        my $pm = Parallel::ForkManager->new($threshold);
           $pm->run_on_start(sub{my ($pid, $sample) = @_; package::main::process_bar_array($sample, $samples)}); # # 进度条
        foreach my $i(0 .. @$samples - 1) {
            $pm->start($samples->[$i]) and next;
			sub_run($hashPara, $repeat_loh_dir, $mutation_dir, $all_library, $species, $step, $window, $requirement, $hash_analysis, $samples->[$i]);
            $pm->finish;    
        }
        $pm->wait_all_children;
    }
    else {
        print "[Note] Process None, check your input files\n";
    }
	# # 合并项目所有位点的注释文件
	combine_samples_annotation($samples, $hash_database, $mutation_dir, $repeat_loh_dir, $mutation);
	# # 输出
	my %hash_report  = read_samples_infos($samples, $mutation_dir, $mutation);
	my %hash_anno    = read_annotation($samples, $mutation_dir);
	my $workbook     = Excel::Writer::XLSX->new("$cnv_dir/$mutation.xlsx");
	my %format       = package::format::run($workbook);
	write_sample_sheet($workbook, \%format, \%hash_report, \%hash_anno); # # 样本表输出
	write_readme($workbook, \%format, "readme_$mutation.txt");
}

sub write_readme {
    my $workbook = shift @_;
	my $format   = shift @_;
	my $file   = shift @_;
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
    my $workbook           = shift @_;
	my $format             = shift @_;
	my $hash_report        = shift @_;
	my $hash_anno          = shift @_;
	my @titles    = ('Chr', 'Start', 'End', 'Length', 'SNV Count', 'Gene', 'Gene Regions', 'Type', 'Mean Ratio');
	foreach my $sample(keys %{$hash_report}) {
		my $row   = 0;
        my $sheet = $workbook->add_worksheet("$sample");
		$sheet->set_row(0, 60);
		foreach my $col(0 .. @titles - 1) { # # 输出标题
    	    $sheet->write($row, $col, $titles[$col], $format->{"title"});
    	}
		foreach my $chr(sort {$a cmp $b} keys %{$hash_report->{$sample}}) {
			foreach my $start(sort {$a <=> $b} keys %{$hash_report->{$sample}{$chr}}) {
			    foreach my $end(sort {$a <=> $b} keys %{$hash_report->{$sample}{$chr}{$start}}) {
				    $row ++;
					my $length        = $hash_report->{$sample}{$chr}{$start}{$end}{'Length'};
					my $type          = $hash_report->{$sample}{$chr}{$start}{$end}{'Type'};
				    my $snv_count     = $hash_report->{$sample}{$chr}{$start}{$end}{'SNV Count'};
					my $average_ratio = $hash_report->{$sample}{$chr}{$start}{$end}{'Mean Ratio'};
					my $gene          = exists $hash_anno->{$sample}{$chr}{$start}{'Gene'} ? $hash_anno->{$sample}{$chr}{$start}{'Gene'} : "";
					my $region        = exists $hash_anno->{$sample}{$chr}{$start}{'Gene Regions'} ? $hash_anno->{$sample}{$chr}{$start}{'Gene Regions'} : "";
					my @infos         = ($chr, $start, $end, $length, $snv_count, $gene, $region, $type, $average_ratio);
					foreach my $col(0 .. @infos - 1) {
    	                $sheet->write($row, $col, $infos[$col], $format->{'normal'});
    	            }
				}
			}
		}
	}
}

sub read_annotation {
    my $samples = shift @_;
	my $mutation_dir = shift @_;
	my %hash_anno;
	foreach my $sample(@$samples) {
	    my $anno_file = "$mutation_dir/$sample/library/library.variant_function";
		open ANNO, $anno_file;
	    while(<ANNO>) {
	        next if($_ !~ /\w/);
	        $_ =~ s/[\r\n]//g;
	    	my @infos = split /\t/, $_;
	    	my ($chr, $start)   = ($infos[2], $infos[3]);
	    	my ($region, $gene) = ($infos[0], $infos[1]);
	    	$hash_anno{$sample}{$chr}{$start}{'Gene Regions'} = $region;
	    	$hash_anno{$sample}{$chr}{$start}{'Gene'}         = $gene;
	    }
	    close ANNO;
	}
	return %hash_anno;
}

sub read_samples_infos {
    my $samples      = shift @_;
	my $mutation_dir = shift @_;
	my $mutation     = shift @_;
	my %hash_report;
	foreach my $sample(@$samples) {
	    my $result_file = "$mutation_dir/$sample/$mutation.txt";
		open RES, $result_file;
		while(<RES>) {
		    next if($_ !~ /\w/);
			$_ =~ s/[\r\n]//g;
			my ($chr, $start, $end, $length, $snv_count, $type, $average_ratio) = split /\t/, $_;
			$hash_report{$sample}{$chr}{$start}{$end}{'Length'}     = $length;
			$hash_report{$sample}{$chr}{$start}{$end}{'SNV Count'}  = $snv_count;
			$hash_report{$sample}{$chr}{$start}{$end}{'Type'}       = $type;
			$hash_report{$sample}{$chr}{$start}{$end}{'Mean Ratio'} = $average_ratio;
		}
		close RES;
	}
	return %hash_report;
}

sub combine_samples_annotation {
    my $samples         = shift @_;
	my $hash_database   = shift @_;
	my $mutation_dir    = shift @_;
	my $repeat_loh_dir  = shift @_;
	my $mutation        = shift @_;
	my %hash_anno       = read_loh_repeat_annotation($samples, $mutation_dir, $mutation);
	my $annotation_file = "$repeat_loh_dir/$mutation\_annotation.txt";
	write_annotation($hash_database, \%hash_anno, $annotation_file, $mutation);
}

sub write_annotation {
    my $hash_database   = shift @_;
	my $hash_anno       = shift @_;
	my $annotation_file = shift @_;
	my $mutation        = shift @_;
	open FILE, ">$annotation_file";
	foreach my $title(keys %{$hash_database}) {
	    my ($chr, $pos, $ref) = split /\|/, $title;
		my $start = (split /\-/, $pos)[0];
		foreach my $alt(keys %{$hash_database->{$title}}) {
		    if(exists $hash_anno->{$chr}{$start}{$ref}{$alt}) {
			    my $to_annotate = $hash_anno->{$chr}{$start}{$ref}{$alt};
				print FILE "$title\t$alt\t$mutation\t$to_annotate\n";
			}
		}
	}
	close FILE;
}

sub read_loh_repeat_annotation {
    my $samples      = shift @_;
	my $mutation_dir = shift @_;
	my $mutation     = shift @_;
	my %hash_anno;
	foreach my $sample(@$samples) {
	    my $anno_file = "$mutation_dir/snp_filter_DB/samples_annotation/library.hg19_$mutation\_$sample";
		open ANNO, $anno_file;
		while(<ANNO>) {
		   next if($_ !~ /\w/);
		   $_ =~ s/[\r\n]//g;
		   my ($mark, $infos, $chr, $start, $end, $ref, $alt, $tmp)= split /\t/, $_;
		   my $to_annotate = (split /=/, $infos)[1];
		   $hash_anno{$chr}{$start}{$ref}{$alt} .= "$to_annotate,";
		}
		close ANNO;
	}
	return %hash_anno;
}

sub sub_run {
    my $hashPara       = shift @_;
	my $repeat_loh_dir = shift @_;
	my $mutation_dir   = shift @_;
	my $all_library    = shift @_;
	my $species        = shift @_;
	my $step           = shift @_;
	my $window         = shift @_;
	my $requirement    = shift @_;
	my $hash_analysis  = shift @_;
	my $sample         = shift @_;
	my $mutation       = $requirement->[-1];
	my $sample_dir      = "$mutation_dir/$sample";
	package::main::make_dir($sample_dir);
	# # 获取窗口的突变信息
	my %regions  = filter_regions($hash_analysis, $step, $window, $requirement, $sample);
	# # 重复区域判定（根据 mutation > normal > unknown 优先级）
	my %regions_rm_overlap = rm_overlap_region($hashPara, $sample_dir, \%regions, $mutation, $sample);	
	# # 按样本进行结果记录、区域注释 及 绘图输入文件生成
	memo_n_annotation($hashPara, $species, $repeat_loh_dir, $sample_dir, $requirement, $hash_analysis, \%regions_rm_overlap, $sample);
	# # 项目snp区域注释 及 整理
	all_snv_annotation($hashPara, $species, $mutation_dir, $all_library, \%regions_rm_overlap, $mutation, $sample);
}

sub all_snv_annotation {
    my $hashPara           = shift @_;
	my $species            = shift @_;
	my $mutation_dir       = shift @_;
    my $all_library        = shift @_;
	my $regions_rm_overlap = shift @_;
	my $mutation           = shift @_;
	my $sample             = shift @_;
	my $AnnovarBuild       = $hashPara->{$species}{'AnnovarBuild'};
	my $snp_filter_DB      = "$mutation_dir/snp_filter_DB";
	my $DB_database        = "$snp_filter_DB/samples_database";
	my $DB_annotation      = "$snp_filter_DB/samples_annotation";
	package::main::make_dir($snp_filter_DB);
	package::main::make_dir($DB_database);
	package::main::make_dir($DB_annotation);
	my $db_file            = "$DB_database/$AnnovarBuild\_$mutation\_$sample.txt";
	open DB, ">$db_file"; # # 建库，用于项目的SNP位点筛选
	my $db_line_num = 1;
	foreach my $chr(sort {$a cmp $b} keys %{$regions_rm_overlap}) {
		foreach my $count(sort {$a <=> $b} keys %{$regions_rm_overlap->{$chr}}) {
			my $type   = $regions_rm_overlap->{$chr}{$count}{'type'};
			next if($type ne $mutation);
			my ($start, $end) = split /,/, $regions_rm_overlap->{$chr}{$count}{'pos'};
			my $to_annotate   = "$sample|$type:chr$chr-$start-$end";
			print DB "$db_line_num\t$chr\t$start\t$end\t$to_annotate\n";
			$db_line_num ++;
		}
	}
	close DB;	
	# # 注释
	my $all_library_new = "$DB_annotation/library";
	`ln -s $all_library $all_library_new` if(not -e $all_library_new);
	my $AnnovarDB    = $hashPara->{$species}{'AnnovarDB'};
	my $annovar_dir  = $hashPara->{'Soft'}{'AnnovarDIR'};
	`$annovar_dir/annotate_variation.pl -regionanno -dbtype $mutation\_$sample --buildver $AnnovarBuild $all_library_new $DB_database`;
}

sub memo_n_annotation {
    my $hashPara           = shift @_;
	my $species            = shift @_;
	my $repeat_loh_dir     = shift @_;
    my $sample_dir         = shift @_;
	my $requirement        = shift @_;
	my $hash_analysis      = shift @_;
	my $regions_rm_overlap = shift @_;
	my $sample             = shift @_;	
	my $mutation           = $requirement->[-1];
	my %mutation_pos; # # 记录被突变区覆盖的位点，用于后续未突变位点的记录
	my %infos_for_plot; # # 记录需要被绘图的位点及区域信息
	my $library_dir        = "$sample_dir/library";
	package::main::make_dir($library_dir);
	my $output_file        = "$sample_dir/$mutation.txt";
	my $library_file       = "$library_dir/library";
	open OUT, ">$output_file";  # # 记录单个样本结果，用于输出读取及合并
	open LIB, ">$library_file"; # # 生成library，用于建库
	foreach my $chr(sort {$a cmp $b} keys %{$regions_rm_overlap}) {
		foreach my $count(sort {$a <=> $b} keys %{$regions_rm_overlap->{$chr}}) {
			my $type   = $regions_rm_overlap->{$chr}{$count}{'type'};
			next if($type ne $mutation);
			my $length = $regions_rm_overlap->{$chr}{$count}{'length'};
			my ($start, $end) = split /,/, $regions_rm_overlap->{$chr}{$count}{'pos'};
			my %snv_infos     = get_covered_snv_infos($hash_analysis, \%mutation_pos, \%infos_for_plot, $sample, $chr, $start, $end, $requirement);
			my $snv_count     = $snv_infos{'mutation_snp_num'};
			my $average_ratio = $snv_infos{'average_ratio'};
			print OUT "$chr\t$start\t$end\t$length\t$snv_count\t$type\t$average_ratio\n";
			print LIB "$chr\t$start\t$end\t0\t-\t$sample\n";
		}
	}
    close OUT;
	close LIB;
	# # 注释
	my $AnnovarBuild = $hashPara->{$species}{'AnnovarBuild'};
	my $annovar_dir  = $hashPara->{'Soft'}{'AnnovarDIR'};
	
	my %hashAnnotationList = package::annotation_db::get_annotation(); # # 获取注释列表
	foreach my $anno_model(keys %hashAnnotationList) {
        next if($hashAnnotationList{$anno_model}{'From'} ne 'Annovar');
        my $AnnovarDB = $hashAnnotationList{$anno_model}{'DBDir'}{$species};
        `$annovar_dir/$hashAnnotationList{$anno_model}{'Para'} --buildver $AnnovarBuild $library_file $AnnovarDB`;
    }
	
	
	# # 绘图输入文件生成
	generate_normal_plot_file($hash_analysis, \%mutation_pos, \%infos_for_plot, $repeat_loh_dir, $sample); # # 获取所有区域及位点信息，返回染色体
	generate_plot_file(\%infos_for_plot, $repeat_loh_dir, $sample, $mutation);
}

sub generate_plot_file {
    my $infos_for_plot = shift @_;
	my $repeat_loh_dir = shift @_;
	my $sample         = shift @_;
	my $mutation       = shift @_;
	my $bed_for_plot_dir   = "$repeat_loh_dir/bed_for_plot";
	my $sample_plot_dir    = "$bed_for_plot_dir/$sample";
	my $original_data_dir  = "$sample_plot_dir/original_data";
	package::main::make_dir($bed_for_plot_dir);
	package::main::make_dir($sample_plot_dir);
	package::main::make_dir($original_data_dir);
	foreach my $chr(keys %{$infos_for_plot}) {
	    my $chr_file     = "$original_data_dir/$chr-$mutation";
		open CHR, ">$chr_file";
	    foreach my $start(sort {$a <=> $b} keys %{$infos_for_plot->{$chr}}) {
			foreach my $end(sort {$a <=> $b} keys %{$infos_for_plot->{$chr}{$start}}) {
			    my $snv_poses = $infos_for_plot->{$chr}{$start}{$end}{'mutation_snp_pos'};
				my $snv_depth = $infos_for_plot->{$chr}{$start}{$end}{'mutation_snp_depth'};
				my $snv_ratio = $infos_for_plot->{$chr}{$start}{$end}{'mutation_snp_log2'};
				my $type      = $infos_for_plot->{$chr}{$start}{$end}{'Type'};	
				print CHR "$chr\t$start\t$end\t$snv_poses\t$snv_depth\t$snv_ratio\t$type\n";
			}
		}
		close CHR;
	}
}

sub generate_normal_plot_file {
    my $hash_analysis  = shift @_;
	my $mutation_pos   = shift @_;
	my $infos_for_plot = shift @_;
	my $repeat_loh_dir = shift @_;
	my $sample         = shift @_;
	my $bed_for_plot_dir   = "$repeat_loh_dir/bed_for_plot";
	my $sample_plot_dir    = "$bed_for_plot_dir/$sample";
	my $original_data_dir  = "$sample_plot_dir/original_data";
	package::main::make_dir($bed_for_plot_dir);
	package::main::make_dir($sample_plot_dir);
	package::main::make_dir($original_data_dir);
	foreach my $chr(keys %{$hash_analysis->{$sample}}) {
	    my %poses; # # 用于多态位点的处理
	    next if(not exists $infos_for_plot->{$chr});              # # 突变区域已经在之前放入hash，不存在突变区域的染色体舍弃
		my $chr_file = "$original_data_dir/$chr-normal";
		open CHR, ">$chr_file";
		print CHR "chr\tpos\tlog2r\n";
	    foreach my $pos(keys %{$hash_analysis->{$sample}{$chr}}) {
		    next if(exists $mutation_pos->{$sample}{$chr}{$pos}); # # 排除突变区域的位点
		    foreach my $ref_alt(keys %{$hash_analysis->{$sample}{$chr}{$pos}}) {
			    next if(exists $poses{$pos});                      # # 多态位点的处理，只保留第一个
				$poses{$pos} ++;
				my $ratio = $hash_analysis->{$sample}{$chr}{$pos}{$ref_alt}{'ratio'};
				$ratio = 0.0009765625 if($ratio < 0.0009765625); # # log2上下限为±10
				my $log2  = log($ratio) / log(2);
			    print CHR "$chr\t$pos\t$log2\n";
			}
		}
		close CHR;
	}
}

sub get_covered_snv_infos {
    my $hash_analysis  = shift @_;
	my $mutation_pos   = shift @_;
	my $infos_for_plot = shift @_;
	my ($sample, $chr, $start, $end, $requirement) = @_;
	my ($greater_than, $less_than, $min_snp_num, $right_snp_ratio, $mutation_type, $mutation) = @$requirement;
	my %snv_infos;
	my %hash_ratio;
	$hash_ratio{'greater'}{'snp_num'} = 0;
	$hash_ratio{'less'}{'snp_num'}    = 0;
	foreach my $pos(sort {$a <=> $b} keys %{$hash_analysis->{$sample}{$chr}}) {
		next if($pos <= $start);
		last if($pos >= $end);
		#突变类型判断
		my ($if_pos_greater_ok, $if_pos_less_ok) = (0, 0); # # 用于多态位点的判断，只要有一个满足条件即可
		my $type  = "";
		my ($ratio, $all_ratio) = (0, 0);
		my ($depth, $all_depth) = (0, 0);
		foreach my $ref_alt(keys %{$hash_analysis->{$sample}{$chr}{$pos}}) {
			my $this_type  = $hash_analysis->{$sample}{$chr}{$pos}{$ref_alt}{'type'};
			my $this_ratio = $hash_analysis->{$sample}{$chr}{$pos}{$ref_alt}{'ratio'};
			my $this_depth = $hash_analysis->{$sample}{$chr}{$pos}{$ref_alt}{'depth'};
			$all_ratio     = $this_ratio;
			$all_depth     = $this_depth;
			if(($mutation_type eq 'ALL') or ($this_type eq $mutation_type)) {
			    $type      = $this_type;
				$ratio     = $this_ratio;
				$depth     = $this_depth;
				$if_pos_greater_ok = 1 if($this_ratio > $greater_than); # # 重复区域的深度比率判断条件
				$if_pos_less_ok    = 1 if($this_ratio < $less_than);
			}
		}
		if(($mutation_type eq 'ALL') or ($type eq $mutation_type)) {
		    $snv_infos{'mutation_snp_num'} ++;
		    if($if_pos_greater_ok == 1) {
			    $hash_ratio{'greater'}{'snp_num'} ++;
				$hash_ratio{'greater'}{'snp_ratio_sum'} += $ratio;				
			}
			if($if_pos_less_ok == 1) {
			    $hash_ratio{'less'}{'snp_num'} ++;
				$hash_ratio{'less'}{'snp_ratio_sum'} += $ratio;
			}
			$ratio = 1024 if($ratio > 1024);
			$ratio = 0.0009765625 if($ratio < 0.0009765625); # # log2上下限为±10
			my $log2 = log($ratio) / log(2);
			$mutation_pos->{$chr}{$pos} ++; # # 记录被突变区覆盖的位点，用于后续未突变位点的记录
			$infos_for_plot->{$chr}{$start}{$end}{'mutation_snp_pos'}   .= "$pos,";
			$infos_for_plot->{$chr}{$start}{$end}{'mutation_snp_log2'}  .= "$log2,";
			$infos_for_plot->{$chr}{$start}{$end}{'mutation_snp_depth'} .= "$depth,";
			$infos_for_plot->{$chr}{$start}{$end}{'Type'}                = $mutation;
		}
	}
	my $chosen_type  = $hash_ratio{'greater'}{'snp_num'} >= $hash_ratio{'less'}{'snp_num'} ? 'greater' : 'less';
	my $chosen_ratio = $hash_ratio{$chosen_type}{'snp_ratio_sum'} / $hash_ratio{$chosen_type}{'snp_num'};
	$snv_infos{'average_ratio'} = $chosen_ratio;
	return %snv_infos;
}

sub rm_overlap_region {
    my $hashPara   = shift @_;
    my $sample_dir = shift @_;
	my $regions    = shift @_;
	my $mutation   = shift @_;
	my $sample     = shift @_;
	my $bed_dir    = "$sample_dir/bed_for_rm_overlap";
	package::main::make_dir($bed_dir);
	my %regions_rm_overlap;
    # # 把三种类型的区域写为bed文件用于去重	   
	my $target_region_file  = "$bed_dir/$sample\_target_ori.bed";
	my $normal_region_file  = "$bed_dir/$sample\_normal_ori.bed";
	my $unknown_region_file = "$bed_dir/$sample\_unknown_ori.bed";
	open TAR, ">$target_region_file";
	open NOR, ">$normal_region_file";
	open UNK, ">$unknown_region_file";
	foreach my $chr(sort {$a cmp $b} keys %{$regions}) {
		foreach my $start(sort {$a <=> $b} keys %{$regions->{$chr}}) {
		    my $end  = $regions->{$chr}{$start}{'end'};
			my $type = $regions->{$chr}{$start}{'type'};
			print TAR "$chr\t$start\t$end\t$type\n" if($type eq $mutation);
			print NOR "$chr\t$start\t$end\t$type\n" if($type eq 'normal');
			print UNK "$chr\t$start\t$end\t$type\n" if($type eq 'unknown');
		}
	}
	close TAR;
	close NOR;
	close UNK;
	# # 合并各类区域
	my $bedtools           = $hashPara->{'Soft'}{'bedtools'};
	my $target_merge_file  = "$bed_dir/$sample\_target_merge.bed";
	my $normal_merge_file  = "$bed_dir/$sample\_normal_merge.bed";
	my $unknown_merge_file = "$bed_dir/$sample\_unknown_merge.bed";
	`$bedtools merge -c 4 -o distinct -i $target_region_file > $target_merge_file`;
	`$bedtools merge -c 4 -o distinct -i $normal_region_file > $normal_merge_file`;
	`$bedtools merge -c 4 -o distinct -i $unknown_region_file > $unknown_merge_file`;
	# # 优先级低的区域去重(重叠区域类型优先级 target > normal > unknown)
	my $normal_rm_target_file         = "$bed_dir/$sample\_normal_rm_target.bed";
	my $unknown_rm_target_file        = "$bed_dir/$sample\_unknown_rm_target.bed";
	my $unknown_rm_target_normal_file = "$bed_dir/$sample\_unknown_rm_target_normal.bed";
	`$bedtools subtract -a $normal_merge_file -b $target_merge_file > $normal_rm_target_file`;
	`$bedtools subtract -a $unknown_merge_file -b $target_merge_file > $unknown_rm_target_file`;
	`$bedtools subtract -a $unknown_rm_target_file -b $normal_merge_file > $unknown_rm_target_normal_file`;
	# # 获取无重叠的区域	
	read_bed(\%regions_rm_overlap, $target_merge_file);
	# # read_bed(\%regions_rm_overlap, $normal_rm_target_file);
	# # read_bed(\%regions_rm_overlap, $unknown_rm_target_normal_file);
	
	return %regions_rm_overlap;
}

sub read_bed {
    my ($regions_rm_overlap, $bed_file) = @_;
	open BED, $bed_file;
	my $last_chr = "";
	my $count    = 0; # # 用于记录顺序，可提取前后两个区域
	while(<BED>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($chr, $start, $end, $type) = split /\t/, $_;
		$count    = 0 if($chr ne $last_chr);
		$last_chr = $chr;
		$regions_rm_overlap->{$chr}{$count}{'type'}    = $type;
		$regions_rm_overlap->{$chr}{$count}{'pos'}     = "$start,$end";
		$regions_rm_overlap->{$chr}{$count}{'length'}  = $end - $start;
		$count ++;
	}
	close BED;
}

sub filter_regions {
    my $hash_analysis = shift @_;
	my $step          = shift @_;
	my $window        = shift @_;
	my $requirement   = shift @_;
	my $sample        = shift @_;
	my ($greater_than, $less_than, $min_snp_num, $right_snp_ratio, $mutation_type, $mutation) = @$requirement;
	my %regions;
	# # 位点分到窗口中
	my %pos_group_window = group_pos2window($hash_analysis, $step, $window, $sample);
	# # 条件判断
	foreach my $chr(keys %pos_group_window) {
	    ##### # my $last_window_type        = "start";
		##### # my ($last_start, $last_end) = (-1, -1);
	    foreach my $count(sort {$a <=> $b} keys %{$pos_group_window{$chr}}) {
		    my @poses = exists $pos_group_window{$chr}{$count}{'poses'} ? (split /,/, $pos_group_window{$chr}{$count}{'poses'}) : ();
			my ($mutation_snp_num, $snp_greater_num, $snp_less_num) = (0, 0, 0);
			foreach my $pos(@poses) {
				my ($if_pos_greater_ok, $if_pos_less_ok)       = (0, 0); # # 用于多态位点的判断，只要有一个满足条件即可
			    my $type = "";
				foreach my $ref_alt(keys %{$hash_analysis->{$sample}{$chr}{$pos}}) {
					my $this_type  = $hash_analysis->{$sample}{$chr}{$pos}{$ref_alt}{'type'};
					my $ratio      = $hash_analysis->{$sample}{$chr}{$pos}{$ref_alt}{'ratio'};
					if(($mutation_type eq 'ALL') or ($this_type eq $mutation_type)) {
					    $type      = $this_type;
						$if_pos_greater_ok = 1 if($ratio > $greater_than); # # 重复区域的深度比率判断条件
						$if_pos_less_ok    = 1 if($ratio < $less_than);
					}
				}
				if(($mutation_type eq 'ALL') or ($type eq $mutation_type)) {
	                $mutation_snp_num ++;
	                $snp_greater_num ++ if($if_pos_greater_ok == 1);
	                $snp_less_num ++ if($if_pos_less_ok == 1);
	            }
			}
			my ($window_start, $window_end) = split /,/, $pos_group_window{$chr}{$count}{'window'};
			# # 窗口区域所属类型判断
			my $window_type;
			if($mutation_snp_num == 0 or $mutation_snp_num < $min_snp_num) {
				$window_type = 'unknown';
			}
			else {
				my $greater_ratio = $snp_greater_num / $mutation_snp_num;
			    my $less_ratio    = $snp_less_num / $mutation_snp_num;
				if(($greater_ratio > $right_snp_ratio) or ($less_ratio > $right_snp_ratio)) {
				    $window_type = $mutation;
				}
				else {
				    $window_type = 'normal';
				}
			}
			$regions{$chr}{$window_start}{'end'}  = $window_end;
			$regions{$chr}{$window_start}{'type'} = $window_type;
			##### # 与上一窗口合并(不能只考虑与上一窗口合并，而要考虑与有重叠的同类型窗口合并，直接用bedtools merge实现)	
			####if($window_type eq $last_window_type) { # # 如果可以合并，则只会改变end
			####    $regions{$sample}{$chr}{$last_start}{'end'}  = $window_end;
			####	$last_end = $window_end;
			####}
			####else {			    
			####	$regions{$sample}{$chr}{$window_start}{'end'}  = $window_end;
			####	$regions{$sample}{$chr}{$window_start}{'type'} = $window_type;
			####	$last_window_type = $window_type;
			####	($last_start, $last_end) = ($window_start, $window_end);
			####}
		}
	}
	return %regions;
}

sub group_pos2window {
    my $hash_analysis = shift @_;
	my $step          = shift @_;
	my $window        = shift @_;
	my $sample        = shift @_;
	my %pos_group_window;
	foreach my $chr(keys %{$hash_analysis->{$sample}}) {
		my @all_poses = sort {$a <=> $b} keys %{$hash_analysis->{$sample}{$chr}};
	    my ($min_pos, $max_pos) = ($all_poses[0], $all_poses[-1]);
	    # # my $max_count = int((($max_pos - $min_pos) - ($window - 1)) / $step; #不确定计算是否正确
		my $start_pos_count = 0;
	    foreach my $count(0 .. $max_pos) {
	        my $window_start = $min_pos + $count * $step;
	    	my $window_end   = $window_start + $window;
	    	foreach my $i($start_pos_count .. @all_poses - 1) {
			    my $pos      = $all_poses[$i];
	    	    if($pos >= $window_start and $pos < $window_end) {
				    $pos_group_window{$chr}{$count}{'poses'} .= "$pos,";
				}
				$pos_group_window{$chr}{$count}{'window'} = "$window_start,$window_end";
				my $next_window_start = $window_start + $step;
				$start_pos_count = $i if($pos <= $next_window_start); # # 下次循环的起始
				last if($pos >= $window_end);
	    	}		
	    	last if($window_end >= $max_pos); # # 循环数按照最大设定（步长为1，起始为0时循环max_pos次），此处才是循环次数限制条件
	    }
	}
	return %pos_group_window;
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

sub read_analysis_all {
    my $analysis_all = shift @_;
    my %hash_analysis;
    open ANA, $analysis_all;
	while(<ANA>) {
	    next if($_ !~ /\w/);
		$_ =~ s/[\r\n]//g;
		my ($title, $alt, $score, $homr, $het, $homa) = split /\t/, $_;
		my ($chr, $pos, $ref) = split /\|/, $title;
		$pos = (split /\-/, $pos)[0];
		my @infos_needed = ((split /,/, $het),(split /,/, $homa));
		foreach my $infos(@infos_needed) {
		    my ($sample, $type, $geno, $ref_depth, $alt_depth) = split /\:/, $infos;
			# # 为了之后做除法分母不为0，由于ref_depth不为0时最小值为1，若alt_depth=0，设为0.000001，
			# # 则ratio=1000000，均满足LOH和重复区域的过滤条件ratio>10和ratio>1.25(HET应没有depth=0这种情况)
			my $depth  = $ref_depth + $alt_depth;
			$alt_depth = ($alt_depth == 0) ? 0.000001 : $alt_depth; 
			my $ratio  = $ref_depth / $alt_depth;
			$hash_analysis{$sample}{$chr}{$pos}{"$ref|$alt"}{'type'}  = $type;
			$hash_analysis{$sample}{$chr}{$pos}{"$ref|$alt"}{'ratio'} = $ratio;
			$hash_analysis{$sample}{$chr}{$pos}{"$ref|$alt"}{'depth'} = $depth;
		}
	}
	close ANA;
	return %hash_analysis;
}