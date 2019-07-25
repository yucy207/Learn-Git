package package::circos;
use strict;
use warnings;

sub run{
	my $hashConfig = shift @_;
	my $hashPara   = shift @_;
	my $plotDatas  = shift @_;
	my $report     = $hashConfig ->{'Report'};
	my $dataDir    = "$report/sv/circosData"; package::utils::make_dir($dataDir);
	my @samples    = package::utils::get_sample($hashConfig, "case", "control");
    my @case       = package::utils::get_sample($hashConfig, "case");
    my @control    = package::utils::get_sample($hashConfig, "control");
	generateData ($plotDatas, $dataDir, \@case, \@control);
	CIRCOSplot ($hashConfig, $hashPara, $dataDir, \@samples);
}

sub CIRCOSplot{
	my $hashConfig = shift @_;
	my $hashPara   = shift @_;
	my $dataDir    = shift @_;
	my $samples    = shift @_;
    my @all        = (@$samples, "case", "control");
	my $reportDir  = $hashConfig -> {'Report'};
	my $sv_dir     = "$reportDir/document/4_SV"; package::utils::make_dir($sv_dir);
	my $svImageDir = "$sv_dir/SVImage"; package::utils::make_dir($svImageDir);
	my $Species    = $hashConfig -> {'Species'};
	my $Circos     = $hashPara   -> {'Soft'}{'Circos'};  
	my $SVCircosKT = (exists($hashPara -> {$Species}{'SVCircosKT'})) ? $hashPara -> {$Species}{'SVCircosKT'} : "";
	if($SVCircosKT eq ''){
		print "No Circos Plot! Lost karyotype File\n";
		return;
	}
  # 准备绘图
	foreach my $sample (@all){
		my $CNVGain = "$dataDir/$sample.gain";
		my $CNVLoss = "$dataDir/$sample.loss";
		my $SV      = "$dataDir/$sample.sv";
		my $config  = "$dataDir/$sample.ini";
        next if (package::utils::is_file_ok($CNVGain)==0 and package::utils::is_file_ok($CNVLoss)==0 and package::utils::is_file_ok($SV)==0);
		open CONFIG,">$config";
		print CONFIG "karyotype = $SVCircosKT\n";# 基因组信息
		# 基因组显示配置
		print CONFIG "<ideogram> \n";
		print CONFIG "\t<spacing> \n";
		print CONFIG "\t\tdefault = 0.005r \n";# 设置圈图中染色体之间的空隙大小，以下设置为每个空隙大小为周长的 0.5% 
		print CONFIG "\t</spacing> \n";
		print CONFIG "\tradius            = 0.85r \n";# 设定 ideograms 的位置，以下设定 ideograms 在图离圆心的 90% 处 
		print CONFIG "\tthickness         = 60p   \n";# 设定 ideograms 的厚度，可以使用 r（比例关系） 或 p（像素）作为单位 
		print CONFIG "\tstroke_color      = black \n";# 设定 ideograms 轮廓的颜色及其厚度。如果没有该参数或设定其厚度为0，则表示没有轮廓
		print CONFIG "\tstroke_thickness  = 2p    \n";# 设定 ideograms 轮廓的颜色及其厚度。如果没有该参数或设定其厚度为0，则表示没有轮廓
		print CONFIG "\tshow_label        = yes   \n";# 设定是否显示 label 。 label 对应着 karyotype 文件的第 4 列
		print CONFIG "\tlabel_font        = bold  \n";# label字体
		print CONFIG "\tlabel_radius      = 1r + 130p \n";# 设定 label 的位置 
		print CONFIG "\tlabel_size        = 60    \n";# 设定 label 的大小 
		print CONFIG "\tlabel_parallel    = yes   \n";# 设定 label 的字体方向，yes 是易于浏览的方向。
		print CONFIG "\tshow_bands        = yes   \n"; 
		print CONFIG "\tfill_bands        = yes   \n"; 
		print CONFIG "\tband_transparency = 0     \n"; 
		print CONFIG "</ideogram>   \n"; 
		# 显示基因组刻度尺
		print CONFIG "chromosomes_units = 1000000   \n"; # 1u代表的长度
		print CONFIG "show_ticks        = yes       \n"; # 是否显示ticks
		print CONFIG "show_tick_labels  = yes       \n"; # 是否显示ticks标签
		print CONFIG "<ticks>   \n";  
		print CONFIG "\tradius           = 1r    \n";  # 设定 ticks 的位置 
		print CONFIG "\tcolor            = black \n";  # 设定 ticks 的颜色 
		print CONFIG "\tthickness        = 2p    \n";  # 设定 ticks 的厚度
		print CONFIG "\tmultiplier       = 1e-6  \n";  # 设定 ticks' label 的值的计算。将该刻度对应位置的值 * multiplier 得到能展示到圈图上的 label 值。
		print CONFIG "\tformat           = %d    \n";  # label 值的格式化方法; %d 表示结果为整数
		print CONFIG "\t<tick> \n";   # 添加一个刻度
		print CONFIG "\t\tspacing        = 30u    \n"; # 没30u（chromosomes_units）添加一个刻度   
		print CONFIG "\t\tsize           = 15p    \n"; #    
		print CONFIG "\t\tshow_label     = yes    \n"; #    
		print CONFIG "\t\tlabel_font     = bold    \n"; #   
		print CONFIG "\t\tlabel_size     = 30p    \n"; #    
		print CONFIG "\t\tlabel_offset   = 10p    \n"; #    
		print CONFIG "\t\tformat         = %d    \n"; #    
		print CONFIG "\t</tick>    \n";   
		print CONFIG "</ticks>  \n";
		#######
		# 绘制CNV结果，直方图
		#######
		# CNV Gain 红色
		print CONFIG "<plots>  \n";
		print CONFIG "\t<plot>   \n";
		print CONFIG "\t\ttype        = histogram    \n";
		print CONFIG "\t\tfile        = $CNVGain   \n";
		print CONFIG "\t\torientation = out   \n";# 方向向外
		print CONFIG "\t\textend_bin  = no   \n"; # 区域不自动连接
		print CONFIG "\t\tthickness   = 2p   \n"; # 边框宽度
		print CONFIG "\t\tr0          = 0.95r   \n"; # 直方图绘制范围
		print CONFIG "\t\tr1          = 0.99r   \n";
		print CONFIG "\t\tfill_color  = vdred   \n";# 填充色
		print CONFIG "\t\tcolor       = vdred   \n";# 边框色
		print CONFIG "\t\tmin         = 0   \n";# 坐标最小值
		print CONFIG "\t\tmax         = 1   \n";# 最表最大值
		print CONFIG "\t</plot>    \n";# 最表最大值
		# CNV Loss 蓝色
		print CONFIG "\t<plot>   \n";
		print CONFIG "\t\ttype        = histogram    \n";
		print CONFIG "\t\tfile        = $CNVLoss   \n";
		print CONFIG "\t\torientation = out   \n";# 方向向外
		print CONFIG "\t\textend_bin  = no   \n"; # 区域不自动连接
		print CONFIG "\t\tthickness   = 2p   \n"; # 边框宽度
		print CONFIG "\t\tr0          = 0.95r   \n"; # 直方图绘制范围
		print CONFIG "\t\tr1          = 0.99r   \n";
		print CONFIG "\t\tfill_color  = vdblue   \n";# 填充色
		print CONFIG "\t\tcolor       = vdblue   \n";# 边框色
		print CONFIG "\t\tmin         = 0   \n";# 坐标最小值
		print CONFIG "\t\tmax         = 1   \n";# 最表最大值
		print CONFIG "\t</plot>    \n";# 最表最大值   	    
		print CONFIG "</plots>  \n"; 
		#######
		# 绘制基因融合
		####### 
		print CONFIG "<links>	  \n"; 
		print CONFIG "\t<link>	  \n"; 
		print CONFIG "\tfile          = $SV 	 \n"; 
		print CONFIG "\tradius        = 0.94r 	 \n"; # 设置 link 曲线的半径
		print CONFIG "\tbezier_radius = 0r 	     \n"; # 设置贝塞尔曲线半径，该值设大后曲线扁平，使图像不太好看。
		print CONFIG "\tcolor         = black_a4 \n"; # 设置 link 曲线的颜色 
		print CONFIG "\tthickness     = 4 	     \n"; # 设置 link 曲线的厚度 
		print CONFIG "\t</link>	  \n"; 
		print CONFIG "</links>	  \n"; 


		print CONFIG "<image> \n"; 
		print CONFIG "\t<<include etc/image.conf>>  \n"; 
		print CONFIG "</image>  \n"; 
		print CONFIG "<<include etc/colors_fonts_patterns.conf>>  \n"; 
		print CONFIG "<<include etc/housekeeping.conf>>  \n"; 
		system("$Circos -conf $config -outputdir $svImageDir -outputfile $sample")	;			
	}
}    


sub generateData{
	my $plotDatas = shift @_;
	my $dataDir   = shift @_;
	my $case      = shift @_;
    my $control   = shift @_;
    my %allRegion = ();
	foreach my $sample(@$case, @$control){
		my $CNVGain = "$dataDir/$sample.gain";
		my $CNVLoss = "$dataDir/$sample.loss";
		my $SV      = "$dataDir/$sample.sv";
		open GAIN, ">$CNVGain";
		open LOSS, ">$CNVLoss";
		open SV,   ">$SV";
		foreach my $svtype (sort keys %{$plotDatas -> {$sample}})
		{
			foreach my $info (sort keys %{$plotDatas -> {$sample}{$svtype}})
			{
                my ($chr1, $pos1, $pos2, $chr2, $pos3, $pos4) = split /\t/, $info;
                $allRegion{'case'}{"gain"}{"hs$chr1\t$pos1\t$pos2\t1"}=1     if($svtype eq "DUP" and $sample ~~ @$case);
                $allRegion{'control'}{"gain"}{"hs$chr1\t$pos1\t$pos2\t1"}=1  if($svtype eq "DUP" and $sample ~~ @$control);
                $allRegion{'case'}{"loss"}{"hs$chr1\t$pos1\t$pos2\t1"}=1     if($svtype eq "DEL" and $sample ~~ @$case);
                $allRegion{'control'}{"loss"}{"hs$chr1\t$pos1\t$pos2\t1"}=1  if($svtype eq "DEL" and $sample ~~ @$control);                
                $allRegion{'case'}{'sv'}{"hs$chr1\t$pos1\t$pos2\ths$chr2\t$pos3\t$pos4"}=1    if (($svtype eq "BND" or $svtype eq "INV") and $sample ~~ @$case);
                $allRegion{'control'}{'sv'}{"hs$chr1\t$pos1\t$pos2\ths$chr2\t$pos3\t$pos4"}=1 if (($svtype eq "BND" or $svtype eq "INV") and $sample ~~ @$control);
				print GAIN "hs$chr1\t$pos1\t$pos2\t1\n" if ($svtype eq 'DUP');
				print LOSS "hs$chr1\t$pos1\t$pos2\t1\n" if ($svtype eq 'DEL');
				print SV   "hs$chr1\t$pos1\t$pos2\ths$chr2\t$pos3\t$pos4\n" if ($svtype eq "BND" or $svtype eq "INV");
			}
		}
		close GAIN;
		close LOSS;
		close SV;
	}
    foreach my $sampleinfo ("case", "control")
    {
        foreach my $svtype ("gain", "loss", "sv")
        {
            my $file = "$dataDir/$sampleinfo.$svtype";
            open OUT,">$file";            
            map {print OUT "$_\n"} sort keys %{$allRegion{$sampleinfo}{$svtype}} if (exists $allRegion{$sampleinfo} and  exists $allRegion{$sampleinfo}{$svtype} and $allRegion{$sampleinfo}{$svtype} =~ /\w/);
            close OUT;            
        }       
    }
}    
		
1