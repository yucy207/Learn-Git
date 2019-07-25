package package::plot;
use Encode;
use strict;
use warnings;
use Statistics::R;
$|=1;
sub run {
    my $hashPara   = shift @_;
	my $hashConfig = shift @_;    
    print "########## Start plot ".package::main::get_time()." ##########\n";
    my $region_orifile = $hashConfig->{'RealignBed'};
    my @cases          = package::main::get_all_sample($hashConfig, 'case');
    my @controls       = package::main::get_all_sample($hashConfig, 'control');
    my @samples        = (@cases, @controls);
    my $region_file    = $hashConfig->{'Report'}."/excavator2/merge_".(split /\//, $region_orifile)[-1]; #合并后的bed文件
    my %mergeBed       = readBed($region_file); 
    my %coverages      = ();
    my $excavator_dir  = $hashConfig->{'Report'}."/excavator2";
	my $document_dir   = $hashConfig->{'Report'}."/document";
	my $cnv_dir        = "$document_dir/3_CNV";
	my $all_plot_dir   = "$cnv_dir/excavator2_plot";
	my $out_dir        = "$all_plot_dir/CNV_DepthRatio_Plot";
	package::main::make_dir($document_dir);
	package::main::make_dir($cnv_dir);
	package::main::make_dir($all_plot_dir);
	package::main::make_dir($out_dir);
    map { readtag("$excavator_dir/status/$_"."_final.tag", $_, \%coverages)}@samples;  
    
    #运行
    my $Threshold      = @cases;
    my $pm             = Parallel::ForkManager->new($Threshold);
    foreach my $sample (@cases)
    {
        $pm -> start and next;
        subrun($sample, $hashConfig, $hashPara, \%mergeBed, \%coverages, $out_dir);
        $pm -> finish;   
    }
    $pm->wait_all_children;   
}

sub subrun{
    my $sample        = shift @_;
    my $hashConfig    = shift @_;
    my $hashPara      = shift @_;
    my $mergeBed      = shift @_; 
    my $coverages     = shift @_;
	my $out_dir       = shift @_;
    my @chrs          = sort keys %{$mergeBed};   
    my $report_dir    = $hashConfig->{'Report'};
    my $excavator_dir = "$report_dir/excavator2";
    my $bed_dir       = "$excavator_dir/bed_for_plot";    
    my $Rbin          = $hashPara -> {'Soft'}{'R'};
    my $RLib          = $hashPara -> {'Soft'}{'RLib'};
    my $R = Statistics::R->new(bin => $Rbin);
    $R -> startR;
    $R -> send(qq`.libPaths("$RLib") \n`);
    $R -> send(qq` library(karyoploteR) \n`) ; 
    $R -> send(qq` library(Gviz) \n`) ; 	
    my %caseCOV          = COV($sample, $coverages);
    while (my $sampleDir = glob("$bed_dir/$sample*"))
    {  
          my $dirName          = (split /\//, $sampleDir)[-1];
          my $controlsName     = (split /\.VS\./, $dirName)[1];
          my %controlCOV       = ();
          if ($controlsName =~ /control\d+/) #如果对照是多个样本的时候测序深度取均值
          { 
              my $controlList  = "$report_dir/document/3_CNV/control_groups.txt";
              meanCOV($controlList, $controlsName, $coverages, \%controlCOV);                        
          }
          else
          {
              %controlCOV = COV($controlsName, $coverages);
          }
          my $controlsNameNEW     = "VS.$controlsName";
          my $resultsDir          = "$excavator_dir/result_final";
          my $cnvFile             = "$resultsDir/$sample.$controlsNameNEW/$sample.$controlsNameNEW\_result_filter.txt";
          my $original_dataDir    = "$sampleDir/original_data";
          my $plot_DataDir        = "$sampleDir/plot_data";
		  my $out_case_suffix_dir = "$out_dir/$sample.$controlsNameNEW";
          mkdir $plot_DataDir if(not -e $plot_DataDir);
          package::main::make_dir($out_case_suffix_dir);		  
          foreach my $chr (@chrs){ print "$sample: $chr\n";
              #生成绘图数据
              my $original   = "$original_dataDir/$chr";
              
              #总图
              my $gainfile   = "$plot_DataDir/$chr-regionGain";
              my $lossfile   = "$plot_DataDir/$chr-regionLoss";
              my $normalfile = "$plot_DataDir/$chr-regionNormal";
              generateData($gainfile, $lossfile, $normalfile, $chr, $original, $mergeBed, \%caseCOV, \%controlCOV); 
              
              # 分区域
              my $regionDir = "$plot_DataDir/$chr"."_region";
              mkdir $regionDir if (not -e $regionDir);
              
              my %region_all = getReagionData($cnvFile, $regionDir, $chr, $mergeBed, \%caseCOV, \%controlCOV); 
              my @regions = sort {$a<=>$b} keys %region_all;  

              # 绘图
              my $pdf = "$sampleDir/$dirName.$chr.pdf";                             
              $R -> send(qq` pdf("$pdf",width=15) \n`);
              
              #总图
              $R -> send(qq` plot.params <- getDefaultPlotParams(plot.type=2) \n`);
              $R -> send(qq` plot.params\$topmargin <- 5 \n`);
              $R -> send(qq` plot.params\$outmargin <- 5 \n`);
              $R -> send(qq` plot.params\$data1height <- 60 \n`);
              $R -> send(qq` plot.params\$ideogramheight <- 20 \n`);
              $R -> send(qq` plot.params\$data2height <- 500 \n`);
              $R -> send(qq` plot.params\$bottommargin <- 50 \n`);
              $R -> send(qq` kp <- plotKaryotype(chromosomes = c("chr$chr"),plot.type = 2, plot.params=plot.params) \n`);
              $R -> send(qq` kpAddBaseNumbers(kp)  \n`);   
 
              $R -> send(qq` gainup = makeGRangesFromDataFrame(read.table("$gainfile", head = T)) \n`) if (counts($gainfile) >= 2);
              $R -> send(qq` lossup = makeGRangesFromDataFrame(read.table("$lossfile", head = T)) \n`) if (counts($lossfile) >= 2);
              
              $R -> send(qq` normalF = read.table("$normalfile", head = T) \n`) if (counts($normalfile) >= 2);
              $R -> send(qq` gainF   = read.table("$gainfile", head = T) \n`) if (counts($gainfile) >= 2);
              $R -> send(qq` lossF   = read.table("$lossfile", head = T) \n`) if (counts($lossfile) >= 2);

              $R -> send(qq` kpPlotRegions(kp, gainup, col="red",data.panel=1, r1=0.5)  \n`) if (counts($gainfile) >= 2);
              $R -> send(qq` kpPlotRegions(kp, lossup, col="blue",data.panel=1, r1=0.5)  \n`) if (counts($lossfile) >= 2);              
              $R -> send(qq` kpAxis(kp, ymin=0, ymax = 5, data.pane=2)  \n`);
              $R -> send(qq` kpSegments(kp, chr = normalF\$chr,  x0 = normalF\$start,  x1= normalF\$end, y0 = normalF\$value, y1 = normalF\$value, ymin = 0, ymax = 5, data.panel=2, col = "gray", border = "gray")  \n`) if (counts($normalfile) >= 2);
              $R -> send(qq` kpSegments(kp, chr = gainF\$chr,    x0 = gainF\$start,    x1= gainF\$end,   y0 = gainF\$value,   y1 = gainF\$value,   ymin = 0, ymax = 5, data.panel=2, col = "red",  border  = "red")  \n`) if (counts($gainfile) >= 2);
              $R -> send(qq` kpSegments(kp, chr = lossF\$chr,    x0 = lossF\$start,    x1= lossF\$end,   y0 = lossF\$value,   y1 = lossF\$value,   ymin = 0, ymax = 5, data.panel=2, col = "blue", border = "blue")  \n`) if (counts($lossfile) >= 2);

              # 区域图             
              if (@regions > 0)
              {   
                  foreach my  $count (@regions)
                 {   
                     my $filename        = $region_all{$count}{'name'};
                     my $cnv             = $region_all{$count}{'cnv'};
                     my $cov             = $region_all{$count}{'CASE_COV'};
                     my $cnvfile         = "$regionDir/$filename\_$cnv.txt";
                     my $normal          = "$regionDir/$filename\_normal.txt";
                     my ($start, $end)   = split /\-/, $filename;
                     my $length          = $end -$start;
                     my $starnew         = $start - $length;
                        $starnew         = 0 if ($starnew < 0);
                     my $endnew          = $end + $length;
                     my $color           = "red";
                        $color           = "blue" if ($cnv =~ /loss/ig);
                     $R -> send(qq` ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr$chr") \n`);
                     $R -> send(qq` axisTrack <- GenomeAxisTrack() \n`);                    
                     $R -> send(qq` data1 = read.table("$cnvfile",header = T) \n`); 
                     $R -> send(qq` data2 = read.table("$normal",header = T) \n`); 
                     $R -> send(qq` dtrack1 = DataTrack(data = data1\$value,start = data1\$start,end = data1\$end, chr = "chr$chr", genome = "hg19",type = "histogram",fill.histogram = "$color", col.histogram = "$color", name = "$filename\_$cnv ($cov)",  col.title = "black", cex.title = 1 ) \n`);#col.title
                     $R -> send(qq` dtrack2 = DataTrack(data = data2\$value,start = data2\$start,end = data2\$end, chr = "chr$chr", genome = "hg19",type = "histogram",fill.histogram = "gray",col.histogram = "gray", background.title = "darkblue", cex.title = 1) \n`) if (counts($normal) >= 2);
                     $R -> send(qq` ot <- OverlayTrack(trackList = list(dtrack1)) \n`);
                     $R -> send(qq` ot <- OverlayTrack(trackList = list(dtrack1,dtrack2)) \n`) if (counts($normal) >= 2);
                     $R -> send(qq` ylims <- NULL \n`);
                     $R -> send(qq` ylims <- extendrange(range(c(values(dtrack1), values(dtrack2)))) \n`) if (counts($normal) >= 2);                
                     $R -> send(qq` plotTracks(list(ideoTrack, axisTrack, ot),from = c($starnew) , to = c($endnew), ylim = ylims, showBandId = TRUE) \n`);               
                 }                           
             }
             $R->send(qq` dev.off() \n`) ;
			 system "cp $pdf $out_case_suffix_dir";
             print "ok\n";
          }
    }
    $R->stopR;	    
}

sub counts{
    my $file = shift @_;
    return 0 if (!(-e $file));
    my $info = `wc -l $file`;
    my $count = (split /\s+/, $info)[0];
    return $count ;
}

sub generateData{
    my $gainfile   = shift @_;
    my $lossfile   = shift @_;
    my $normalfile = shift @_;
    my $chr        = shift @_;
    my $original   = shift @_;
    my $mergeBed   = shift @_;
    my $caseCOV    = shift @_;
    my $controlCOV = shift @_;
    my %cnv        = readcnv ($original);   
    open GAIN,">$gainfile";
    open LOSS,">$lossfile";
    open NORMAL,">$normalfile";
    print GAIN "chr\tstart\tend\tvalue\n";
    print LOSS "chr\tstart\tend\tvalue\n";
    print NORMAL "chr\tstart\tend\tvalue\n";
 
    foreach my $start (sort {$a<=>$b} keys %{$mergeBed->{$chr}})
    {
        foreach my $end (sort {$a<=>$b} keys %{$mergeBed->{$chr}{$start}})
        {
            my $cnvinfo  = exists $cnv{"$chr|$start|$end"}{'INFO'} ? $cnv{"$chr|$start|$end"}{'INFO'} : "normal";  
            my $copynum  = exists $cnv{"$chr|$start|$end"}{'CPN'} ? $cnv{"$chr|$start|$end"}{'CPN'} : "NA";
            my $case_coverage    = $caseCOV    -> {"$chr|$start|$end"};             
            my $control_coverage = $controlCOV -> {"$chr|$start|$end"};          
            my $ratio            = $case_coverage / $control_coverage if ($control_coverage != 0);
               $ratio            = 0  if ($control_coverage == 0);
               $ratio            = 5 if ($ratio > 5);                    
            print NORMAL "chr$chr\t$start\t$end\t$ratio\n" if ($cnvinfo =~ /normal/ig);
            print GAIN "chr$chr\t$start\t$end\t$ratio\n" if ($cnvinfo =~ /gain/ig);
            print LOSS "chr$chr\t$start\t$end\t$ratio\n" if ($cnvinfo =~ /loss/ig);
         }           
     }
     close NORMAL;
     close GAIN;
     close LOSS;
}

sub getReagionData{
    my $cnvFile    = shift @_;  
    my $regionDir  = shift @_;
    my $chrmark    = shift @_;   
    my $mergeBed   = shift @_;    
    my $caseCOV    = shift @_;
    my $controlCOV = shift @_;
    my %hash = ();
    return %hash if(not package::main::is_file_ok($cnvFile));
	open IN,$cnvFile;
    my $count = 0;
    while (<IN>)
    {
           $_ =~ s/[\r\n]//g;
           my ($chr, $start_cnv, $end_cnv, $cpn, $cnvinfo) = (split /\t/,$_)[0,1,2,3,4,];
		   next if($cpn eq '2');
           next if ($chr ne $chrmark);
           $count ++ ;
           my $length = $end_cnv - $start_cnv;
           my $min = $start_cnv - $length;
              $min = 0 if ($min < 0);
           my $max = $end_cnv +$length;
           $hash{$count}{'name'} = "$start_cnv-$end_cnv";
           $hash{$count}{'cnv'}  = "$cnvinfo";
           my %cov = ();
           my $cnvfile    = "$regionDir/$start_cnv\-$end_cnv\_$cnvinfo.txt";
           my $normalfile = "$regionDir/$start_cnv\-$end_cnv\_normal.txt";
           
           open CNV,">$cnvfile";
           open NORMAL,">$normalfile";
           print CNV "start\tend\tvalue\n";
           print NORMAL "start\tend\tvalue\n";
           my $caseSUM   = 0;
           my $caseCount = 0;
           foreach my $start (sort {$a<=>$b} keys %{$mergeBed->{$chrmark}})
           {
               foreach my $end (sort {$a<=>$b} keys %{$mergeBed->{$chrmark}{$start}})
               {
                   if ($start >=$min and $end <= $max)
                   {
                       my $case_coverage    = $caseCOV    -> {"$chr|$start|$end"};             
                       my $control_coverage = $controlCOV -> {"$chr|$start|$end"};          
                       my $ratio            = $case_coverage / $control_coverage if ($control_coverage != 0);
                       $ratio               = 0  if ($control_coverage == 0);
                       $ratio               = 5 if ($ratio > 5);
                       $cov{$ratio}         = 1; 
                       if ($start >= $start_cnv and $end <= $end_cnv)
                       {
                           $caseSUM += $case_coverage; 
                           $caseCount ++;                            
                           print CNV "$start\t$end\t$ratio\n";                          
                       }
                       else
                       {
                          print NORMAL "$start\t$end\t$ratio\n";                                                                                                            
                       }                                      
                   }                            
               }
            }
            close CNV;
            close NORMAL;
            $hash{$count}{'CASE_COV'} = sprintf("%.2f", $caseSUM / $caseCount);
    }
    close IN;
    return %hash;
}

sub readcnv{
    my $file = shift @_;
    my %hash = ();
    if (-e $file)
    {
       open IN,$file;
       while (<IN>)
       {
              $_ =~ s/[\r\n]//g;
              my ($chr, $start, $end, $cpn, $info) = split /\t/, $_;
              $hash{"$chr|$start|$end"}{'CPN'}  = $cpn;
              $hash{"$chr|$start|$end"}{'INFO'} = $info;
       }close IN;
    }
    return %hash;
}

sub COV{
    my $sample    = shift @_;
    my $coverages = shift @_;
    my %hash      = ();
    map { $hash{$_}= $coverages ->{$_}{$sample} } keys %{$coverages};
    return %hash;
}
sub meanCOV{
    my $controlList  = shift @_;
    my $controlsName = shift @_;
    my $coverages    = shift @_;
    my $controlCOV   = shift @_;   
    my @samples      = ();
    open IN,$controlList;
    while (<IN>)
    {
          $_ =~ s/[\r\n\s+]//g;
          my ($sample, $controls) =split /=/, $_;
          @samples = (split /,/,$controls) if ($sample eq $controlsName);
    }
    close IN;
    foreach my $info (keys %{$coverages})
    {
        my $sum =0;
        map { $sum += $coverages ->{$info}{$_} } @samples;
        $controlCOV -> {$info} = $sum / @samples ;          
    }   
}

sub readtag{
    my $file      = shift @_;
    my $sample    = shift @_;
    my $coverages = shift @_;
    open IN,$file;
    my $count = 0;
    while (<IN>)
    {
          $count ++;
          next if ($count==1);          
          $_ =~ s/[\r\n]//g;
          my ($chr ,$start, $end, $normalized_cov) = (split /\t/,$_)[0,1,2,7,];
          $coverages->{"$chr|$start|$end"}{$sample}= $normalized_cov;
    }    
    close IN;
}

sub readBed{
    my $file     = shift @_;
    my %mergeBed = ();
    open IN,$file;
    while (<IN>)
    {
           $_ =~ s/[\r\n]//g;   
           my ($chr, $start, $end) = split /\t/, $_;
           $mergeBed{$chr}{$start}{$end}=1;
    }
    close IN;
    return %mergeBed;
}

1