package package::filter;
use strict;
use warnings;
use Encode;
use Excel::Writer::XLSX;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "##########Start output Excel and Ciocos ".package::utils::get_time()." ##########\n";
    my $report_dir        = $hashConfig->{'Report'};
    my $SVdbtype          = exists $hashConfig->{'GenomeBed'} ? "SVdbWGS" : "SVdbWES"; #判断是WES 还是WGS
    my $document          = "$report_dir/document";
    my $sv_dir            = "$document/4_SV";
    my $excel             = "$sv_dir/SV.xlsx";
    # if(package::utils::is_file_ok($excel))
    # {
    #     print "[Note] Process None, for reProcess Delete Result\n";
    #     return;
    # }  
    my $SVreportDir       = "$report_dir/sv";
    my $svFile            = "$SVreportDir/sample.all.final.sv.vcf";
    my %hashSV            = package::utils::readSV($svFile, $hashConfig, $SVdbtype);  
    my %hashAnnotation    = annotation(\%hashSV, $hashPara, $hashConfig); #注释信息提取 
    datafilter (\%hashSV, \%hashAnnotation, $hashPara, $hashConfig, $SVdbtype);#数据过滤，给 hashSV 增加过滤项
    writeSV (\%hashSV, \%hashAnnotation, $hashConfig, $hashPara, $SVdbtype); # EXCEL 结果输出
}

sub writeSV {
    my $hashSV          = shift @_;
    my $hashAnnotation  = shift @_;
    my $hashConfig      = shift @_;
    my $hashPara        = shift @_;
    my $SVdbtype        = shift @_;
    my $SUMin           = (exists $hashConfig->{'MinReads'}) ? $hashConfig->{'MinReads'} : 5;# 对变异结果的证据Reads数不得低于5重(>=5)
    my $is_output_raw   = (exists $hashConfig->{'OUTPUT_RAW'}) ? $hashConfig->{'OUTPUT_RAW'} : 'TRUE';# 是否输出原始数据sheet
    my $is_circos_plot  = (exists $hashConfig->{'PLOT_CIRCOS'}) ? $hashConfig->{'PLOT_CIRCOS'} : 'TRUE';# 是否绘制圈图
    my @samples         = package::utils::get_sample($hashConfig, "case", "control");
    my @cases           = package::utils::get_sample($hashConfig, "case");
    my $cnvSeqDir       = $hashConfig -> {'Report'}."/cnvseq";
    my $plotmark        = (package::utils::is_dir_ok($cnvSeqDir) == 1) ? 1 : 0;
    my $plotDirTemp     = $hashConfig -> {'Report'}."/sv/region_plot"; package::utils::make_dir($plotDirTemp) if ($plotmark ==1);  #region plot 临时文件存储目录，后续需要删除
    my %hashRegion;
    my %sample_project;
    map{$sample_project{$_} ++} @samples;
    my @DatabaseSamples = package::utils::GetSVDBSample($hashConfig, $hashPara, \%sample_project, $SVdbtype);
    my $OutputDir       = $hashConfig->{Report}."/document"; package::utils::make_dir($OutputDir);
    my $sv_dir          = "$OutputDir/4_SV"; package::utils::make_dir($sv_dir);
    my $excel           = "$sv_dir/SV.xlsx"; print "Output $excel ... ";
    my $workbook        = Excel::Writer::XLSX->new($excel); 
    my %format          = package::utils::formatRun($workbook);    
    my @sheetNames      = ("SV Filter","CNV Filter","SV Original","CNV Original");
    my %hashsheet       = (); #sheet
    my %hashrow         = (); #sheet的row
    my %hashTitleDatas  = (); #value color
    my %plotDatas       = (); #circos data
    map { $hashrow{$_} ++ } @sheetNames;
    my @SVtitles       = ("Chr 1","Position 1","Gene 1","Gene Region 1","Chr 2","Position 2","Gene 2","Gene Region 2",'Cover Region','Gene 1 HGMD','Gene 2 HGMD','Gene 1 MalaCards',
                          'Gene 2 MalaCards','Gene 1 OMIM','Gene 2 OMIM','Strands Info','ControlDBFreq','CaseFreq','ControlFreq','CaseSample', 'ControlSample',"Start Homology", "End Homology"
                          );
    my @CNVtitles      = ("Chr","Start","End","SVlength","SVtype","Gene","Gene Region",'HGMD','MalaCards','OMIM','ControlDBFreq','CaseFreq','ControlFreq','CaseSample', 'ControlSample',"Start Homology","End Homology","dbVar",
                          "cytoBand","dgvFreq","iscaPathGainCum","iscaPathLossCum","iscaLikelyPathogenic","iscaPathogenic","iscaCuratedPathogenic","CNVD","DECIPHER"
                         );
    foreach my $name (@sheetNames)
    {
           $hashsheet{$name} =  $workbook->add_worksheet("$name");
           $hashsheet{$name} -> set_row(0, 30);           
           $hashsheet{$name} -> set_column(1,6,15); 
           map { $hashTitleDatas{$name}{'data'} .= "$_!!!!!" ; $hashTitleDatas{$name}{'color'} .= "title!!!!!" } (@SVtitles, @samples) if ($name =~ /SV/);           
           map { $hashTitleDatas{$name}{'data'} .= "$_!!!!!" ; $hashTitleDatas{$name}{'color'} .= "title!!!!!" } (@CNVtitles, @samples) if ($name =~ /CNV/);           
    }
    sheetRawPrint (\%hashsheet, \%hashTitleDatas, \%hashrow, \%format);
    foreach my $title (sort {$hashSV->{$a}{'Order'} <=> $hashSV->{$b}{'Order'}} keys %$hashSV){ 
        my $FilterFlag = $hashSV -> {$title}{'FilterFlag'};
        my %hashdatas  = ();        
        my ($chr1, $pos1, $chr2, $pos2, $svtype) = split /\|/, $title;
        my $title2 = "$chr2|$pos2|$chr1|$pos1|$svtype";
        my @tempTitle = ($svtype eq "BND" or $svtype eq "INV") ?  @SVtitles : @CNVtitles;
        
        # 注释
        my ($value_string, $color_string) = ("", "");
        foreach (@tempTitle){
               my $value = " ";
               $value = exists ($hashSV->{$title}{$_})                     ? $hashSV->{$title}{$_}                     : " " if ($value eq " ");
               $value = exists ($hashAnnotation->{$title}{$_})             ? $hashAnnotation->{$title}{$_}             : " " if ($value eq " ");                   
               $value = exists ($hashAnnotation->{$title}{'Gene'})         ? $hashAnnotation->{$title}{'Gene'}         : " " if ($_ eq 'Gene 1');
               $value = exists ($hashAnnotation->{$title2}{'Gene'})        ? $hashAnnotation->{$title2}{'Gene'}        : " " if ($_ eq 'Gene 2');
               $value = exists ($hashAnnotation->{$title}{'Gene Region'})  ? $hashAnnotation->{$title}{'Gene Region'}  : " " if ($_ eq 'Gene Region 1');
               $value = exists ($hashAnnotation->{$title2}{'Gene Region'}) ? $hashAnnotation->{$title2}{'Gene Region'} : " " if ($_ eq 'Gene Region 2');
               $value = exists ($hashAnnotation->{$title}{'HGMD'})  ? $hashAnnotation->{$title}{'HGMD'}  : " " if ($_ eq 'Gene 1 HGMD');
               $value = exists ($hashAnnotation->{$title2}{'HGMD'}) ? $hashAnnotation->{$title2}{'HGMD'} : " " if ($_ eq 'Gene 2 HGMD');
               $value = exists ($hashAnnotation->{$title}{'MalaCards'})  ? $hashAnnotation->{$title}{'MalaCards'}  : " " if ($_ eq 'Gene 1 MalaCards');
               $value = exists ($hashAnnotation->{$title2}{'MalaCards'}) ? $hashAnnotation->{$title2}{'MalaCards'} : " " if ($_ eq 'Gene 2 MalaCards');               
               $value = exists ($hashAnnotation->{$title2}{'OMIM'}) ? $hashAnnotation->{$title2}{'OMIM'} : " " if ($_ eq 'Gene 1 OMIM');               
               $value = exists ($hashAnnotation->{$title2}{'OMIM'}) ? $hashAnnotation->{$title2}{'OMIM'} : " " if ($_ eq 'Gene 2 OMIM');               
               $value = exists ($hashAnnotation->{"$title|cover"}{'Gene Region'}) ? $hashAnnotation->{"$title|cover"}{'Gene Region'} : " " if ($_ eq 'Cover Region');
               $value_string .= "$value!!!!!";
               $color_string .= "normal!!!!!";
         }
         foreach (@samples)
         {   
             my $color;
             my $value = exists $hashSV ->{$title}{'SampleOutputInfo'}{$_}  ? $hashSV->{$title}{'SampleOutputInfo'}{$_} : " ";  
             my ($su, $cnv) = split /:/, $value;  
             $color = ( $su =~ /\d/ and $su >= $SUMin) ? 'orange' : 'normal';
             $color = ( $cnv =~ /\d/ and $cnv < 1.6 ) ? "skyblue" : 'normal' if ($svtype eq 'DEL' and (defined $cnv));
             $color = ( $cnv =~ /\d/ and $cnv > 2.4 ) ? "red"     : 'normal' if ($svtype eq 'DUP' and (defined $cnv));
             $value_string .= "$value!!!!!";
             $color_string .= "$color!!!!!";
             
             # circos 绘图数据保存
             $plotDatas{$_}{$svtype}{"$chr1\t$pos1\t$pos1\t$chr2\t$pos2\t$pos2"} = 1 if ($FilterFlag == 0 and ($svtype eq "BND" or $svtype eq "INV") and $su =~ /\d/ and $su >= $SUMin and $chr1 !~ /MT/i and $chr2 !~ /MT/i);            
             $plotDatas{$_}{$svtype}{"$chr1\t$pos1\t$pos2"}                      = 1 if ($FilterFlag == 0 and ($svtype ne "BND" and $svtype ne "INV") and $su =~ /\d/ and $su >= $SUMin and $chr1 !~ /MT/i );  
         } 

 
        # 结果保存 
        
        if ($svtype eq "BND" or $svtype eq "INV"){
            $hashdatas{'SV Original'}{'data'}   = $value_string if($is_output_raw eq 'TRUE');
            $hashdatas{'SV Original'}{'color'}  = $color_string if($is_output_raw eq 'TRUE');
            if ($FilterFlag == 0){
                $hashdatas{'SV Filter'}{'data'}     = $value_string;
                $hashdatas{'SV Filter'}{'color'}    = $color_string;
            }
        }else{        
            $hashdatas{'CNV Original'}{'data'}  = $value_string if($is_output_raw eq 'TRUE');
            $hashdatas{'CNV Original'}{'color'} = $color_string if($is_output_raw eq 'TRUE');
            if ($FilterFlag ==0){
                $hashdatas{'CNV Filter'}{'data'}    = $value_string;
                $hashdatas{'CNV Filter'}{'color'}   = $color_string;
                if ($plotmark == 1){  
                    my $length = $pos2 - $pos1 + 1;
                    next if ($chr1 =~ /\D/);   # 除去性染色体和非常规染色体
                    foreach my $sample (@samples) {
                         my ($su2, $pe2, $sr2, $cnv2) = split /:/, $hashSV->{$title}{'SampleSVInfo'}{$sample};
                         $cnv2  .= "_filtered_out" if($su2 !~ /\d/ or $cnv2 !~ /\d/ or $su2 <= 10 or ( ($svtype eq 'DEL' and $cnv2 >= 1.6) or ($svtype eq 'DUP' and $cnv2 <= 2.4) ) );                                   
                         my $startFinal = $pos1 - $length > 0 ? $pos1 - $length : 0;  ##区域左右扩展相同长度
                         my $endFinal   = $pos2 + $length;
                         $hashRegion{$sample} .= "$chr1\t$startFinal\t$endFinal\t$pos1\t$pos2\t$svtype\t$cnv2\n";
                    } 
               }                    
            }
        }
               
        # 结果输出
        sheetRawPrint (\%hashsheet, \%hashdatas, \%hashrow, \%format);     
    }                    
    writeReadme($workbook, \%format, "readmeWGS.txt") if ($SVdbtype eq "SVdbWGS");
    writeReadme($workbook, \%format, "readmeWES.txt") if ($SVdbtype eq "SVdbWES");
    print "ok\n";
    package::circos::run($hashConfig, $hashPara, \%plotDatas) if($is_circos_plot eq 'TRUE'); 
    regionPlotBedOutput($plotDirTemp, \@samples, \%hashRegion) if ($plotmark ==1);   
}

sub datafilter{
    my $hashSV          = shift @_;
    my $hashAnnotation  = shift @_;
    my $hashPara        = shift @_;
    my $hashConfig      = shift @_;
    my $SVdbtype        = shift @_;
    my $SUMin           = (exists($hashConfig->{'MinReads'}))   ? $hashConfig->{'MinReads'}   : 5;# 对变异结果的证据Reads数不得低于5重(>=5)
    my @cases           = package::utils::get_sample($hashConfig,"case");
    my @controls        = package::utils::get_sample($hashConfig,"control");
    my @samples         = (@cases,@controls);
    my %sample_project;
    map{$sample_project{$_} ++} @samples;
    my @DatabaseSamples = package::utils::GetSVDBSample($hashConfig, $hashPara, \%sample_project, $SVdbtype);   
    foreach my $title (keys %$hashSV)
    {
        my ($chr1, $pos1, $chr2, $pos2, $svtype) = split /\|/, $title;
        my $title2 = "$chr2|$pos2|$chr1|$pos1|$svtype";
        # 区域过滤标记 1表示需要过滤
        $hashSV->{$title}{'FilterFlag'} = 0 ;         
        # reads , CNV过滤条件 至少一个样本需满足条件 则保留此区域 只针对 WGS
        foreach my $sample(@samples) {
            my $svtype = $hashSV->{$title}{'SVtype'};
            my ($su, $pe, $sr, $cnv) = split /:/, $hashSV->{$title}{'SampleSVInfo'}{$sample};
            if (($SVdbtype eq "SVdbWGS") and $su =~ /^\d+$/ and $su > 10 and (($svtype eq 'DEL' and $cnv =~ /\d/ and $cnv < 1.6) or ($svtype eq 'DUP' and $cnv =~ /\d/ and $cnv > 2.4) or $svtype eq 'BND' or $svtype eq 'INV') )
            {
               $hashSV->{$title}{'FilterFlag'} = 0;
               last;   
            }
            elsif (($SVdbtype eq "SVdbWES") and $su =~ /^\d+$/ and $su > 10)
            {
               $hashSV->{$title}{'FilterFlag'} = 0;
               last; 
            }
            else
            {
               $hashSV->{$title}{'FilterFlag'} = 1;             
            }
        }     
        # 同源性
        my $start_hom                   = $hashAnnotation -> {$title}{'Start Homology'};
        my $end_hom                     = $hashAnnotation -> {$title}{'End Homology'};
        $hashSV->{$title}{'FilterFlag'} = 1 if ($start_hom > 10 and $end_hom > 10); 
      
        # # 功能性设定  priority
        my $geneRegion = $hashAnnotation -> {$title}{'Gene Region'};
        my $priority   = ($geneRegion eq 'exonic') ? 'First' : 'Second';# 外显子变异，造成功能性改变           
        if ($svtype =~ /BND/ or $svtype =~ /INV/){
            my $geneRegion2   = $hashAnnotation -> {$title2}{'Gene Region'};
            my $gene2         = $hashAnnotation -> {$title2}{'Gene'};
            my $gene1         = $hashAnnotation -> {$title}{'Gene'};
            my $coverRegion   = $hashAnnotation -> {"$title|cover"}{'Gene Region'};
            if($chr1 eq $chr2 and $coverRegion eq 'exonic' and $gene1 eq $gene2) # BND位于相同染色体相同基因,看覆盖区域内是否有外显子
            { 
                $priority='First';
            }
            elsif( ( ($chr1 ne $chr2) or ($chr1 eq $chr2 and $gene1 ne $gene2) ) and ($geneRegion ne 'intergenic' or $geneRegion2 ne 'intergenic')) # BND位于不同染色体，或者同一染色体不同基因，倒易位点不是同时在基因间区，则会造成功能性改变
            {
                $priority='First';
            }   
        }      
        # Control频率设定 ：对照数据库频率
        my @projectDepthOKSamples;# 项目样本中，结果测序深度合格的样本
        my @projectDepthOKCaseSamples;# 项目样本中，结果测序深度合格的case样本
        my @projectDepthOKControlSamples;# 项目样本中，结果测序深度合格的control样本
        my @DatabaseDepthOKSamples;# 对照数据库中，结果测序深度合格的样本
        foreach my $sample((@samples, @DatabaseSamples))
        {
            my $su = (exists $hashSV->{$title}{'SampleOutputInfo'}{$sample}) ? $hashSV->{$title}{'SampleOutputInfo'}{$sample} : 0;
               ($su) = (split /:/, $su)[0];
            my $is_depth_ok = ($su=~/\d/ and $su >= $SUMin) ? 1 : 0;
            next if($is_depth_ok == 0);
            push @projectDepthOKSamples, $sample        if($sample ~~ @samples);
            push @projectDepthOKCaseSamples, $sample    if($sample ~~ @cases);
            push @projectDepthOKControlSamples, $sample if($sample ~~ @controls);
            push @DatabaseDepthOKSamples, $sample       if($sample ~~ @DatabaseSamples);
        }
        my $DatabasePerc = (scalar(@DatabaseSamples) == 0) ? 0 : @DatabaseDepthOKSamples / scalar(@DatabaseSamples);# 对照数据库中频率
        my $casePerc     = (scalar(@cases) == 0)           ? 0 : @projectDepthOKCaseSamples / scalar(@cases);       # case样本中频率
        my $controlPerc  = (scalar(@controls) == 0)        ? 0 : @projectDepthOKControlSamples / scalar(@controls); # control样本中频率
        $hashSV          -> {$title}{'ControlDBFreq'} = $DatabasePerc;
        $hashSV          -> {$title}{'CaseFreq'}      = $casePerc;
        $hashSV          -> {$title}{'CaseSample'}    = join ",", @projectDepthOKCaseSamples;
        $hashSV          -> {$title}{'ControlFreq'}   = $controlPerc;
        $hashSV          -> {$title}{'ControlSample'} = join ",", @projectDepthOKControlSamples;
        $hashSV          -> {$title}{'FilterFlag'}    = 1 if($DatabasePerc > 0.1 or $priority ne 'First'); 
    }        
}  
         
sub annotation{
    my $hashSV             = shift @_;  
    my $hashPara           = shift @_;
    my $hashConfig         = shift @_;
    my %hashAnnotation     = ();
    my %hashAnnotationList = package::utils::get_annotation_list($hashPara, $hashConfig); # 获取当前物种的注释列表
    my %hashKeyWord        = package::annotation_db::get_annotation_key_word($hashPara, $hashConfig); # 自定义关键词
    foreach my $need_model(sort keys %hashAnnotationList)
    {   
        my $method_get_result = $hashAnnotationList{$need_model}{'GetResult'};
        package::utils::read_general_snp_annotation(\%hashAnnotation, $hashAnnotationList{$need_model}, \%hashKeyWord, $need_model) if($method_get_result eq 'General Pos'); # 提取常规SNP注释
    } 
    my $species      = package::utils::get_species($hashConfig);    
    my $annovarBuild = $hashPara -> {$species}{'AnnovarBuild'};    

   # 同源性注释信息提取
    my $HomologyFile = $hashConfig -> {'Report'}."/sv/homology/homology.anno";
    Homology ($hashSV, \%hashAnnotation, $HomologyFile);  


   # 基因疾病注释信息提取
   foreach my $need_model(sort keys %hashAnnotationList)    
   { 
       my $method_get_result = $hashAnnotationList{$need_model}{'GetResult'};    
       next if($method_get_result ne 'General Gene');       
       my %temp              = package::utils::read_general_gene_annotation($hashAnnotationList{$need_model}, \%hashKeyWord, $need_model, $species); # 提取gene注释
       foreach my $snptitle (keys %hashAnnotation)
       {
           my $gene = $hashAnnotation{$snptitle}{'Gene'} if (exists $hashAnnotation{$snptitle}{'Gene'});
           foreach (split /,/, $gene)
           {
               next if (not exists $temp{$_});
               foreach my $info_name(keys %{$temp{$_}})
               {
                   my $value = $temp{$_}{$info_name}; 
                   $hashAnnotation{$snptitle}{"$info_name"} .="$_:[$value]; ";                              
               }       
           }                
       } 
   }   
   return %hashAnnotation;        
}

sub Homology{
    my $hashSV         = shift @_;
    my $hashAnnotation = shift @_;
    my $HomologyFile   = shift @_;
    my %hits           = ();
    open HITS, $HomologyFile;
    while (<HITS>)
    {
          $_ =~ s/[\r\n]//g;
          my ($chr, $pos, $count) =split /\t/, $_;
          $hits {"$chr|$pos"} = $count;
    }
    close HITS;
    foreach my $title (keys %$hashSV)
    {                   
            my ($chr1, $pos1, $chr2, $pos2, $svtype) = split /\|/, $title;
            my $hits1 = exists $hits{"$chr1|$pos1"} ? $hits{"$chr1|$pos1"} : 0;
            my $hits2 = exists $hits{"$chr2|$pos2"} ? $hits{"$chr2|$pos2"} : 0;
            $hashAnnotation -> {$title}{'Start Homology'} = $hits1;
            $hashAnnotation -> {$title}{'End Homology'}   = $hits2;
    }
}

sub getGeneAnno{
    my $hashAnnotation = shift @_;
    my $temp           = shift @_;

}    

sub sheetRawPrint{
    my $hashsheet   = shift @_;
    my $hashdatas   = shift @_;
    my $hashrow     = shift @_;
    my $format      = shift @_;
    foreach my $name (keys %$hashdatas)
    {
        my $datas  =  $hashdatas -> {$name}{'data'} ;
        my $colors =  $hashdatas -> {$name}{'color'};        
        my @data  = split /!!!!!/, $datas;
        my @color = split /!!!!!/, $colors;
        my $row   = ($hashrow -> {$name}) - 1;
        foreach my $col  (0..$#data)
        {
          ($hashsheet -> {$name}) -> write($row, $col, $data[$col], $format -> {$color[$col]});            
        }
        $hashrow -> {$name} ++ ;
    }
}

sub writeReadme{
    my $workbook = shift @_;
    my $format   = shift @_;
    my $file     = shift @_;
    open IN, "$file";
    my $readme = $workbook -> add_worksheet("ReadMe");
    my $row=0;
    $readme -> set_column(0,0,15);
    $readme -> set_column(1,1,50);
    while(<IN>){
        $_=~s/[\r\n]//g;
        my @data=split /\t/,$_;
        $readme -> write($row,0,decode("gb2312", $data[0]), $format->{"title"});
        $readme -> write($row,1,decode("gb2312", $data[1]), $format->{"normal"});
        $row++;     
    }
    close IN;
}

sub regionPlotBedOutput{
    my $plotDirTemp = shift @_;
    my $samples     = shift @_;
    my $hashRegion  = shift @_;
    foreach my $sample (@$samples)
    {
        my $regionBed = "$plotDirTemp/$sample\_region.bed";
        open BED, ">$regionBed";
        my $info = $hashRegion -> {$sample};
        print BED "$info\n";
    }
    close BED; 
}


1