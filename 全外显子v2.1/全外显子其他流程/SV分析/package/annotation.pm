package package::annotation;
use strict;
use warnings;

# 获取注释关键字
sub run{
    my $hashPara    = shift @_;
    my $hashConfig  = shift @_;
    my $report_dir  = $hashConfig->{'Report'};
    my $SVdbtype    = (exists $hashConfig->{'GenomeBed'}) ? "SVdbWGS" : "SVdbWES"; 
    my $SVreportDir = "$report_dir/sv";
    my $library_dir = "$SVreportDir/library"; package::utils::make_dir($library_dir);
    my $finishFlag  =  "$library_dir/FinishFlag";
    print "########## Start Annotation ".package::utils::get_time()." ##########\n";
    if(package::utils::is_file_ok($finishFlag))
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    } 
    my $species            = package::utils::get_species($hashConfig);
    my $AnnovarDir         = $hashPara->{'Soft'}{'AnnovarDIR'};
    my $AnnovarBuild       = $hashPara->{$species}{'AnnovarBuild'};
    my $AnnovarDB          = $hashPara->{$species}{'AnnovarDB'};
    my $genome             = $hashPara->{$species}{'Genome'};
    my %hashAnnotationList = package::utils::get_annotation_list($hashPara, $hashConfig); # 获取当前物种的注释列表
    my %hashKeyWord        = package::annotation_db::get_annotation_key_word($hashPara, $hashConfig); # 自定义关键词

   #####
   # 开始注释
   #####
    
    # 1 建立注释库文件
    print "1.Build annovar input\n";
    my $library_file      = "$library_dir/library"; 
    my $svFile            = "$SVreportDir/sample.all.final.sv.vcf";
    my %hashSV            =  package::utils::readSV($svFile, $hashConfig, $SVdbtype);
    create_annovar_input(\%hashSV, $library_file);  # 创建annovar输入文件
    
    # 2 注释
    print "2.annovar annotation\n";
    my $threshold = exists $hashPara->{"Process"}{"Annotation"} ? $hashPara->{"Process"}{"Annotation"} : 10;
       $threshold = $hashConfig->{"Process_Annotation"} if(exists $hashConfig->{"Process_Annotation"});
    my $pm = Parallel::ForkManager->new($threshold);
    foreach my $anno_model(keys %hashAnnotationList)
    {
        next if($hashAnnotationList{$anno_model}{'From'} ne 'Annovar');
        $pm->start() and next;
        my $AnnovarDB = $hashAnnotationList{$anno_model}{'DBDir'}{$species};
        system("$AnnovarDir/$hashAnnotationList{$anno_model}{'Para'} --buildver $AnnovarBuild $library_file $AnnovarDB");
        process_dgvFreq($library_file, $AnnovarBuild) if ($anno_model eq 'dgvFreq'); # 重新生成 dgvFreq 注释文件
        process_dbVar($library_file, $AnnovarBuild)   if ($anno_model eq 'dbVar');   # 重新生成 dbVar 注释文件
        $pm->finish; 
    }
    $pm->wait_all_children;
     
    # 完成
    open OUT, ">$finishFlag";
    print OUT "Annotation Finish\n";
    close OUT;
}

sub create_annovar_input{
    my $hashSV     = shift @_;
    my $library    = shift @_;
    my %titleCount = ();
    open LIBRARY,">$library";
    foreach my $title(sort keys %$hashSV){
        my $svType     = $hashSV->{$title}{'SVtype'};
        my ($chr1, $pos1, $chr2, $pos2, $svtype)=split /\|/,$title;
        if ( $svtype =~ /BND/ or $svtype =~ /INV/)
        {
            my $title2 = "$chr2|$pos2|$chr1|$pos1|$svtype";       
            my $title3 = "$chr1|$pos1|$chr2|$pos2|$svtype|cover";       
            print LIBRARY "$chr1\t$pos1\t$pos1\t0\t-\t$title\n";
            print LIBRARY "$chr2\t$pos2\t$pos2\t0\t-\t$title2\n";
            print LIBRARY "$chr1\t$pos1\t$pos2\t0\t-\t$title3\n" if ($chr1 eq $chr2);
        }
        else
        {
            print LIBRARY "$chr1\t$pos1\t$pos2\t0\t-\t$title\n";
        }
    }
    close LIBRARY; 
}

sub process_dgvFreq{
    my $library_file = shift @_;
    my $AnnovarBuild = shift @_;
    my %hashAnno     = ();
    my $annoFile     = "$library_file.$AnnovarBuild\_dgvFreq";
    open DGVFREQ, $annoFile;
    while (<DGVFREQ>)
    {
          $_=~ s/[\r\n]//g;       
          my ($database, $DBannoInfo, $chr, $pos1, $pos2, $ref, $alt, $title) = split /\t/,$_;
          my $svtype = (split /\|/, $title)[-1];
          next if ($svtype eq "BND" or $svtype eq "INV");
          $DBannoInfo   =~ s/^Name\=//; 
          my ($gainAdd, $lossAdd, $allAdd) = (0, 0, 0);
          foreach my $key(split /,/, $DBannoInfo)
          {
                  my ($gain, $loss, $all) = split /\;/, $key;
                  $gainAdd += $gain;
                  $lossAdd += $loss;
                  $allAdd  += $all;
         }
         my $value = ($allAdd > 0) ? (($gainAdd + $lossAdd) / $allAdd) : 0;
         $hashAnno{$title} = $value;          
    }
    close DGVFREQ;    
  
    open OUT, ">$annoFile";
    foreach my $title (sort keys %hashAnno)
    {
        my $value = $hashAnno{$title};
        print OUT "$title\t$value\n";
    }
    close OUT;
}


sub process_dbVar{
    my $library_file = shift @_;
    my $AnnovarBuild = shift @_;
    my %hashAnno     = ();
    my $annoFile     = "$library_file.$AnnovarBuild\_dbVar";
    open DBVAR, $annoFile;
    while (<DBVAR>)
    {
          $_=~s/[\r\n]//g;
          my ($database, $DBannoInfo, $chr, $pos1, $pos2, $ref, $alt, $title) = split /\t/,$_;
          my $svtype = (split /\|/, $title)[-1];
          next if ($svtype eq "BND" or $svtype eq "INV");
          $DBannoInfo   =~ s/^Name\=//;
          my @annoInfos = split /[,;]/,$DBannoInfo;
          my %hashCNV   = ();
          foreach my $annoInfo (@annoInfos)
          {
              next if($annoInfo!~/\d/);
            my ($CNVType,$count) = split /\=/, $annoInfo;
            $hashCNV{$CNVType}  += $count;
          }
          my $CNVstring    = "";
          map{ $CNVstring .= "$_=$hashCNV{$_};" }sort keys %hashCNV;
          $hashAnno{$title} = $CNVstring;        
    }
    close DBVAR;
    
    open OUT, ">$annoFile";
    foreach my $title (sort keys %hashAnno)
    {
        my $value = $hashAnno{$title};
        print OUT "$title\t$value\n";
    }
    close OUT;
}







1
