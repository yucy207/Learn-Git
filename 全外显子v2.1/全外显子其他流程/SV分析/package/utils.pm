package package::utils;
use strict;
use warnings;
$|=1;

# 读取配置文件
sub read_config{
    my $config_file = shift @_;
    my %hashConfig;
    open CONFIG, $config_file or die "[ERR] Loss Config File,$config_file\n";
    while(my $_=<CONFIG>){
        $_=~s/[\s\t\r\n]//g;
        next if($_ =~/^#/ or $_!~/\w/);
        my ($name, $value) = split /=/, $_;
        next if(!defined $value);

        if($name eq 'case' or $name eq 'control')
        {
            $hashConfig{$name} .= "$value,";
        }
        elsif($name eq 'Force_Run' and $value =~ /\w/) 
        {
            map{ $hashConfig{'Force_Run'}{$_}++; } split /,/, $value;
        }
        else
        {
            $hashConfig{$name} = $value;
        }
    }
    close CONFIG;
    return %hashConfig;
}

# 获取默认参数
sub get_para{
    my %hashPara;

    # 软件配置
    $hashPara{'Soft'}{'Tmp'}        = "/home/tmp";
    $hashPara{'Soft'}{'BWA'}        = "/home/genesky/software/bwa/0.7.17/bwa";
    $hashPara{'Soft'}{'SamTools'}   = "/home/genesky/software/samtools/1.9/samtools";
    $hashPara{'Soft'}{'samblaster'} = "/home/genesky/software/samblaster/0.1.24/samblaster";
    $hashPara{'Soft'}{'sambamba'}   = "/home/genesky/software/sambamba/0.6.7/sambamba";
    $hashPara{'Soft'}{'PicardDIR'}  = "/home/genesky/software/picard/2.18.29/picard.jar";
    $hashPara{'Soft'}{'bedtools'}   = "/home/genesky/software/bedtools/2.28.0/bin/bedtools";
    $hashPara{'Soft'}{'svtools'}    = "/home/genesky/software/python/2.7.13/bin/svtools"; # 0.3.2版本，0.4.0有bug，lmerge之后，vcf头多了VARIOUS字符
    $hashPara{'Soft'}{'AnnovarDIR'} = "/home/genesky/software/annovar/2018Apr16";
    $hashPara{'Soft'}{'BlastDIR'}   = "/home/genesky/software/blast/20110130/bin";
    $hashPara{'Soft'}{'Java'}       = "/home/genesky/software/java/1.8.0_181/bin/java";
    $hashPara{'Soft'}{'speedseq'}   = "/home/genesky/software/speedseq/0.1.2/bin";
    $hashPara{'Soft'}{'Circos'}     = "/home/genesky/software/circos/0.69-6/bin/circos";
    $hashPara{'Soft'}{'InterVar'}   = "/home/genesky/software/intervar/2.1.2";
    $hashPara{'Soft'}{'Python27'}   = "/home/genesky/software/python/2.7.13/bin/python";
    $hashPara{'Soft'}{'R'}          = "/home/genesky/software/r/3.5.1/bin/R";
    $hashPara{'Soft'}{'RLib'}       = "/home/genesky/software/r/3.5.1/lib64/R/library";

    # 参考基因组配置

    $hashPara{'Human_hg19'}{'Genome'}         = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.fa";
    $hashPara{'Human_hg19'}{'blast_idx'}      = "/home/genesky/database/ucsc/hg19_modify/genome/blast+_idx/hg19_modify.fa";  
    $hashPara{'Human_hg19'}{'DBsnp'}          = "/home/genesky/database/gatk/hg19_modify/dbsnp/dbsnp_138.hg19.vcf";
    $hashPara{'Human_hg19'}{'InDel'}          = "/home/genesky/database/gatk/hg19_modify/1000G_phase1_indels/1000G_phase1.indels.hg19.vcf";
    $hashPara{'Human_hg19'}{'InDelGold'}      = "/home/genesky/database/gatk/hg19_modify/mills_and_1000g_gold_standard_indels/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf";
    $hashPara{'Human_hg19'}{'AnnovarBuild'}   = "hg19";
    $hashPara{'Human_hg19'}{'Excludebed'}     = "/home/genesky/software/speedseq/0.1.2/annotations/ceph18.b37.exclude.2014-01-15.bed";
    $hashPara{'Human_hg19'}{'SVdbWGS'}        = "/home/genesky/database/self_build_database/SV_VCF_DB_genome";
    $hashPara{'Human_hg19'}{'SVdbWES'}        = "/home/genesky/database/self_build_database/SV_VCF_DB_exon/hg19";
    $hashPara{'Human_hg19'}{'SVCircosKT'}     = "/home/genesky/database/self_build_database/circos-karyotype/karyotype.human.hg19.rmMT.txt";
    
    return %hashPara;
}

# 读取SV文件
sub readSV{
    my $SVFile     = shift @_;
    my $hashConfig = shift @_;
    my $SVdbtype   = shift @_;
    my %hashPara   = get_para();
    my $species    = get_species($hashConfig);
    my $dict       = (split /\./, $hashPara{$species}{'Genome'})[0].".dict";
    my %orders     = readdict($dict); #染色体顺序
    my $SUMin      = (exists($hashConfig -> {'MinReads'}))   ? $hashConfig -> {'MinReads'}   : 5;# 对变异结果的证据Reads数不得低于5重(>=5)
    my @samples    = get_sample($hashConfig, "case","control");
    my $rmChrWith_ = (exists($hashConfig -> {'rmChrWith_'}) and $hashConfig -> {'rmChrWith_'} eq 'TRUE') ? 1 : 0;# 是否排除含有_的染色体
    print "Reading $SVFile ...  ";
    open SV, $SVFile;
    my @heads;
    my @VCFsamples;
    my %hashSV = ();
    my %titles = ();
    my %dups   = (); #判断重复的位置区域
    ##
    while (<SV>){
        $_ =~ s/[\r\n]//g;
        my @datas = split /\t/, $_;
        if ($_ =~ /^#CHROM/){
           @heads      = @datas;
           $heads[0]   =~ s/^#//;
           @VCFsamples = splice @datas, 9;
        }
        next if ($_ =~ /^#/);
        my %hashInfos = (); 
        map { $hashInfos{$heads[$_]} = $datas[$_] } (0..$#heads);
        readAnnoINFO(\%hashInfos);
        next if ($hashInfos{'CHROM'} =~ /\_/ and $rmChrWith_ == 1);
        # 对本项目样本检测覆盖是否合格
        my $nextSU = 0;##SU限制,只要有一个样本满足条件，则输出
        foreach my $sample (@samples){
            my ($gt, $su, $pe, $sr)=split /:/,$hashInfos{$sample};
            if($su =~ /\d/ and $su >= $SUMin){
               $nextSU = 1;
               last;
            }
        }
        next if ($nextSU == 0);
        my $chr      = $hashInfos{'CHROM'};
        my $pos      = $hashInfos{'POS'};
        my $ID       = $hashInfos{'ID'};
        my $svtype   = $hashInfos{'SVTYPE'};
        my ($chr1, $pos1, $chr2, $pos2) = ($chr, $pos, $chr, $hashInfos{'END'});  # DUP DEL INV
        if ($svtype eq "BND"){
            my ($empty1, $BNDchr, $BNDposition, $empty2)=split /[:\[\]]/,$hashInfos{'ALT'}; 
            ($pos1, $pos2)  = ($pos, $BNDposition) ;            
            ($chr1, $chr2) = ($chr, $BNDchr) if ($chr ne $BNDchr);        
        }
        ($pos1, $pos2) = sort {$a<=>$b} ($pos1, $pos2) if($chr1 eq $chr2);
        next if (($chr1 =~ /\_/ or $chr2 =~ /\_/) and $rmChrWith_ == 1);
        # 跳过的重复的记录，默认保留首条记录
        my $duptitle = join "|", sort ("$chr1|$pos1", "$chr2|$pos2");
        next if (exists $dups{$duptitle});
        $dups{$duptitle} = 1;
        #        
        my $title = "$chr1|$pos1|$chr2|$pos2|$svtype";
        $titles{$chr1}{$pos1}{$chr2}{$pos2}{$svtype} = $title;  # # 若键值不包含 $svtype，则同一位置不同类型的信息会被覆盖掉，导致一些title没有排序count      
        $hashSV{$title}{'Chr 1'}     = $chr1;
        $hashSV{$title}{'Chr 2'}     = $chr2;
        $hashSV{$title}{'Position 1'}= $pos1;
        $hashSV{$title}{'Position 2'}= $pos2;
        $hashSV{$title}{'Chr'}       = $chr;
        $hashSV{$title}{'Start'}     = $pos;
        $hashSV{$title}{'End'}       = $pos2;
        $hashSV{$title}{'SVlength'}  = $hashInfos{'SVLEN'};
        $hashSV{$title}{'SVtype'}    = $svtype;
        my $infos  = $hashInfos{'ALT'};
        my $posnew = "$chr1:$pos1";
        $infos =~ s/N/$posnew/;
        $hashSV{$title}{'Strands Info'}   = $infos if ($svtype eq "BND");
        foreach my $sample(@VCFsamples){
            my ($gt, $su, $pe, $sr, $cnvOriginal)=split /:/,$hashInfos{$sample};
            my @tmp = split /:/,$hashInfos{$sample};
            my $cnv = (defined($cnvOriginal)) ? $cnvOriginal : "#";
            $hashSV{$title}{'SampleSVInfo'}{$sample}     = "$su:$pe:$sr:$cnv";
            $hashSV{$title}{'SampleOutputInfo'}{$sample} = ($SVdbtype eq "SVdbWGS") ? "$su:$cnv" : "$su"; 
            $hashSV{$title}{'SampleOutputInfo'}{$sample} = $su if ($SVdbtype eq "SVdbWES" or ($hashInfos{'SVTYPE'} eq "BND" or $hashInfos{'SVTYPE'} eq "INV"));
        }
    }
    close SV;
    ordersTitle(\%hashSV, \%titles, \%orders); # 排序    
    print "OK\n";
    return %hashSV;
}

sub readAnnoINFO{
    my $hashInfos = shift @_;
    my $annoInfo  = $hashInfos -> {'INFO'};
    my @datas     = split /;/,$annoInfo;
    foreach my $anno (@datas)
    {
        my ($content, $value)    = split /=/, $anno;
        $hashInfos -> {$content} = $value;
    }
}

# 获取注释列表
sub get_annotation_list{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    my $species     = package::utils::get_species($hashConfig);
    my $report_dir  = $hashConfig->{'Report'}; 
    my $library_dir = "$report_dir/sv/library";
    my $anno_log    = "$library_dir/anno.log";
    my %hashKeyWord = package::annotation_db::get_annotation_key_word($hashPara, $hashConfig); # 获取关键字哈希

    open ANNOLOG, ">$anno_log";
    my %hashAnnotation_all = package::annotation_db::get_annotation(); # 所有注释数据库
    my %hashAnnotation; # 当前物种最终可以注释的数据库
    foreach my $anno_model(sort keys %hashAnnotation_all)
    {
        my $model_species = $hashAnnotation_all{$anno_model}{'Species'};

        # 1 物种检测
        next if($model_species ne '[all]' and $model_species !~/\[$species\]/); # 过滤，非当前物种

        # 2 分析必要文件检测
        my $result = check_annotation_necessary_file($hashAnnotation_all{$anno_model}, \%hashKeyWord, $species);

        # 结论
        if($result eq '') # 正常
        {
            print ANNOLOG "annotation $anno_model is good to run\n";
            $hashAnnotation{$anno_model} = $hashAnnotation_all{$anno_model};
        }
        else # 有错误
        {
            print ANNOLOG "annotation $anno_model is Error: $result\n";
        }
    }
    return %hashAnnotation;
}

# 读取通用SNP注释
sub read_general_snp_annotation{
    my $hashAnnotation     = shift @_;
    my $hashAnnotationList = shift @_;
    my $hashKeyWord        = shift @_;
    my $need_model         = shift @_;

    foreach my $out_index(keys %{$hashAnnotationList->{'Output'}})
    {   
        my %hashFile_Info = split_annotation_info($hashAnnotationList->{'Output'}{$out_index}{'File'}); # 文件信息
        my %hashNeed_Info = split_general_snp_need_info($hashAnnotationList->{'Output'}{$out_index}{'Need_Info'}); # 需要的信息坐标
        # 文件检测
        my $prefix = (exists($hashKeyWord->{$hashFile_Info{'prefix'}})) ? $hashKeyWord->{$hashFile_Info{'prefix'}} : $hashFile_Info{'prefix'};
        my $suffix = $hashFile_Info{'suffix'};
        my $output_file = $prefix.$suffix; # 当前注释内容的输出文件
        if(is_file_ok($output_file)==0)
        {
            print "Lost $need_model : $output_file\n";
            next;
        }

        # 文件信息读取
        print "Process $need_model : $output_file ... ";
        open OUTPUTFILE, $output_file;
        while(<OUTPUTFILE>)
        {
            $_=~s/[\r\n]//g;
            next if($_=~/^#/);
            my @datas = split /\t/, $_;
            foreach my $info_name(keys %hashNeed_Info)
            {
                my $snp_title_index = $hashNeed_Info{$info_name}{'snp_title_index'}; 
                my $value_index     = $hashNeed_Info{$info_name}{'value_index'}; 
                my $value           = '';
                my $snp_title       = '';
                my @snp_titles      = ();
                foreach my $index (split /,/, $snp_title_index)
                {
                     my $st = exists($datas[$index]) ? $datas[$index] : '';
                     push @snp_titles, $st;
                }
                $snp_title = join "|", @snp_titles;
                
                if($value_index=~/\//) # 复合式结果，需要拆分提取
                {
                    my ($index1, $split_mark, $index2) = split /\//, $value_index;
                    my $value1 = exists($datas[$index1]) ? $datas[$index1] : '' ;
                    my @datas2 = split /$split_mark/, $value1;
                    $value = exists($datas2[$index2]) ? $datas2[$index2] : '';
                }
                else
                {
                    $value = exists($datas[$value_index]) ? $datas[$value_index] : '' ;
                }
                if ($info_name eq "Gene") #对UTS2(dist=9486),类型基因名进行特殊保留
                {
                    my @temp = ();
                    foreach (split /\)\,/, $value)
                    {
                        $_=~ s/\(.*//;
                        next if ($_ eq "NONE");
                        push @temp, $_;
                    }
                    $value = join ",", @temp; 
                    $value = "NONE" if ($value eq "");                    
                }
                $value =~ s/Name=// if ($info_name eq "iscaPathGainCum" or $info_name eq "iscaPathLossCum");
                $value =~ s/\\x2c/,/g; # annovar 注释内容会把‘，’或‘#’替换为特殊字符，防止冲突，所以使用时需要再替换回来
                $value =~ s/\\x23/#/g;

                $hashAnnotation->{$snp_title}{$info_name} = $value if($value  ne "");
            }
        }
        print "OK\n";
    }
}

# 读取通用GENE注释
sub read_general_gene_annotation{
    my $hashAnnotationList = shift @_;
    my $hashKeyWord        = shift @_;
    my $need_model         = shift @_;    
    my $species            = shift @_;
    my %temp = ();
    foreach my $out_index(keys %{$hashAnnotationList->{'NecessaryFile'}})
    {
        my %hashFile_Info = split_annotation_info($hashAnnotationList->{'NecessaryFile'}{$out_index}); # 文件信息
        my %hashNeed_Info = split_general_gene_need_info($hashAnnotationList->{'Output'}{$out_index}{'Need_Info'}); # 需要的信息坐标

        # 文件检测
        my $anno_dir = (exists $hashAnnotationList->{'DBDir'}{$species}) ? $hashAnnotationList->{'DBDir'}{$species}."/": "";
        my $prefix   = (exists($hashKeyWord->{$hashFile_Info{'prefix'}})) ? $hashKeyWord->{$hashFile_Info{'prefix'}} : $hashFile_Info{'prefix'};
        my $suffix   = $hashFile_Info{'suffix'};
        my $output_file = $anno_dir.$prefix.$suffix; # 当前注释内容的输出文件
        if(is_file_ok($output_file)==0)
        {
            print "Lost $need_model : $output_file\n";
            next;
        }

        # 文件信息读取
        print "Process $need_model : $output_file ... ";
        open OUTPUTFILE, $output_file;
        while(<OUTPUTFILE>)
        {
            $_=~s/[\r\n]//g;
            next if($_=~/^#/);
            my @datas = split /\t/, $_;
            foreach my $info_name(keys %hashNeed_Info)
            {
                my $gene        = $datas[$hashNeed_Info{$info_name}{'gene_index'}];
                my $value_index = $hashNeed_Info{$info_name}{'value_index'}; # 取值index
                my $value       = '';

                if($value_index=~/\//) # 复合式结果，需要拆分提取
                {
                    my ($index1, $split_mark, $index2) = split /\//, $value_index;
                    my $value1 = exists($datas[$index1]) ? $datas[$index1] : '' ;
                    my @datas2 = split /$split_mark/, $value1;

                    my @values = map{ my $value_tmp = exists($datas2[$_]) ? $datas2[$_] : '' ; $value_tmp} (split /,/, $index2); # 应对需要获得多个结果的情况
                    $value = join "|", @values;
                }
                else
                {   
                    my @values = map{ my $value_tmp = exists($datas[$_]) ? $datas[$_] : '' ; $value_tmp} (split /,/, $value_index); # 应对需要获得多个结果的情况
                    $value = join "|", @values;
                }
                $temp{$gene}{$info_name} = $value if($value ne "");
            }

        }
        print "OK\n";
    }
    return %temp; 
}

# 拆分注释信息
sub split_annotation_info{
    my $info = shift @_;
    my %hashInfo;
    foreach(split /\t/, $info)
    {
        my ($name, $value) = split /=/, $_;
        $hashInfo{$name} = $value;
    }
    return %hashInfo;
}

sub split_general_snp_need_info{
    my $need_info = shift @_;
    my %hashNeed_Info;
    foreach my $info_detail(split /\t/, $need_info)
    {
        my ($info_name, $value_index, $snp_title_index) = split /\=/, $info_detail;
        $hashNeed_Info{$info_name}{'value_index'}       = $value_index;
        $hashNeed_Info{$info_name}{'snp_title_index'}   = $snp_title_index;
    }
    return %hashNeed_Info;
}

sub split_general_gene_need_info{
    my $need_info = shift @_;
    my %hashNeed_Info;
    foreach my $info_detail(split /\t/, $need_info)
    {
        my ($info_name, $value_index, $gene_index) = split /\=/, $info_detail;
        $hashNeed_Info{$info_name}{'value_index'} = $value_index;
        $hashNeed_Info{$info_name}{'gene_index'}   = $gene_index;
    }
    return %hashNeed_Info;
}

# 检测必要文件
sub check_annotation_necessary_file{
    my $hashAnnotation = shift @_;
    my $hashKeyWord    = shift @_;
    my $species        = shift @_;

    my $isOK = 1;
    my @errors = ();
    foreach my $file_index(sort keys %{$hashAnnotation->{'NecessaryFile'}})
    {
        my $anno_dir = (exists $hashAnnotation->{'DBDir'}{$species}) ? $hashAnnotation->{'DBDir'}{$species}."/": "";
        my %hashInfo = split_annotation_info($hashAnnotation->{'NecessaryFile'}{$file_index});
        my $prefix = (exists($hashKeyWord->{$hashInfo{'prefix'}})) ? $hashKeyWord->{$hashInfo{'prefix'}} : $hashInfo{'prefix'};
        my $suffix = $hashInfo{'suffix'};
        my $need_file = $anno_dir.$prefix.$suffix;
        my $condition = (package::utils::is_file_ok($need_file)==1) ? 'OK' : 'Lost';
        next if($condition eq 'OK');
        push @errors, "Lost $need_file";
    }
    
    return (join ",", @errors);
}

# 获得control 数据库样本
sub GetSVDBSample{
    my $hashConfig     = shift @_;
    my $para           = shift @_;
    my $sample_project = shift @_;
    my $SVdbtype       = shift @_;
    my $species        = get_species($hashConfig);
    my @DatabaseSamples;
    my $SVDatabase_dir = (exists($para -> {$species}{$SVdbtype})) ? $para->{$species}{$SVdbtype} : "";
    return @DatabaseSamples if($SVDatabase_dir !~ /\w/);
    my @SVDatabase_vcf = glob "$SVDatabase_dir/*.sv.vcf";
    foreach my $vcf(@SVDatabase_vcf) {
        $vcf =~ /.*\/(.*).sv.vcf/;
        my $sample = $1;
        next if(exists $sample_project->{$sample}); # # 与项目有重叠的样本不做
        next if(is_file_ok($vcf) == 0);
        push @DatabaseSamples, $sample;
    }
    return @DatabaseSamples;
}

sub GetSVDBFiles{
    my $hashConfig     = shift @_;
    my $para           = shift @_;
    my $sample_project = shift @_;
    my $SVdbtype       = shift @_;
    my $species        = get_species($hashConfig);
    my %hashSVDB;
    my $SVDatabase_dir     = (exists($para->{$species}{$SVdbtype})) ? $para->{$species}{$SVdbtype} : "";
    return %hashSVDB if($SVDatabase_dir!~/\w/);
    my @SVDatabase_vcf = glob "$SVDatabase_dir/*.sv.vcf";
    my @sv_vcfs;
    foreach my $vcf(@SVDatabase_vcf) {
        $vcf =~ /.*\/(.*).sv.vcf/;
        my $sample = $1;
        next if(exists $sample_project->{$sample}); # # 与项目有重叠的样本不做
        next if(is_file_ok($vcf) == 0);
        push @sv_vcfs, $vcf;
    }
    $hashSVDB{'sv_vcfs'}=join ",", @sv_vcfs;
    return %hashSVDB;
}

# 获取当前项目分析的样本
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

# 读取参照基因组
sub readgenome{
    my $genome     = shift @_;
    my %hashGenome = ();
    my $chr        = "";
    open FILE,$genome;
    while(my $line = <FILE>){
        $line =~ s/[\r\n]//g;
        if($line =~ /^\>/){
            my @tmps = split /\s/, $line;
            $chr     = $tmps[0];
            $chr     =~ s/\>//g;
            print "\rReading $genome ... $chr ";
            next;
        }
        $hashGenome{$chr} .= $line if($chr ne "");
    }
    close FILE;
    print "\rReading $genome ... OK         \n";
    return %hashGenome;
}

# title 排序
sub ordersTitle{
    my $hashSV = shift @_;
    my $titles = shift @_;
    my $orders = shift @_;
    my $count = 0;
    foreach my $chr1 (sort {$orders -> {$a} <=> $orders -> {$b}} keys %$titles)
    {
        foreach my $pos1 (sort {$a <=> $b} keys %{$titles -> {$chr1}})
        {
            foreach my $chr2 (sort {$orders -> {$a} <=> $orders -> {$b}} keys %{$titles -> {$chr1}{$pos1}})
            {
                foreach my $pos2 (sort {$a <=> $b} keys %{$titles -> {$chr1}{$pos1}{$chr2}})
                {
                    foreach my $svtype (keys %{$titles -> {$chr1}{$pos1}{$chr2}{$pos2}})
                    {
                        $count ++;
                        my $title = $titles -> {$chr1}{$pos1}{$chr2}{$pos2}{$svtype};
                        $hashSV -> {$title}{'Order'} = $count;
                    }
                }            
            }            
        
        }    
    }
   
}

# 读取参照基因组dict信息，将title排序
sub readdict{
    my $file   = shift @_;
    my %orders = ();
    open IN, $file;
    my $count = 0;
    while (<IN>)
    {
          next if ($_ !~ /^\@/);
          $count ++;        
          my $chr = (split /\t/, $_)[1];
          $chr =~ s/SN\://g;
          $orders {$chr} = $count;
    }
    close IN;
    return %orders;
}

# 获取时间
sub get_time{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = localtime(time());
    $year += 1900;
    $mon  += 1;
    return "$year-$mon-$day $hour:$min:$sec";
}

sub show_config{
    my $hashConfig = shift @_;
    my @case_samples    = get_sample($hashConfig, "case");
    my @control_samples = get_sample($hashConfig, "control");

    my $case_num    = @case_samples;
    my $control_num = @control_samples;
    my $gender_file = (exists $hashConfig->{"GenderFile"}) ? $hashConfig->{"GenderFile"} : "[NA]";
    my $fastq_dir   = (exists $hashConfig->{"Fastq"})      ? $hashConfig->{"Fastq"}      : "[NA]";
    my $output_dir  = (exists $hashConfig->{"Output"})     ? $hashConfig->{"Output"}     : "[NA]";
    my $report_dir  = (exists $hashConfig->{"Report"})     ? $hashConfig->{"Report"}     : "[NA]";
    my $target_bed  = (exists $hashConfig->{"TargetBed"})  ? $hashConfig->{"TargetBed"}  : "[NA]";
    my $realign_bed = (exists $hashConfig->{"RealignBed"}) ? $hashConfig->{"RealignBed"} : "[NA]";
    my $type        = (exists $hashConfig->{"GenomeBed"})  ? "WGS"                       : "WES";
    my $species     = (exists $hashConfig->{"Species"})    ? $hashConfig->{"Species"}    : "[NA]";
 
    print "case:        $case_num\n";
    print "control:     $control_num\n";
    print "GenderFile:  $gender_file\n";
    print "Fastq:       $fastq_dir\n";
    print "Output:      $output_dir\n";
    print "Report:      $report_dir\n";
    print "TargetBed:   $target_bed\n";
    print "RealignBed:  $realign_bed\n";
    print "Species:     $species\n";
    print "ProjectType: $type\n";
    print "[WARNINGS]:  Distinguish project type by judging the existence of \"GenomeBed\" \n";
}

# 获取物种信息
sub get_species{
    my $hashConfig = shift @_;
    my $species    = (exists $hashConfig->{'Species'}) ? $hashConfig->{'Species'} : 'Lost Species';
    return $species;
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK  = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir  = shift @_;
    my $isOK = 0;
    $isOK    = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 创建目录
sub make_dir{
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
}

#是否继续
sub is_continue{
    print "\n[Option] Confirm and Start?[y/n]";
    my $input = get_input();
    exit if($input ne "y" and $input ne 'Y');
}

# 从屏幕获得输入
sub get_input{
    my $input=<STDIN>;
    $input=~ s/[\r\n]//g;
    return $input;
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
    my @finish_samples = exists $hashCondition->{'Finish'} ? sort keys %{$hashCondition->{'Finish'}} : ();
    print LOG 'Finish   : '. (join ",", @finish_samples) . "\n"; 
    # Good2Run
    my @good2run_samples = exists $hashCondition->{'Good2Run'} ? sort keys %{$hashCondition->{'Good2Run'}} : (); 
    print LOG 'Good2Run : '. (join ",", @good2run_samples) . "\n"; 
    # Error
    my @error_samples = exists $hashCondition->{'Error'} ? sort keys %{$hashCondition->{'Error'}} : ();
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
    my ($workbook) = @_;
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

    $format{'red'} = $workbook->add_format();
    $format{'red'} ->set_align('vcenter');
    $format{'red'} ->set_size(12);
    $format{'red'} ->set_font("Times New Roman");
    $format{'red'} ->set_bg_color("red");
    $format{'red'} ->set_border();  
    
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