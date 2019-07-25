package package::mRNA_SNP::utils;
use strict;
use warnings;
$|=1;

# 创建目录
sub make_dir
{
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
}
# 获取路径文件名
sub get_basename
{
    my $path = shift @_;
    my $path_curf = File::Spec->rel2abs($path);
    my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
    return $file;
}

# 获取文件路径名
sub get_dirname
{
    my $path = shift @_;
    my $path_curf = File::Spec->rel2abs($path);
    my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
    return $dirs;
}

# 获取注释列表
sub get_annotation_list
{
    my $metadata    = shift @_;
    my $base        = shift @_;
    my $organ       = qq{$metadata->{organ}};
    my $species     = qq{Human\_$organ};
    my $library_dir = qq{$metadata->{project}/mRNA/snp/result/library};
    my $anno_log    = qq{$library_dir/anno.log};

    my %hashKeyWord = package::mRNA_SNP::annotation_db::get_annotation_key_word($metadata, $base); # 获取关键字哈希

    open ANNOLOG, ">$anno_log";
    my %hashAnnotation_all = package::mRNA_SNP::annotation_db::get_annotation(); # 所有注释数据库
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

# 检测必要文件
sub check_annotation_necessary_file
{
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
        my $condition = (package::mRNA_SNP::utils::is_file_ok($need_file)==1) ? 'OK' : 'Lost';
        next if($condition eq 'OK');
        push @errors, "Lost $need_file";
    }
    
    return (join ",", @errors);
}

# 拆分注释信息
sub split_annotation_info
{
    my $info = shift @_;
    my %hashInfo;
    foreach(split /\t/, $info)
    {
        my ($name, $value) = split /=/, $_;
        $hashInfo{$name} = $value;
    }
    return %hashInfo;
}

# 检验文件是否为空
sub is_file_ok
{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# 提取fasta序列
sub read_fasta
{
    my $fasta = shift @_;
    my %hashFasta;
    print "read $fasta ... ";
    my $IN = Bio::SeqIO->new(-file => $fasta, -format=>'Fasta') or die "Could not open up file $fasta: $!";
    while(my $inSeq = $IN->next_seq)
    {
        my $id  = $inSeq->display_id;
        my $seq = $inSeq->seq;
        $hashFasta{$id}= $seq;
    }
    print "\n";
    return %hashFasta;
}

# 读取通用SNP注释
sub read_general_snp_annotation
{
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
                my $chr   = $datas[$hashNeed_Info{$info_name}{'chr_index'}];
                my $start = $datas[$hashNeed_Info{$info_name}{'start_index'}];
                my $end   = $datas[$hashNeed_Info{$info_name}{'end_index'}];
                my $ref   = $datas[$hashNeed_Info{$info_name}{'ref_index'}];
                my $alt   = $datas[$hashNeed_Info{$info_name}{'alt_index'}];
                my $snp_title   = "$chr|$start|$end|$ref|$alt"; # SNP表头
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
                
                $value =~ s/\\x2c/,/g; # annovar 注释内容会把‘，’或‘#’替换为特殊字符，防止冲突，所以使用时需要再替换回来
                $value =~ s/\\x23/#/g;

                $hashAnnotation->{$snp_title}{$info_name} = $value if($value  ne "");
            }
        }
        print "OK\n";
    }
 
}

sub split_general_snp_need_info
{
    my $need_info = shift @_;
    my %hashNeed_Info;
    foreach my $info_detail(split /\t/, $need_info)
    {
        my ($info_name, $value_index, $snp_title_index) = split /\=/, $info_detail;
        my ($chr_index, $start_index, $end_index, $ref_index, $alt_index) = split /,/, $snp_title_index;
        $hashNeed_Info{$info_name}{'value_index'} = $value_index;
        $hashNeed_Info{$info_name}{'chr_index'}   = $chr_index;
        $hashNeed_Info{$info_name}{'start_index'} = $start_index;
        $hashNeed_Info{$info_name}{'end_index'}   = $end_index;
        $hashNeed_Info{$info_name}{'ref_index'}   = $ref_index;
        $hashNeed_Info{$info_name}{'alt_index'}   = $alt_index;
    }
    return %hashNeed_Info;
}

# 检查VCF标题是否缺失样本
sub check_vcf_sample_isOK
{
    my $vcf_file = shift @_;
    my $samples = shift @_;
    my $is_gzip = ($vcf_file=~/\.gz$/) ? 1 : 0;# VCF是否是压缩的

    my @losts;
    open VCF, "gzip -cd $vcf_file|" if($is_gzip == 1);
    open VCF, $vcf_file             if($is_gzip == 0);
    while(<VCF>)
    {
        next if($_=~/^##/);
        if($_=~/^#CHROM/)
        {
            $_=~s/[\r\n]//g;
            my @vcf_heads = split /\t/, $_;
            foreach my $sample(@$samples)
            {
                next if($sample ~~ @vcf_heads);
                push @losts, $sample;
            }
            last;
        }
    }
    close VCF;
    return @losts;
}

# 读取通用GENE注释
sub read_general_gene_annotation
{
    my $hashAnnotation     = shift @_;
    my $hashAnnotationList = shift @_;
    my $hashKeyWord        = shift @_;
    my $need_model         = shift @_;
    my $species            = shift @_;
    foreach my $out_index(keys %{$hashAnnotationList->{'NecessaryFile'}})
    {
        my %hashFile_Info = split_annotation_info($hashAnnotationList->{'NecessaryFile'}{$out_index}); # 文件信息
        my %hashNeed_Info = split_general_gene_need_info($hashAnnotationList->{'Output'}{$out_index}{'Need_Info'}); # 需要的信息坐标

        # 文件检测
        my $anno_dir = (exists $hashAnnotationList->{'DBDir'}{$species}) ? $hashAnnotationList->{'DBDir'}{$species}."/": "";
        my $prefix = (exists($hashKeyWord->{$hashFile_Info{'prefix'}})) ? $hashKeyWord->{$hashFile_Info{'prefix'}} : $hashFile_Info{'prefix'};
        my $suffix = $hashFile_Info{'suffix'};
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

                $hashAnnotation->{$gene}{$info_name} = $value if($value ne "");
            }

        }
        print "OK\n";
    }
 
}

sub split_general_gene_need_info
{
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


1