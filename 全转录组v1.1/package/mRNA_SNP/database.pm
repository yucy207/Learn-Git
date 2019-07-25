package package::mRNA_SNP::database;
use strict;
use warnings;

sub run
{
    
    my $metadata = shift;
    my $base     = shift;
    my @samples  = split /,/, $metadata->{'samples'};

    my $result        = qq{$metadata->{project}/mRNA/snp}; 
    my $database_dir  = qq{$result/result/database};
    system qq{mkdir -p $database_dir} if not -d $database_dir;

    # 输出文件
    my $database_snv  = qq{$database_dir/database.snv};
    my $analysis_all  = qq{$database_dir/analysis_all.txt};
    my $analysis_func = qq{$database_dir/analysis.txt};
    if(package::mRNA_SNP::utils::is_file_ok($database_snv, $analysis_all, $analysis_func) == 1)
    {
        print "[Note] Database Process None, for reProcess Delete Result\n";
        return;
    }

    # 输入文件
    my $gender_log    = qq{$result/result/gender/gender.log};
    my $library       = qq{$result/result/library/library};
    my $library_final = qq{$result/result/library/library_final};
    my $vcf_final     = qq{$result/result/vcf/sample.all.final.vcf.gz};
    die "Lost : $gender_log\n"    if(package::mRNA_SNP::utils::is_file_ok($gender_log) == 0);
    die "Lost : $library\n"       if(package::mRNA_SNP::utils::is_file_ok($library) == 0);
    die "Lost : $library_final\n" if(package::mRNA_SNP::utils::is_file_ok($library_final) == 0);
    die "Lost : $vcf_final\n"     if(package::mRNA_SNP::utils::is_file_ok($vcf_final) == 0);
    my @losts    = package::mRNA_SNP::utils::check_vcf_sample_isOK($vcf_final, \@samples); # VCF文件样本检测，防止丢失
    die "[ERR] Lost Sample in VCF FILE: @losts\n" if(@losts > 0);

    # 开始处理
    my %hashGender       = package::mRNA_SNP::gender::read_gender_file($gender_log);
    my %hashLibraryVCF   = read_library_vcf($library);
    my %hashLibraryFinal = read_library_final($library_final);
    open DATABASESNV, ">$database_snv";
    open ANALYSISALL, ">$analysis_all";
    open ANALYSISFUNC, ">$analysis_func";
    my %hashAnnovarFinished; # VCF转annovar后，存在重复的结果，这里记录已经处理过的位点，防止重复
    my @heads;
    open VCF, "gzip -cd $vcf_final|";
    print "Process $vcf_final ... ";
    while(<VCF>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/ or $_=~/^##/);
        if($_=~/^#CHROM/)
        {
            $_=~s/^#+//;
            $_=~s/[\r\n]//g;
            @heads = split /\t/, $_;
            next;
        }
        # 数据放入临时哈希
        my %hashTmp;
        my @datas = split /\t/, $_;
        foreach my $col(0..$#heads)
        {
            $hashTmp{$heads[$col]} = $datas[$col];
        }
        my %hashVCFData   = format_vcf(\%hashTmp, \@samples, \%hashGender); # 数据格式转换
        my $vcf_title     = $hashVCFData{'VCF_TITLE'};
        my $annovar_title = exists $hashLibraryVCF{$vcf_title} ? $hashLibraryVCF{$vcf_title} : "";
        my $alt_freq      = exists($hashVCFData{'ALT_FREQ'}) ? $hashVCFData{'ALT_FREQ'} : 0;
        next if($alt_freq==0 or $annovar_title eq '');
        $hashAnnovarFinished{$annovar_title}++;
        next if($hashAnnovarFinished{$annovar_title} >=2);
 
        # 输出
        my ($chr, $start, $end, $ref, $alt) = split /\|/, $annovar_title;
        my $pos          = ($ref eq '-' or $alt eq '-') ? "$start-$end" : "$start";
        my $output_title = "$chr|$pos|$ref";
        my $is_func      = $hashLibraryFinal{$annovar_title}{'ISFUNC'}; # 突变是否有功能

        # (1) database.snv 输出
        print DATABASESNV "$hashLibraryFinal{$annovar_title}{'Data'}\n";

        # (2) database.genotyping 输出，并为analysis文件输出准备
        my ($homr_info, $het_info, $homa_info) = ('', '', '');
        my $homr_annovar = "$ref/$ref";
        my $het_annovar  = "$ref/$alt";
        my $homa_annovar = "$alt/$alt";
        my $geno_count   = 0; # 分型样本数
        my $high_quality_count = 0; # 高质量分型数量
        foreach my $sample(@samples)
        {
            my $geno        = exists($hashVCFData{'DATA'}{$sample}{'Geno'})       ? $hashVCFData{'DATA'}{$sample}{'Geno'}       : '';
            my $type        = exists($hashVCFData{'DATA'}{$sample}{'Type'})       ? $hashVCFData{'DATA'}{$sample}{'Type'}       : '';
            my $ref_depth   = exists($hashVCFData{'DATA'}{$sample}{'REF_DEPTH'})  ? $hashVCFData{'DATA'}{$sample}{'REF_DEPTH'}  : '';
            my $alt_depth   = exists($hashVCFData{'DATA'}{$sample}{'ALT_DEPTH'})  ? $hashVCFData{'DATA'}{$sample}{'ALT_DEPTH'}  : '';
            my $GQ          = exists($hashVCFData{'DATA'}{$sample}{'GQ'})         ? $hashVCFData{'DATA'}{$sample}{'GQ'}         : 0;
            my $PL          = exists($hashVCFData{'DATA'}{$sample}{'PL'})         ? $hashVCFData{'DATA'}{$sample}{'PL'}         : 0;
            next if($geno eq '');
            $PL=~s/,/:/g;
            $geno = $homr_annovar if($type eq 'HOMR');# 转换成annovar格式
            $geno = $het_annovar  if($type eq 'HET');
            $geno = $homa_annovar if($type eq 'HOMA');
            
            # analysis 文件格式
            my $analysis_format = "$sample:$type:$geno:$ref_depth:$alt_depth:$GQ:$PL"; 
            $homr_info .= "$analysis_format," if($type eq 'HOMR');
            $het_info  .= "$analysis_format," if($type eq 'HET');
            $homa_info .= "$analysis_format," if($type eq 'HOMA');
            $high_quality_count++ if($GQ>20);
            $geno_count++;
        }
        
        # （3）analysis_all.txt,analysis.txt输出
        $homr_info = ',' if($homr_info eq '');
        $het_info  = ',' if($het_info  eq '');
        $homa_info = ',' if($homa_info eq '');
        my $high_quality_perc = sprintf "%0.4f", $high_quality_count/$geno_count;
        print ANALYSISALL  "$output_title\t$alt\t$high_quality_perc\t$homr_info\t$het_info\t$homa_info\n";
        print ANALYSISFUNC "$output_title\t$alt\t$high_quality_perc\t$homr_info\t$het_info\t$homa_info\n" if($is_func==1);
    }
    close VCF;
    close DATABASESNV;
    close ANALYSISALL;
    close ANALYSISFUNC;
    print "OK\n";
    
    print "database分析完成\n";

}

sub format_vcf
{
    my $hashTmp    = shift @_;
    my $samples    = shift @_;
    my $hashGender = shift @_;
    my %hashVCFData;

    my $chr = $hashTmp->{'CHROM'};
    my $pos = $hashTmp->{'POS'};
    my $ref = $hashTmp->{'REF'};
    my $alt = $hashTmp->{'ALT'};
    my $vcf_title = "$chr|$pos|$ref|$alt";
    $hashVCFData{'VCF_TITLE'} = $vcf_title;
    my ($geno0, $geno1, $geno2) = (0, 0, 0);
    my @formats = split /\:/, $hashTmp->{'FORMAT'};
    foreach my $sample(@$samples)
    {   
        my @datas = split /\:/, $hashTmp->{$sample};
        my %hashSampleResult;
        foreach my $col(0..$#formats)
        {
            $hashSampleResult{$formats[$col]} = $datas[$col];
        }
        my $geno_vcf = $hashSampleResult{'GT'};
        $geno_vcf =~ s/\|/\//g;
        my $depth    = $hashSampleResult{'AD'};
        my $GQ       = $hashSampleResult{'GQ'};
        my $PL       = $hashSampleResult{'PL'};
        my @depths    = split /,/, $depth;
        my $ref_depth = $depths[0];
        my $alt_depth = $depths[1];
        my $depth_sum = $ref_depth + $alt_depth;

        next if($geno_vcf!~/\d/);
        next if($depth_sum<5); # 测序深度限制
        next if($chr eq "Y" and $hashGender->{$sample} eq "female"); # 女性Y染色体跳过

        my $geno = '';
        my $type = '';
        if($geno_vcf eq '0/0')
        {
            $geno = "$ref/$ref";
            $type = "HOMR";
            $geno0++;
        }
        if($geno_vcf eq '0/1' or $geno_vcf eq '1/0')
        {
            $geno = "$ref/$alt";
            $type = "HET";
            $geno1++;
        }
        if($geno_vcf eq '1/1')
        {
            $geno = "$alt/$alt";
            $type = "HOMA";
            $geno2++;
        }

                
        if(($chr eq "X" or $chr eq "Y") and $hashGender->{$sample} eq "male" and $type eq "HET") # 男性的性染色体，如果是杂合突变，改成纯合突变
        {
            $geno = "$alt/$alt";
            $type = "HOMA";
            $geno1--;
            $geno2++;
        }
        $hashVCFData{'DATA'}{$sample}{'Geno'}       = $geno;
        $hashVCFData{'DATA'}{$sample}{'Type'}       = $type;
        $hashVCFData{'DATA'}{$sample}{'REF_DEPTH'}  = $ref_depth;
        $hashVCFData{'DATA'}{$sample}{'ALT_DEPTH'}  = $alt_depth;
        $hashVCFData{'DATA'}{$sample}{'GQ'}         = $GQ;
        $hashVCFData{'DATA'}{$sample}{'PL'}         = $PL;
    }
    my $genosum = $geno0 + $geno1 + $geno2;
    $hashVCFData{'ALT_FREQ'} = ($geno1 + 2*$geno2) / (2*$genosum) if($genosum > 0);
    return %hashVCFData;
}

# 读取library_final
sub read_library_final
{
    my $library_final = shift @_;
    print "Read $library_final ... ";
    my %hashLibraryFinal;
    open LIBRARYFINAL, $library_final;
    while(<LIBRARYFINAL>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        my ($chr, $pos, $ref, $alt)    = ($datas[6], $datas[7], $datas[4], $datas[5]);
        my ($region, $function, $hgvs) = ($datas[9], $datas[10], $datas[11]);
        my $is_func = (($region eq "exonic" or $region eq "splicing" or $region eq "exonic;splicing") and $function ne "synonymous SNV" and $hgvs=~ /\w/ ) ? 1 : 0;# 有功能
        my @positions = split /-/, $pos;
        my $start = $positions[0];
        my $end   = $positions[$#positions];
        my $annovar_title = "$chr|$start|$end|$ref|$alt";
        my $db_title      = "$chr|$pos|$ref"; 
        $datas[0] = $db_title;
        $hashLibraryFinal{$annovar_title}{'Data'}   = join "\t", @datas;
        $hashLibraryFinal{$annovar_title}{'ISFUNC'} = $is_func;
    }
    close LIBRARYFINAL;
    print "OK\n";
    return %hashLibraryFinal;
}

# 读取VCF与annovar格式的对应关系
sub read_library_vcf
{
    my $library = shift @_;
    print "Read $library ... ";
    my %hashLibraryVCF;
    open LIBRARYVCF, $library;
    while(<LIBRARYVCF>)
    {
        $_=~s/[\r\n]//g;
        my ($chr, $start, $end, $ref, $alt, $qual, $gatk, $vcf_titles) = split /\t/, $_;
        my $annovar_title = "$chr|$start|$end|$ref|$alt";
        foreach my $vcf_title(split /;/, $vcf_titles)
        {   
            $hashLibraryVCF{$vcf_title} = $annovar_title;
        }       
    }
    close LIBRARYVCF;
    print "OK\n";
    return %hashLibraryVCF;
}


1