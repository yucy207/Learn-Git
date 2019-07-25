package package::annotation;
use strict;
use warnings;

# 获取注释关键字
sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start Annotation ".package::utils::get_time()." ##########\n";
    my $report_dir = $hashConfig->{'Report'}; 
    my $species    = package::utils::get_species($hashConfig);
    my $annotation = "$report_dir/annotation.txt";
    #####
    # 写日志
    #####
    package::utils::write_log_simple("$report_dir/run.log", "Annotation");
    if(package::utils::is_file_ok($annotation) == 1 and package::utils::is_force_sample($hashConfig) == 0)
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }

    
    my $database_snv       = "$report_dir/database.snv";
    my %hashAnnotationList = package::utils::get_annotation_list($hashPara, $hashConfig); # 获取当前物种的注释列表
    my %hashKeyWord        = package::annotation_db::get_annotation_key_word($hashPara, $hashConfig); # 自定义关键词

    my %hashDatabaseSNV = read_database_snv($database_snv);
    my @snp_titles      = sort {$hashDatabaseSNV{$a}{'SortOrder'} <=> $hashDatabaseSNV{$b}{'SortOrder'}} keys %hashDatabaseSNV;
    my @no_need_models  = ('refGene', 'SNPID', '1000g', 'Homology', 'STR'); # 排除掉的注释内容
 
    
    open ANNOTATION, ">$annotation";
    foreach my $need_model(sort keys %hashAnnotationList)
    {   
        next if($need_model ~~ @no_need_models);
        my %hashAnnotation;
        my $method_get_result = $hashAnnotationList{$need_model}{'GetResult'};
        package::utils::read_general_snp_annotation(\%hashAnnotation, $hashAnnotationList{$need_model}, \%hashKeyWord, $need_model) if($method_get_result eq 'General SNP'); # 提取常规SNP注释
        package::utils::read_general_gene_annotation(\%hashAnnotation, $hashAnnotationList{$need_model}, \%hashKeyWord, $need_model, $species) if($method_get_result eq 'General Gene'); # 提取常规SNP注释
        foreach my $snp_title(@snp_titles)
        { 
            my $chr  = $hashDatabaseSNV{$snp_title}{'Chr'};
            my $pos  = $hashDatabaseSNV{$snp_title}{'Pos'};
            my $ref  = $hashDatabaseSNV{$snp_title}{'Ref'};
            my $alt  = $hashDatabaseSNV{$snp_title}{'Alt'};
            my $gene = $hashDatabaseSNV{$snp_title}{'Gene'};
            next if($method_get_result eq 'General SNP' and !exists($hashAnnotation{$snp_title}));
            next if($method_get_result eq 'General Gene' and !exists($hashAnnotation{$gene}));
            my $key_title = ($method_get_result eq 'General SNP') ? $snp_title : $gene;

            foreach my $anno_name(keys %{$hashAnnotation{$key_title}})
            {   
                print ANNOTATION "$chr|$pos|$ref\t$alt\t$anno_name\t$hashAnnotation{$key_title}{$anno_name}\n" if($hashAnnotation{$key_title}{$anno_name} ne "");
            }
        }
    }  
    close ANNOTATION;

}

sub read_database_snv{
    my $database_snv = shift @_;

    print "Read $database_snv ... ";
    my %hashDatabaseSNV;
    my $count = 0;
    open DATABASESNV, $database_snv;
    while(<DATABASESNV>)
    {
        $_=~s/[\r\n]//g;
        $count++;
        my ($db_title, $qual, $snpid, $g1000, $ref, $alt, $chr, $pos, $gene, $tmp) = split /\t/, $_, 10;
        my @positions = split /-/, $pos;
        my $start = $positions[0];
        my $end   = $positions[$#positions];
        my $snp_title = "$chr|$start|$end|$ref|$alt";
        $hashDatabaseSNV{$snp_title}{'Gene'} = $gene;
        $hashDatabaseSNV{$snp_title}{'Chr'} = $chr;
        $hashDatabaseSNV{$snp_title}{'Pos'} = $pos;
        $hashDatabaseSNV{$snp_title}{'Ref'} = $ref;
        $hashDatabaseSNV{$snp_title}{'Alt'} = $alt;
        $hashDatabaseSNV{$snp_title}{'SortOrder'} = $count;
    }
    close DATABASESNV;
    print "OK\n";
    return %hashDatabaseSNV;
}

1