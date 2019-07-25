package package::mRNA_SNP::annotation;
use strict;
use warnings;

# 获取注释关键字
sub run
{
    my $metadata = shift;
    my $base     = shift;
    my @samples  = split /,/, $metadata->{'samples'}; 
    my $organ    = qq{$metadata->{organ}};
    my $species  = qq{Human\_$organ};

    my $result     = qq{$metadata->{project}/mRNA/snp};
    my $annotation = qq{$result/result/annotation.txt};
 
    if (-e qq{$result/result/annotation.txt}){

        print qq{annotation分析已经完成!\n};
        exit;
    }
    
    my $database_snv       = "$result/result/database/database.snv";
    my %hashAnnotationList = package::mRNA_SNP::utils::get_annotation_list($metadata, $base); # 获取当前物种的注释列表
    my %hashKeyWord        = package::mRNA_SNP::annotation_db::get_annotation_key_word($metadata, $base); # 自定义关键词

    my %hashDatabaseSNV = read_database_snv($database_snv);
    my @snp_titles      = sort {$hashDatabaseSNV{$a}{'SortOrder'} <=> $hashDatabaseSNV{$b}{'SortOrder'}} keys %hashDatabaseSNV;
    my @no_need_models  = ('refGene', 'SNPID', '1000g', 'Homology', 'STR'); # 排除掉的注释内容
 
    
    open ANNOTATION, ">$annotation";
    foreach my $need_model(sort keys %hashAnnotationList)
    {   
        next if($need_model ~~ @no_need_models);
        my %hashAnnotation;
        my $method_get_result = $hashAnnotationList{$need_model}{'GetResult'};
        package::mRNA_SNP::utils::read_general_snp_annotation(\%hashAnnotation, $hashAnnotationList{$need_model}, \%hashKeyWord, $need_model) if($method_get_result eq 'General SNP'); # 提取常规SNP注释
        package::mRNA_SNP::utils::read_general_gene_annotation(\%hashAnnotation, $hashAnnotationList{$need_model}, \%hashKeyWord, $need_model, $species) if($method_get_result eq 'General Gene'); # 提取常规SNP注释
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

    print qq{Annotation分析完成!\n};

}

sub read_database_snv
{
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