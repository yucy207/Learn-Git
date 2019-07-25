package package::mRNA_SNP::library;
use strict;
use warnings;
$|=1;

sub run
{
    my $metadata = shift;
    my $base     = shift;
    my @samples  = split /,/, $metadata->{'samples'};
    my $organ    = qq{$metadata->{organ}};
    my $species  = qq{Human\_$organ};

    my $AnnovarDir   = qq{$base->{AnnovarDir}};
    my $AnnovarBuild = qq{$base->{$organ}{AnnovarBuild}};
    my $genome       = qq{$base->{$organ}{genome_fasta}};
    my $blast_idx    = qq{$base->{$organ}{blast_idx}};
    my $MegaBlast    = qq{$base->{MegaBlast}};
    
    my $result       = qq{$metadata->{project}/mRNA/snp}; 
    my $tmp_dir      = qq{/home/tmp};
    my $library_dir  = qq{$result/result/library};
    system qq{mkdir -p $library_dir} if not -d $library_dir;  

    # 获取注释配置信息列表
    my %hashAnnotationList = package::mRNA_SNP::utils::get_annotation_list($metadata, $base); # 获取当前物种的注释列表
    my %hashKeyWord        = package::mRNA_SNP::annotation_db::get_annotation_key_word($metadata, $base); # 自定义关键词

    ########################################### 开始注释 ###################################################  
    if (-e qq{$library_dir/library_final}){

        print qq{library分析已经完成!\n};
        exit;
    }

    # 1 建立注释库文件
    print "1.Build annovar input\n";
    my $vcf_final        = "$result/result/vcf/sample.all.final.vcf.gz";
    my $library_file     = "$library_dir/library"; # annovar snv 库
    my $vcf2annovar_file = "$library_dir/vcf2annovar";
    my $vcf_first        = "$library_dir/first.alt.vcf";
    keep_vcf_first_alt_allele($vcf_final, $vcf_first);
    system("$AnnovarDir/convert2annovar.pl  --format vcf4old --includeinfo --comment --outfile $vcf2annovar_file $vcf_first");    
    create_annovar_input($vcf2annovar_file, $library_file);  # 创建annovar输入文件
    
    # 2 注释
    print "2.annovar annotation\n";
    # (2-1) annovar 注释
    my $max_threads = 10;
    my $pm = Parallel::ForkManager->new($max_threads);
    foreach my $anno_model(keys %hashAnnotationList)
    {
        next if($hashAnnotationList{$anno_model}{'From'} ne 'Annovar');
        $pm->start() and next;
        my $AnnovarDB = $hashAnnotationList{$anno_model}{'DBDir'}{$species};
        system("$AnnovarDir/$hashAnnotationList{$anno_model}{'Para'} --buildver $AnnovarBuild $library_file $AnnovarDB");
        $pm->finish; 
    }
    $pm->wait_all_children;

    #ACMG 注释
    Intervar_ACMG($metadata, $base, \%hashAnnotationList, \%hashKeyWord) if(exists $hashAnnotationList{'Intervar'});

    #snpEFF 注释
    snpEFF_anno($metadata, $base, $library_file) if(exists $hashAnnotationList{'snpEFF'});

    # #序列/同源性注释
    print "3.homology check\n";
    annotate_homology($library_file, $genome, $blast_idx, $MegaBlast);

    # #STR注释
    # # 功能：(1)indel是否是STR(2)snp上下50bp范围内是否有STR；或者snp上下3bp范围内是否有INDEL
    print "4.STR check\n";
    annotate_STR($library_file, $genome, $AnnovarDir);
  
    # # 3 library_final 基础注释结果生成
    create_library_final($library_file, \%hashAnnotationList, \%hashKeyWord);

    print qq{library分析完成!\n};

}

#snpEFF 注释
sub snpEFF_anno
{
    my $metadata      = shift @_;
    my $base          = shift @_;
    my $library_file  = shift @_;

    my $organ         = qq{$metadata->{organ}};
    my $species       = qq{Human\_$organ};
    my $Java          = qq{$base->{Java}};
    my $snpEFF        = qq{$base->{snpEFF}};
    my $snpEFF_config = qq{$base->{$organ}{snpEFF_config}};
    my $snpEFF_build  = qq{$base->{$organ}{snpEFF_build}};

    #####
    # (1)转vcf格式文件
    #####
    my $library_vcf = "$library_file.vcf";
    open LIBRARY, $library_file;
    open LIBRARYVCF, ">$library_vcf";
    while(<LIBRARY>)
    {
        $_=~s/[\r\n]//g;
        my ($chr, $start, $end, $ref, $alt, $quality, $soft, $vcf_title_list) = split /\t/, $_;
        my $vcf_title = (split /;/, $vcf_title_list)[0];
        my ($chr_vcf, $pos_vcf, $ref_vcf, $alt_list_vcf) = split /\|/, $vcf_title;
        my $alt_vcf = (split /,/, $alt_list_vcf)[0];
        print LIBRARYVCF "$chr_vcf\t$pos_vcf\t.\t$ref_vcf\t$alt_vcf\t.\t$chr|$start|$end|$ref|$alt\t.\n";
    }
    close LIBRARYVCF;
    close LIBRARY;

    #####
    #（2）snpEFF注释
    #####
    my $snpEFF_anno = "$library_vcf.snpEFF.anno";
    my $snpEFF_final = "$snpEFF_anno.final";
    system("$Java -jar $snpEFF -noStats -c $snpEFF_config $snpEFF_build $library_vcf > $snpEFF_anno");
    open ANNO, $snpEFF_anno;
    open FINAL, ">$snpEFF_final";
    while(<ANNO>)
    {
        $_=~s/[\r\n]//g;
        next if($_=~/^#/);
        my ($chr_vcf, $pos_vcf, $id, $ref_vcf, $alt_vcf, $quality, $annovar_title, $snpEFF_anno_info) = split /\t/, $_;
        my ($chr, $start, $end, $ref, $alt) = split /\|/, $annovar_title;
        print FINAL "$chr\t$start\t$end\t$ref\t$alt\t$snpEFF_anno_info\n";
    }
    close ANNO;
    close FINAL;
}

# InterVar ACMG 生成
sub Intervar_ACMG
{
    my $metadata           = shift @_;
    my $base               = shift @_;
    my $hashAnnotationList = shift @_;
    my $hashKeyWord        = shift @_;

    my $organ        = qq{$metadata->{organ}};
    my $species      = qq{Human\_$organ};
    my $report_dir   = qq{$metadata->{project}/mRNA/snp/result}; 
    my $library_dir  = qq{$report_dir/library};
    my $AnnovarDir   = qq{$base->{AnnovarDir}};
    my $Python27     = qq{$base->{Python27}};
    my $InterVar     = qq{$base->{InterVar}};
    my $AnnovarBuild = qq{$base->{$organ}{AnnovarBuild}};

    # 已知疾病注释PP4
    my @find_diseases = split/;/, $metadata->{"Disease"} if(exists $metadata->{"Disease"} and $metadata->{"Disease"}=~/\w+/); # 已知要查找疾病
    my $ACMG_evidence = 'None';
    if(exists $hashAnnotationList->{'HGMD'} and @find_diseases > 0)
    {
        $ACMG_evidence = "$library_dir/ACMG_evidence.txt";# 证据文件
        my %hashHGMD;
        package::mRNA_SNP::utils::read_general_snp_annotation(\%hashHGMD, $hashAnnotationList->{'HGMD'}, $hashKeyWord, 'HGMD') if(exists $hashAnnotationList->{'HGMD'});
        open EVIDENCE, ">$ACMG_evidence";
        foreach my $snp_title(keys %hashHGMD)
        {
            my ($chr, $start, $end, $ref, $alt) = split /\|/, $snp_title;
            my $hgmd_disease = $hashHGMD{$snp_title}{'HGMD_site'}; # HGMD注释疾病

            my $find = 0; # 当前注释里是否包含关注的基因
            foreach my $find_disease (@find_diseases)
            {
                next if($hgmd_disease!~/$find_disease/);
                $find = 1;
            }
            print EVIDENCE "$chr\t$start\t$ref\t$alt\tPP4=1\n" if($find == 1);
        }
        close EVIDENCE;
    }

    # InterVar 运行
    my $intervar_result = "$library_dir/intervar.$AnnovarBuild\_multianno.txt.intervar";
    my $intervar_acmg   = "$library_dir/library.intervar.ACMG";
    my $intervar_db     = "/home/genesky/database/intervar/$AnnovarBuild/2.1.2";
    system("$Python27 $InterVar/Intervar.py -b $AnnovarBuild -i $library_dir/library --input_type=AVinput -o $library_dir/intervar -t $InterVar/intervardb -s $ACMG_evidence --table_annovar=$AnnovarDir/table_annovar.pl --convert2annovar=$AnnovarDir/convert2annovar.pl --annotate_variation=$AnnovarDir/annotate_variation.pl -d $intervar_db");
    # 结果提取，方便后续处理
    open INTERVAR, $intervar_result;
    open ACMG, ">$intervar_acmg";
    while(<INTERVAR>)
    {
        $_=~s/[\r\n]//g;
        next if($_=~/^#/);
        my ($chr, $start, $end, $ref, $alt, $gene, $region, $func, $ensgene, $avsnp, $aachange_ens, $aachange_ref, $clinvar, $ACMG_result, $tmp) = split /\t/, $_, 15;
        my ($class, $evidence) = $ACMG_result =~/InterVar: (.*) (PVS1=.*)/;
        print ACMG "$chr\t$start\t$end\t$ref\t$alt\t$class\t$evidence\n";
    }
    close ACMG;
    close INTERVAR;
}

sub create_library_final
{
    my $library_file       = shift @_;
    my $hashAnnotationList = shift @_;
    my $hashKeyWord        = shift @_;
    my @need_models = ('refGene', 'SNPID', '1000g', 'Homology', 'STR');
    my %hashAnnotation;
    foreach my $need_model(@need_models)
    {
        package::mRNA_SNP::utils::read_general_snp_annotation(\%hashAnnotation, $hashAnnotationList->{$need_model}, $hashKeyWord, $need_model); # 提取常规SNP注释
    }

    my $library_final = "$library_file\_final";
    print "Generate $library_final ... ";
    open LIBRARY, $library_file;
    open FINAL, ">$library_final";
    my $count = 0;
    while(<LIBRARY>)
    {
        $_=~s/[\r\n]//g;
        $count++;
        my ($chr, $start, $end, $ref, $alt, $qual, $tmp) = split /\t/, $_, 7;
        my $snp_title    = "$chr|$start|$end|$ref|$alt";
        my $position     = ($ref eq '-' or $alt eq '-') ? "$start-$end" : $start;
        my $g1000        = (exists $hashAnnotation{$snp_title}{'Freq_Alt (1000g)'})           ? $hashAnnotation{$snp_title}{'Freq_Alt (1000g)'}           : '';
        my $snpid        = (exists $hashAnnotation{$snp_title}{'SNP ID'})                     ? $hashAnnotation{$snp_title}{'SNP ID'}                     : '';
        my $gene         = (exists $hashAnnotation{$snp_title}{'Gene'})                       ? $hashAnnotation{$snp_title}{'Gene'}                       : '';
        my $gene_region  = (exists $hashAnnotation{$snp_title}{'Gene Region'})                ? $hashAnnotation{$snp_title}{'Gene Region'}                : '';
        my $function     = (exists $hashAnnotation{$snp_title}{'Function'})                   ? $hashAnnotation{$snp_title}{'Function'}                   : '';
        my $AAChange     = (exists $hashAnnotation{$snp_title}{'Predicted Protein Variants'}) ? $hashAnnotation{$snp_title}{'Predicted Protein Variants'} : '';
        my $gene_detail  = (exists $hashAnnotation{$snp_title}{'GeneDetail'})                 ? $hashAnnotation{$snp_title}{'GeneDetail'}                 : '';
        my $str_check    = (exists $hashAnnotation{$snp_title}{'STR Check'})                  ? $hashAnnotation{$snp_title}{'STR Check'}                  : '';
        my $seq5         = (exists $hashAnnotation{$snp_title}{'5\'Flanking Sequence'})       ? $hashAnnotation{$snp_title}{'5\'Flanking Sequence'}       : '';
        my $seq3         = (exists $hashAnnotation{$snp_title}{'3\'Flanking Sequence'})       ? $hashAnnotation{$snp_title}{'3\'Flanking Sequence'}       : '';
        my $homhit       = (exists $hashAnnotation{$snp_title}{'Homology Hits'})              ? $hashAnnotation{$snp_title}{'Homology Hits'}              : 0;
        $qual-- if($str_check=~/\w/); # 质量较低，减-
        
        my %hashAAChange_gene = get_AAChange_gene($AAChange); # AAChange基因列表
        my $gene_first        = (split /;/, $gene)[0]; # 第一个基因        
        $gene                 = ($AAChange eq '' or $AAChange =~ /unknown/i or exists $hashAAChange_gene{$gene_first}) ? $gene_first : (sort {$hashAAChange_gene{$a} <=> $hashAAChange_gene{$b}} keys %hashAAChange_gene)[0];
        # 基因的这种选取方式的原因在于，
        # 当两个基因CDS重叠时，同一个突变会有多个基因注释，此时，需要根据AAChange选择基因，而不能直接拿gene_first。因为gene_first可能对应的注释是unknown/synonymous SNV
        
        $gene_detail = '' if($gene_detail =~ /dist/);   # 距离类型的注释不要
        $AAChange    = ($AAChange =~ /\w/) ? "$AAChange;$gene_detail" : $gene_detail;
        $AAChange    =~ s/;$//;
        $snpid      .= $str_check;
        print FINAL "SNV$count\t$qual\t$snpid\t$g1000\t$ref\t$alt\t$chr\t$position\t$gene\t$gene_region\t$function\t$AAChange\t$seq5\t$seq3\t$homhit\n";
    }
    close FINAL;
    close LIBRARY;
    print "OK\n";
}

# 提取AAChange中的基因
sub get_AAChange_gene
{
    my $AAChange = shift @_;
    my %hashAAChange_gene;
    my $count = 0;
    foreach my $AAChange_info(split /,/, $AAChange)
    {   
        $count++;
        my ($gene) = split /:/, $AAChange_info;
        $hashAAChange_gene{$gene} = $count;
    }    
    return %hashAAChange_gene;
}

# STR检查
sub annotate_STR
{
    my $library_file = shift @_;
    my $genome       = shift @_;
    my $AnnovarDir   = shift @_;
    ######
    # 1 提取突变位点前后50bp序列
    ######    
    my $genome_seq50 = "$library_file.seq50bp";
    get_genome_seq50($library_file, $genome_seq50, $genome); # 提取annovarInput格式突变，前后50bp序列

    ######
    # 2 STR注释
    ######
    # 2-1 拆分成多份，每个文件5万个位点，加速
    
    my $library_dir   = package::mRNA_SNP::utils::get_dirname($library_file);
    my $str_split_dir = "$library_dir/str_tmp_dir"; # 拆分临时目录
    package::mRNA_SNP::utils::make_dir($str_split_dir);
    my $str_anno_file = "$str_split_dir/str.summary";

    my %hashSplit = split_file($genome_seq50, $str_split_dir, 50000); # 拆分，每份5万
    my $pm = Parallel::ForkManager->new(10); #固定10线程并行
    foreach my $split_count(sort {$a<=>$b}keys %hashSplit)
    {
        $pm->start() and next;
        my $file = $hashSplit{$split_count};
        my $str_out = "$file.str";
        check_str($file, $str_out);
        $pm->finish; 
    }
    $pm->wait_all_children;
    my @str_results = map{"$hashSplit{$_}.str"}keys %hashSplit;
    system("cat @str_results > $str_anno_file"); # 汇总

    ######
    # 3 SNP附近INDEL检测(INDEL上下3bp有SNP)
    ######
    my $library_dir2 = package::mRNA_SNP::utils::get_dirname($library_file);
    my $indel_db     = "$library_dir2/myIndel_around3bp.txt";
    my $snp_lib      = "$library_dir2/mySNP_library";
    my $snp_around3Indel = "$snp_lib.myIndel_around3bp"; # 注释结果

    open LIBRARY, $library_file;
    open INDEL, ">$indel_db"; # INDEL作为注释数据库
    open SNP, ">$snp_lib"; # SNP作为查找对象
    my $indel_count = 0;
    while(<LIBRARY>)
    {
        $_=~s/[\r\n]//g;
        my ($chr, $start, $end, $ref, $alt, $tmp) = split /\t/, $_;
        if($ref eq '-' or $alt eq '-')
        {
            $indel_count++;
            my $start_down3bp = $start-3;
            my $end_up3bp = $end+3;
            print INDEL "$indel_count\t$chr\t$start_down3bp\t$end_up3bp\t$chr|$start|$end|$ref|$alt\n";
        }
        else
        {
            print SNP "$_\n";
        }
    }
    close LIBRARY;
    close INDEL;
    close SNP;
    system("$AnnovarDir/annotate_variation.pl -regionanno -dbtype around3bp --buildver myIndel $snp_lib $library_dir2");
 
    ######
    # 4 结果汇总
    ######
    my %hashSTR;
    my $str_result_file = "$library_file.str.anno"; # 当前模块最终结果
    # (1) 读取STR注释
    open STR, $str_anno_file;
    while(<STR>)
    {
        $_=~s/[\r\n]//g;
        my ($title, $str_mark) = split /\t/, $_;
        $hashSTR{$title} .= "$str_mark" if(defined $str_mark);
    }
    close STR;
    # (2) 读取INDEL附近SNP 注释
    open INDELAROUND, $snp_around3Indel;
    while(<INDELAROUND>)
    {
        $_=~s/[\r\n]//g;
        my ($name, $anno, $chr, $start, $end, $ref, $alt, $tmp) = split /\t/, $_, 8;
        $hashSTR{"$chr|$start|$end|$ref|$alt"} .= "(INDEL_around)";
    }
    close INDELAROUND;
    # (3) 输出
    open LIBRARY, $library_file;
    open STRRESULT, ">$str_result_file";
    while(<LIBRARY>)
    {
        $_=~s/[\r\n]//g;
        my ($chr, $start, $end, $ref, $alt, $tmp) = split /\t/, $_;
        my $title = "$chr|$start|$end|$ref|$alt";
        print STRRESULT "$chr\t$start\t$end\t$ref\t$alt\t$hashSTR{$title}\n" if(exists($hashSTR{$title}));;
    }
    close LIBRARY;
    close STRRESULT;
}

# STR检测
sub check_str
{
    my $file = shift @_;
    my $str_out = shift @_;
    open IN, $file;
    open OUT, ">$str_out";
    while(<IN>)
    {
        $_=~s/[\r\n]//g;
        my ($title, $seq5, $seq3) = split /\t/, $_;
        my ($chr, $start, $end, $ref, $alt) = split /\|/, $title;
        my $is_snv_str = '';

        if($ref eq '-' or $alt eq '-')
        {
            my $str_find_indel = str_annotation_indel("$seq5$seq3", 6); # 统计2/3/4/5/6重复出现>=6次，且覆盖范围起点<=52,终点>=48的结果数量与对应的motif构成的序列
            my ($indel_str_find, $str_seq) = split /\|/, $str_find_indel;
            if($seq5=~/(([ATCG]{1})\2{7,})$/ig){ # 5seq末端是单碱基重复>=8次，或者3seq首端单碱基重复>=8次
                $str_seq .= ";$1"; 
                $indel_str_find++; 
            }
            if($seq3=~/^(([ATCG]{1})\2{7,})/ig){
                $str_seq .= ";$1"; 
                $indel_str_find++; 
            }
            $is_snv_str = "(STR)" if($indel_str_find>0 &&  ($str_seq=~/$ref/i || $str_seq=~/$alt/i  )); # 如果统计到str，且重复序列中包含了插入、缺失碱基              
        }
        else
        {
            my $str_find_snp = str_annotation_snp("$seq5$seq3", 8); # 100bp序列中1/2个碱基至少重复连续出现8次的数量
            $is_snv_str = "(STR_around)" if($str_find_snp > 0);
        }
        
        $is_snv_str .= "(repeats_or_low_complexity_regions)" if($seq5 =~ /[atcg]$/ and $seq3 =~ /^[atcg]/); # 突变处于高重复或低复杂度区域，小写字符检测

        print OUT "$title\t$is_snv_str\n" if($is_snv_str ne '');
    }
    close IN;
    close OUT;
}
# INDEL STR检测
sub str_annotation_indel
{
    my ($seq, $num) = @_;
    my $str_find = 0;
    my @specs = ([2,$num],
        [3,$num],
        [4,$num],
        [5,$num],
        [6,$num],
    ); 
    my $i = 0;
    my $str_seq_all = "";
    for($i=0; $i<scalar(@specs); $i++){
        my $motiflength = $specs[$i]->[0];
        my $minreps = $specs[$i]->[1] - 1;
        while($seq =~/(([gatc]{$motiflength})\2{$minreps,})/ig){
            my $str_seq = $1;
            my $positon_end = pos($seq);
            my $positon_start = $positon_end - length($str_seq);
            if($positon_start<=52 && $positon_end>=48 ){
                $str_find++;
                $str_seq_all .= "$str_seq;";
            }
        }
    }
    my $result = "$str_find|$str_seq_all";
    return  $result;       
}

sub str_annotation_snp
{
    my ($seq, $num) = @_;
    my $str_find = 0;

    my @specs = ([1,$num],
        [2,$num],
    );
    my $i = 0;
    for($i=0; $i<scalar(@specs); $i++){
        my $motiflength = $specs[$i]->[0];
        my $minreps = $specs[$i]->[1] - 1;
        my $regexp = "(([gatc]{$motiflength})\\2{$minreps,})";
        $str_find++ if($seq =~/(([gatc]{$motiflength})\2{$minreps,})/i);
    }  
    return  $str_find;       
}

# 提取突变前后50bp序列并计算同源性
sub annotate_homology
{
    my $library_file = shift @_;
    my $genome       = shift @_;
    my $blast_idx    = shift @_;
    my $MegaBlast    = shift @_;
    
    ######
    # 1 提取突变位点前后50bp序列
    ######    
    my $genome_seq50 = "$library_file.seq50bp";
    get_genome_seq50($library_file, $genome_seq50, $genome); # 提取annovarInput格式突变，前后50bp序列
    
    ######
    # 2 同源性注释
    ######
    # 2-1 生成fasta序列，并拆分成多份，每个文件5万个位点，加速
    my $genome_seq50_fasta = "$genome_seq50.fa"; # fasta文件
    my $library_dir        = package::mRNA_SNP::utils::get_dirname($genome_seq50_fasta);
    my $fasta_split_dir    = "$library_dir/homhit_tmp_dir"; # 拆分临时目录
    package::mRNA_SNP::utils::make_dir($fasta_split_dir);
    my $homhit_result = "$fasta_split_dir/homhit.summary"; # 临时汇总

    system("more $genome_seq50 |awk '{print \">\"\$1\"\\n\"\$2\$3}' > $genome_seq50.fa"); # 生成fasta
    my %hashSplitFasta = split_fasta($genome_seq50_fasta, $fasta_split_dir, 50000); # 拆分，每份5万
    # 2-2 比对获取同源性
    my $pm = Parallel::ForkManager->new(10); #固定10线程并行
    foreach my $split_count(sort {$a<=>$b}keys %hashSplitFasta)
    {
        $pm->start() and next;
        my $fasta = $hashSplitFasta{$split_count};
        my $blast_out = "$fasta.blastout";
        my $homhit_out = "$fasta.homhit";# 每个文件的同源结果

        system("$MegaBlast -a 2 -d $blast_idx -i $fasta -o $blast_out -F F -D3 -p 80 -s 50 ");# 5万个序列，大约运行15分钟
        get_homology_hit($blast_out, $homhit_out);
        $pm->finish; 
    }
    $pm->wait_all_children; 
    my @homhit_results = map{"$hashSplitFasta{$_}.homhit"} keys %hashSplitFasta;
    system("cat @homhit_results > $homhit_result"); # 汇总

    ######
    # 3 所有结果合并
    ######
    # 3-1 提取同源信息
    my %hashHomhit;
    open HOMHIT, $homhit_result;
    while(<HOMHIT>)
    {
        $_=~s/[\r\n]//g;
        my ($title, $homhit) = split /\t/, $_;
        $hashHomhit{$title} = $homhit;
    }
    close HOMHIT;
    # 3-2 生成最终汇总结果
    my $homology_anno = "$library_file.homology.anno";
    open HOMOLOGYANNO, ">$homology_anno";
    open GENOMESEQ, $genome_seq50;
    while(<GENOMESEQ>)
    {
        $_=~s/[\r\n]//g;
        my ($title, $seq5, $seq3) = split /\t/, $_;
        my ($chr, $start, $end, $ref, $alt) = split /\|/, $title;
        my $homhit = (exists($hashHomhit{$title})) ? $hashHomhit{$title} : 0;
        print HOMOLOGYANNO "$chr\t$start\t$end\t$ref\t$alt\t$seq5\t$seq3\t$homhit\n";
    }
    close GENOMESEQ;
    close HOMOLOGYANNO;
}

# fasta 文件拆分
sub split_fasta
{
    my $fasta_in   = shift @_; # 输入fasta文件
    my $output_dir = shift @_; # 输出路径
    my $limit      = shift @_; # 每个文件最大序列数
    my $fasta_name = package::mRNA_SNP::utils::get_basename($fasta_in);

    open FASTA, $fasta_in;     
    my %hashSplitFasta; # 记录每一个拆分后的fasta文件
    my %hashFastaHandle;# 临时句柄
    my $seq_count = 0;
    my $split_count = 0;
    while(my $name = <FASTA>)
    {
        my $seq = <FASTA>;
        $seq_count++;
        if($seq_count == 1) # 当前组第一条序列
        {   
            my $fasta_split_file = "$output_dir/$fasta_name.$split_count";# 拆分文件名称
            open $hashFastaHandle{'FASTA'}, ">$fasta_split_file";
            $hashSplitFasta{$split_count} = $fasta_split_file;
        }
        my $fasta_handle = $hashFastaHandle{'FASTA'};
        print $fasta_handle "$name$seq";

        if($seq_count == $limit) # 到达上限
        {
            $seq_count = 0;
            $split_count++;
            close $hashFastaHandle{'FASTA'};
            delete $hashFastaHandle{'FASTA'};
        }
    }
    close FASTA;
    close $hashFastaHandle{'FASTA'} if(exists($hashFastaHandle{'FASTA'}));
    return %hashSplitFasta;
}

# 文件拆分
sub split_file
{
    my $file_in    = shift @_; # 输入文件
    my $output_dir = shift @_; # 输出路径
    my $limit      = shift @_; # 每个文件最大序列数
    my $file_name  = package::mRNA_SNP::utils::get_basename($file_in);

    open FILE, $file_in;     
    my %hashSplit; # 记录每一个拆分后的fasta文件
    my %hashFileHandle;# 临时句柄
    my $seq_count = 0;
    my $split_count = 0;
    while(my $info = <FILE>)
    {
        $seq_count++;
        if($seq_count == 1) # 当前组第一条序列
        {   
            my $split_file = "$output_dir/$file_name.$split_count";# 拆分文件名称
            open $hashFileHandle{'FILE'}, ">$split_file";
            $hashSplit{$split_count} = $split_file;
        }
        my $fasta_handle = $hashFileHandle{'FILE'};
        print $fasta_handle "$info";

        if($seq_count == $limit) # 到达上限
        {
            $seq_count = 0;
            $split_count++;
            close $hashFileHandle{'FILE'};
            delete $hashFileHandle{'FILE'};
        }
    }
    close FILE;
    close $hashFileHandle{'FILE'} if(exists($hashFileHandle{'FILE'}));
    return %hashSplit;
}
# 提取annovarInput格式突变，前后50bp序列
sub get_genome_seq50
{
    my $aviInput    = shift @_;
    my $output_file = shift @_;
    my $genome      = shift @_;
    my %hashGenome  = package::mRNA_SNP::utils::read_fasta($genome);
    open AVINPUT, $aviInput;
    open GENOMESEQ, ">$output_file";
    while(<AVINPUT>)
    {
        $_=~s/[\r\n]//g;
        my ($chr, $start, $end, $ref, $alt, $tmp) = split /\t/, $_;
        my $title = "$chr|$start|$end|$ref|$alt";
 
        my $lenth_forword = ($ref eq "-") ? 50 :51; # 前50bp提取时，需要向前移动的长度，注substr第一个bp坐标为0
        my $seq5=substr($hashGenome{$chr}, $start-$lenth_forword, 50);
        my $seq3=substr($hashGenome{$chr}, $end, 50);
        print GENOMESEQ "$title\t$seq5\t$seq3\n";
    }
    close AVINPUT;
    close GENOMESEQ;    
}

# 从比对结果获取同源信息
sub get_homology_hit
{
    my $blast_out = shift @_;
    my $homhit_out = shift @_;

    # 计算同源数量
    my %hashHomhit;
    open BLAST,$blast_out;
    while(<BLAST>){
        next if(/^#/);
        my ($title, $chr, $identity, $align_length, $temp)=split(/\t/,$_,5);
        $hashHomhit{$title}++ if($align_length > 90);# 比对长度大于90，则增加一个同源记录
    }
    close BLAST;
    # 输出到文件
    open HOMHIT, ">$homhit_out";
    foreach my $title(sort keys %hashHomhit)
    {
        print HOMHIT "$title\t$hashHomhit{$title}\n";
    }
    close HOMHIT;  
}

# 取vcf ALT列的第一个等位基因注释
sub keep_vcf_first_alt_allele
{
    my $vcf_file  = shift @_;
    my $vcf_first = shift @_;
    my $is_gzip   = ($vcf_file=~/\.gz$/) ? 1 : 0;# VCF是否是压缩的

    open VCF, "gzip -cd $vcf_file|" if($is_gzip == 1);
    open VCF, $vcf_file             if($is_gzip == 0);
    open VCFFIRST, ">$vcf_first";
    print VCFFIRST "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    while(<VCF>)
    {
        next if($_=~/^#/);
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter,$info,$tmp) = split /\t/, $_, 9;
        my $alt_first = "";
        map{$alt_first = $_ if($_ =~ /[ATCG]/i and $alt_first eq ""); } (split /,/, $alt); # 保留第一个非“*”的突变
        next if($alt_first eq "");

        print VCFFIRST "$chr\t$pos\t$id\t$ref\t$alt_first\t$qual\t$filter\t$chr|$pos|$ref|$alt\n";
    }
    close VCF;
    close VCFFIRST;
}

# 创建annovar输入文件
sub create_annovar_input
{
    my $vcf2annovar_file = shift @_;
    my $library_file     = shift @_;

    # 提取格式转化结果
    my %hashLibrary;
    open VCF2ANNOVAR, $vcf2annovar_file;    
    my $row = 0;
    while(<VCF2ANNOVAR>)
    {
        $_=~s/[\r\n]//g;
        next if($_=~/^#/);
        $row++;
        my ($chr_annovar, $start_annovar, $end_annovar, $ref_annovar, $alt_annovar, $chr_vcf, $pos_vcf, $id_vcf, $ref_vcf, $alt_vcf, $qual, $filter, $vcf_title_original, $tmp) = split /\t/, $_, 14;
        $end_annovar = $start_annovar + length($ref_annovar) - 1; # 尾部位置注释有时会有问题，所以重新计算
        $alt_annovar = '-' if($alt_annovar!~/\w/ or $alt_annovar=~/\d/); # 存在*数值的情况，此时默认是缺失
        my $quality_minus = ($filter eq 'PASS') ? 0 : 1;

        my $annovar_title = "$chr_annovar\t$start_annovar\t$end_annovar\t$ref_annovar\t$alt_annovar";
        my $vcf_title = "$chr_vcf|$pos_vcf|$ref_vcf|$alt_vcf";
        $hashLibrary{$annovar_title}{'VCF_TITLE'} .= "$vcf_title_original;"; # 会存在重复位点
        $hashLibrary{$annovar_title}{'SortOrder'} = $row;
        $hashLibrary{$annovar_title}{'Quality'}   = 3 - $quality_minus;
        $hashLibrary{$annovar_title}{'Filter'}    = $filter; # VCF过滤标志
    }
    close VCF2ANNOVAR;

    # 建库
    open LIBRARY, ">$library_file";
    open VCFFILTER, ">$library_file.vcf.filter";
    foreach my $annovar_title(sort {$hashLibrary{$a}{'SortOrder'} <=> $hashLibrary{$b}{'SortOrder'}} keys %hashLibrary)
    {
        $hashLibrary{$annovar_title}{'VCF_TITLE'}=~s/;$//;
        print LIBRARY "$annovar_title\t$hashLibrary{$annovar_title}{'Quality'}\tGATK\t$hashLibrary{$annovar_title}{'VCF_TITLE'}\n";
        print VCFFILTER "$annovar_title\t$hashLibrary{$annovar_title}{'Filter'}\n";
    }
    close LIBRARY;
    close VCFFILTER;
}


1