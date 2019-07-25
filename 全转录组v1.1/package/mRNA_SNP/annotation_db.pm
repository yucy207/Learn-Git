package package::mRNA_SNP::annotation_db;
use strict;
use warnings;

# 获取注释关键字
sub get_annotation_key_word
{

    my $metadata = shift;
    my $base     = shift;
    my $organ    = qq{$metadata->{organ}};
    my $species  = qq{Human\_$organ};

    my $AnnovarBuild = qq{$base->{$organ}{AnnovarBuild}};
    my $report_dir   = qq{$metadata->{project}/mRNA/snp/result};

    my %hashKeyWord;
    $hashKeyWord{'AnnovarBuild'}                 = qq{$AnnovarBuild};
    $hashKeyWord{'outfile_library'}              = qq{$report_dir/library/library};
    $hashKeyWord{'outfile_library.AnnovarBuild'} = qq{$report_dir/library/library.$AnnovarBuild};
    $hashKeyWord{'GenomeFasta'}                  = qq{$base->{$organ}{genome_fasta}};
    return %hashKeyWord;
}


# 获取注释列表
sub get_annotation{
    my %hashAnnotationList;
    # 关键字说明
    # （1）Species       : 该注释对应的物种，用"[物种名]"表示,如果对应多个物种，用","分割
    # （2）From          : 当前注释的来源方式，Annovar表示需要用annovar进行注释（后面的Para必须提供）;File表示需要写代码从文件提取
    # （3）Para          : Annovar注释时，该参数提供注释参数
    # （4）DBDir         : 当前注释文件所在路径
    # （5）NecessaryFile : 注释时必须存在的文件，由"prefix=...;suffix=..."构成,关键字有"AnnovarBuild",代码会识别该关键字进行相应路径替换(部分注释会需要多个依赖文件，用不同的数字编号key区分)
    # （6）GetResult     : 结果获取方式，"General SNP"表示通用SNP处理代码即可处理，"General Gene"表示同游GENE处理代码即可处理
    # （7）Output:File   : 结果输出文件，由"prefix=...;suffix=..."构成,关键字有"outfile.AnnovarBuild"
    # （8）Output:Need_Info: 该注释需要的信息，两种构成方式：（1）SNP处理模式：名称=信息位置=SNP识别位置（2）基因处理模式：名称=位置信息=Gene识别位置;(1,2)中如果结果是复合结果，可用类似1/,/0提取,表示第一列用‘，’分割，并获取分割后的第0列信息
    # （9）支持多列数据合并在一起，并用字符"|"分隔.例如HGMD   HGMD_site=1/,/1,2=2,3,4,5,6   表示，把第一列用","拆开，然后把1,2两个坐标的数据都保留下来用"|"分隔，并命名为HGMD_site。同时也可以支持类似下面的格式： SNP ID=0,1=2,3,4,5,6。基因注释也同样支持这种操作。  注意，不支持 HGMD_site=0,1/,/1,2=2,3,4,5,6 这种类型
    ###
    # 基础注释
    ###
    $hashAnnotationList{'refGene'}{'Species'}                 = '[all]';
    $hashAnnotationList{'refGene'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'refGene'}{'Para'}                    = 'table_annovar.pl -remove -protocol refGene -operation g --argument \'--splicing_threshold 15 --hgvs\'';
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/ucsc/hg19/gene";
    # $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/refgene";
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/ucsc/hg19/gene";
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/ucsc/hg38/gene";
    $hashAnnotationList{'refGene'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_refGene.txt";
    $hashAnnotationList{'refGene'}{'NecessaryFile'}{2}        = "prefix=AnnovarBuild\tsuffix=_refGeneMrna.fa";
    $hashAnnotationList{'refGene'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'refGene'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_multianno.txt";# 该注释要包含的注释内容
    $hashAnnotationList{'refGene'}{'Output'}{1}{'Need_Info'}  = "Gene Region=5=0,1,2,3,4\tGene=6=0,1,2,3,4\tGeneDetail=7=0,1,2,3,4\tFunction=8=0,1,2,3,4\tPredicted Protein Variants=9=0,1,2,3,4";# 该注释要包含的注释内容

    $hashAnnotationList{'SNPID'}{'Species'}                   = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'SNPID'}{'From'}                      = 'Annovar';
    $hashAnnotationList{'SNPID'}{'Para'}                      = 'annotate_variation.pl -filter -dbtype avsnp150';
    $hashAnnotationList{'SNPID'}{'DBDir'}{'Human_hg19'}       = "/home/genesky/database/annovar/hg19/avsnp";
    $hashAnnotationList{'SNPID'}{'DBDir'}{'Human_hg19_HLA'}   = "/home/genesky/database/annovar/hg19/avsnp";
    $hashAnnotationList{'SNPID'}{'DBDir'}{'Human_hg38'}       = "/home/genesky/database/annovar/hg38/avsnp";
    $hashAnnotationList{'SNPID'}{'NecessaryFile'}{1}          = "prefix=AnnovarBuild\tsuffix=_avsnp150.txt";
    $hashAnnotationList{'SNPID'}{'GetResult'}                 = 'General SNP';
    $hashAnnotationList{'SNPID'}{'Output'}{1}{'File'}         = "prefix=outfile_library.AnnovarBuild\tsuffix=_avsnp150_dropped";
    $hashAnnotationList{'SNPID'}{'Output'}{1}{'Need_Info'}    = "SNP ID=1=2,3,4,5,6\t";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列

    $hashAnnotationList{'1000g'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'1000g'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'1000g'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype 1000g2015aug_all ';
    $hashAnnotationList{'1000g'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/1000g";
    $hashAnnotationList{'1000g'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/1000g";
    $hashAnnotationList{'1000g'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/1000g";
    $hashAnnotationList{'1000g'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'1000g'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_ALL.sites.2015_08_dropped";
    $hashAnnotationList{'1000g'}{'Output'}{1}{'Need_Info'}   = "Freq_Alt (1000g)=1=2,3,4,5,6\t";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列

    $hashAnnotationList{'Homology'}{'Species'}                = '[all]';
    $hashAnnotationList{'Homology'}{'From'}                   = 'File';
    $hashAnnotationList{'Homology'}{'NecessaryFile'}{1}       = "prefix=GenomeFasta\tsuffix=";
    $hashAnnotationList{'Homology'}{'GetResult'}              = 'General SNP';
    $hashAnnotationList{'Homology'}{'Output'}{1}{'File'}      = "prefix=outfile_library\tsuffix=.homology.anno";
    $hashAnnotationList{'Homology'}{'Output'}{1}{'Need_Info'} = "5'Flanking Sequence=5=0,1,2,3,4\t3'Flanking Sequence=6=0,1,2,3,4\tHomology Hits=7=0,1,2,3,4";

    $hashAnnotationList{'STR'}{'Species'}                     = '[all]';
    $hashAnnotationList{'STR'}{'From'}                        = 'File';
    $hashAnnotationList{'STR'}{'NecessaryFile'}{1}            = "prefix=GenomeFasta\tsuffix=";
    $hashAnnotationList{'STR'}{'GetResult'}                   = 'General SNP';
    $hashAnnotationList{'STR'}{'Output'}{1}{'File'}           = "prefix=outfile_library\tsuffix=.str.anno";
    $hashAnnotationList{'STR'}{'Output'}{1}{'Need_Info'}      = "STR Check=5=0,1,2,3,4";

    $hashAnnotationList{'VCF_FILTER'}{'Species'}                 = '[all]';
    $hashAnnotationList{'VCF_FILTER'}{'From'}                    = 'File';
    $hashAnnotationList{'VCF_FILTER'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'VCF_FILTER'}{'Output'}{1}{'File'}       = "prefix=outfile_library\tsuffix=.vcf.filter";
    $hashAnnotationList{'VCF_FILTER'}{'Output'}{1}{'Need_Info'}  = "VCF Filter=5=0,1,2,3,4";
    
    $hashAnnotationList{'Intervar'}{'Species'}                   = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'Intervar'}{'From'}                      = 'Intervar';
    $hashAnnotationList{'Intervar'}{'GetResult'}                 = 'General SNP';
    $hashAnnotationList{'Intervar'}{'Output'}{1}{'File'}         = "prefix=outfile_library\tsuffix=.intervar.ACMG";
    $hashAnnotationList{'Intervar'}{'Output'}{1}{'Need_Info'}    = "InterVar=5=0,1,2,3,4\tInterVar_evidence=6=0,1,2,3,4"; 
    
    ###
    # 频率注释
    ###
    $hashAnnotationList{'ExAC03'}{'Species'}                     = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'ExAC03'}{'From'}                        = 'Annovar';
    $hashAnnotationList{'ExAC03'}{'Para'}                        = 'annotate_variation.pl -filter -dbtype exac03 --otherinfo';
    $hashAnnotationList{'ExAC03'}{'DBDir'}{'Human_hg19'}         = "/home/genesky/database/annovar/hg19/exac03";
    $hashAnnotationList{'ExAC03'}{'DBDir'}{'Human_hg19_HLA'}     = "/home/genesky/database/annovar/hg19/exac03";
    $hashAnnotationList{'ExAC03'}{'DBDir'}{'Human_hg38'}         = "/home/genesky/database/annovar/hg38/exac03";
    $hashAnnotationList{'ExAC03'}{'NecessaryFile'}{1}            = "prefix=AnnovarBuild\tsuffix=_exac03.txt";
    $hashAnnotationList{'ExAC03'}{'GetResult'}                   = 'General SNP';
    $hashAnnotationList{'ExAC03'}{'Output'}{1}{'File'}           = "prefix=outfile_library.AnnovarBuild\tsuffix=_exac03_dropped";
    $hashAnnotationList{'ExAC03'}{'Output'}{1}{'Need_Info'}      = "ExAC03=1/,/0=2,3,4,5,6\tExAC03_EAS=1/,/3=2,3,4,5,6\t";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列

    $hashAnnotationList{'1000g_chbs'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'1000g_chbs'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'1000g_chbs'}{'Para'}                    = 'annotate_variation.pl -filter -dbtype CHBS1000g_130502';
    $hashAnnotationList{'1000g_chbs'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/chbs1000g";
    $hashAnnotationList{'1000g_chbs'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/annovar/hg19/chbs1000g";
    $hashAnnotationList{'1000g_chbs'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_CHBS1000g_130502.txt";
    $hashAnnotationList{'1000g_chbs'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'1000g_chbs'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_CHBS1000g_130502_dropped";
    $hashAnnotationList{'1000g_chbs'}{'Output'}{1}{'Need_Info'}  = "1000g_chbs=1=2,3,4,5,6"; 

    $hashAnnotationList{'esp6500'}{'Species'}                    = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'esp6500'}{'From'}                       = 'Annovar';
    $hashAnnotationList{'esp6500'}{'Para'}                       = 'annotate_variation.pl -filter -dbtype esp6500siv2_all';
    $hashAnnotationList{'esp6500'}{'DBDir'}{'Human_hg19'}        = "/home/genesky/database/annovar/hg19/esp6500siv2_all";
    $hashAnnotationList{'esp6500'}{'DBDir'}{'Human_hg19_HLA'}    = "/home/genesky/database/annovar/hg19/esp6500siv2_all";
    $hashAnnotationList{'esp6500'}{'DBDir'}{'Human_hg38'}        = "/home/genesky/database/annovar/hg38/esp6500siv2_all";
    $hashAnnotationList{'esp6500'}{'NecessaryFile'}{1}           = "prefix=AnnovarBuild\tsuffix=_esp6500siv2_all.txt";
    $hashAnnotationList{'esp6500'}{'GetResult'}                  = 'General SNP';
    $hashAnnotationList{'esp6500'}{'Output'}{1}{'File'}          = "prefix=outfile_library.AnnovarBuild\tsuffix=_esp6500siv2_all_dropped";
    $hashAnnotationList{'esp6500'}{'Output'}{1}{'Need_Info'}     = "esp6500=1=2,3,4,5,6";# 注释名称=注释内容在第几列（0开头）= SNP信息在第几列

    $hashAnnotationList{'Hrcr1'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'Hrcr1'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'Hrcr1'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype hrcr1';
    $hashAnnotationList{'Hrcr1'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/hrcr1";
    $hashAnnotationList{'Hrcr1'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/hrcr1";
    $hashAnnotationList{'Hrcr1'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/hrcr1";
    $hashAnnotationList{'Hrcr1'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_hrcr1.txt";
    $hashAnnotationList{'Hrcr1'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'Hrcr1'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_hrcr1_dropped";
    $hashAnnotationList{'Hrcr1'}{'Output'}{1}{'Need_Info'}   = "Hrcr1=1=2,3,4,5,6"; 

    $hashAnnotationList{'Kaviar'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'Kaviar'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'Kaviar'}{'Para'}                    = 'annotate_variation.pl -filter -dbtype kaviar_20150923';
    $hashAnnotationList{'Kaviar'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/kaviar";
    $hashAnnotationList{'Kaviar'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/annovar/hg19/kaviar";
    $hashAnnotationList{'Kaviar'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/annovar/hg38/kaviar";
    $hashAnnotationList{'Kaviar'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_kaviar_20150923.txt";
    $hashAnnotationList{'Kaviar'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'Kaviar'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_kaviar_20150923_dropped";
    $hashAnnotationList{'Kaviar'}{'Output'}{1}{'Need_Info'}  = "Kaviar=1=2,3,4,5,6"; 

    $hashAnnotationList{'NCI60'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'NCI60'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'NCI60'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype nci60';
    $hashAnnotationList{'NCI60'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/nci60";
    $hashAnnotationList{'NCI60'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/nci60";
    $hashAnnotationList{'NCI60'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/nci60";
    $hashAnnotationList{'NCI60'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_nci60.txt";
    $hashAnnotationList{'NCI60'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'NCI60'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_nci60_dropped";
    $hashAnnotationList{'NCI60'}{'Output'}{1}{'Need_Info'}   = "NCI60=1=2,3,4,5,6";

    $hashAnnotationList{'gnomAD_exome'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'gnomAD_exome'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'gnomAD_exome'}{'Para'}                    = 'annotate_variation.pl -filter -dbtype gnomad2.1_exome --otherinfo';
    $hashAnnotationList{'gnomAD_exome'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/gnomad_exome";
    $hashAnnotationList{'gnomAD_exome'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/annovar/hg19/gnomad_exome";
    $hashAnnotationList{'gnomAD_exome'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/annovar/hg38/gnomad_exome";
    $hashAnnotationList{'gnomAD_exome'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_gnomad2.1_exome.txt";
    $hashAnnotationList{'gnomAD_exome'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'gnomAD_exome'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_gnomad2.1_exome_dropped";
    $hashAnnotationList{'gnomAD_exome'}{'Output'}{1}{'Need_Info'}  = "gnomAD_exome=1/,/0=2,3,4,5,6\tgnomAD_exome_EAS=1/,/4=2,3,4,5,6\t";

    $hashAnnotationList{'gnomAD_genome'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'gnomAD_genome'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'gnomAD_genome'}{'Para'}                    = 'annotate_variation.pl -filter -dbtype gnomad211_genome --otherinfo';
    $hashAnnotationList{'gnomAD_genome'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/gnomad_genome";
    $hashAnnotationList{'gnomAD_genome'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/annovar/hg19/gnomad_genome";
    $hashAnnotationList{'gnomAD_genome'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/annovar/hg38/gnomad_genome";
    $hashAnnotationList{'gnomAD_genome'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_gnomad211_genome.txt";
    $hashAnnotationList{'gnomAD_genome'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'gnomAD_genome'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_gnomad211_genome_dropped";
    $hashAnnotationList{'gnomAD_genome'}{'Output'}{1}{'Need_Info'}  = "gnomAD_genome=1/,/0=2,3,4,5,6\tgnomAD_genome_EAS=1/,/8=2,3,4,5,6\t";

    $hashAnnotationList{'GeneskyExonDB'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'GeneskyExonDB'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'GeneskyExonDB'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype geneskyExon_300_GATK4 -otherinfo';
    $hashAnnotationList{'GeneskyExonDB'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/genesky_exon";
    $hashAnnotationList{'GeneskyExonDB'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/genesky_exon";
    $hashAnnotationList{'GeneskyExonDB'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_geneskyExon_300_GATK4.txt";
    $hashAnnotationList{'GeneskyExonDB'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'GeneskyExonDB'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_geneskyExon_300_GATK4_dropped";
    $hashAnnotationList{'GeneskyExonDB'}{'Output'}{1}{'Need_Info'}   = "GeneskyExonDB_Freq=1/,/0=2,3,4,5,6\tGeneskyExonDB_SampeCount=1/,/2=2,3,4,5,6";
    
    $hashAnnotationList{'GeneskyGenomeDB'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'GeneskyGenomeDB'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'GeneskyGenomeDB'}{'Para'}                    = 'annotate_variation.pl -filter -dbtype geneskyGenome_128_GATK -otherinfo';
    $hashAnnotationList{'GeneskyGenomeDB'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/genesky_genome";
    $hashAnnotationList{'GeneskyGenomeDB'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/annovar/hg19/genesky_genome";
    $hashAnnotationList{'GeneskyGenomeDB'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_geneskyGenome_128_GATK.txt";
    $hashAnnotationList{'GeneskyGenomeDB'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'GeneskyGenomeDB'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_geneskyGenome_128_GATK_dropped";
    $hashAnnotationList{'GeneskyGenomeDB'}{'Output'}{1}{'Need_Info'}  = "GeneskyGenomeDB_Freq=1/,/0=2,3,4,5,6\tGeneskyGenomeDB_SampeCount=1/,/2=2,3,4,5,6";
    
    ###
    # 突变危险性注释
    ###     
    $hashAnnotationList{'DBNSFP'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'DBNSFP'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'DBNSFP'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype dbnsfp35a --otherinfo';
    $hashAnnotationList{'DBNSFP'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/dbnsfp";
    $hashAnnotationList{'DBNSFP'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/dbnsfp";
    $hashAnnotationList{'DBNSFP'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/dbnsfp";
    $hashAnnotationList{'DBNSFP'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_dbnsfp35a.txt";
    $hashAnnotationList{'DBNSFP'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'DBNSFP'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_dbnsfp35a_dropped";
    $hashAnnotationList{'DBNSFP'}{'Output'}{1}{'Need_Info'}   = "SIFT Score=1/,/0=2,3,4,5,6\tSIFT Score Pred=1/,/2=2,3,4,5,6\tPOLYPhen V2 Score=1/,/3=2,3,4,5,6\tPOLYPhen V2 Score Pred=1/,/5=2,3,4,5,6\tMutationTaster=1/,/12=2,3,4,5,6\tMutationTaster Pred=1/,/14=2,3,4,5,6\t";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列

    $hashAnnotationList{'CADD'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'CADD'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'CADD'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype cadd13 -otherinfo ';
    $hashAnnotationList{'CADD'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/cadd";
    $hashAnnotationList{'CADD'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/cadd";
    $hashAnnotationList{'CADD'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/cadd";
    $hashAnnotationList{'CADD'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_cadd13.txt";
    $hashAnnotationList{'CADD'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'CADD'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_cadd13_dropped";
    $hashAnnotationList{'CADD'}{'Output'}{1}{'Need_Info'}   = "Cadd_Raw=1/,/0=2,3,4,5,6\tCadd_Phred=1/,/1=2,3,4,5,6\t";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列
 
    $hashAnnotationList{'Dann'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'Dann'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'Dann'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype dann ';
    $hashAnnotationList{'Dann'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/dann";
    $hashAnnotationList{'Dann'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/dann";
    $hashAnnotationList{'Dann'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/dann";
    $hashAnnotationList{'Dann'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_dann.txt";
    $hashAnnotationList{'Dann'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'Dann'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_dann_dropped";
    $hashAnnotationList{'Dann'}{'Output'}{1}{'Need_Info'}   = "Dann=1=2,3,4,5,6";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列
 
    $hashAnnotationList{'Eigen'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'Eigen'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'Eigen'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype eigen ';
    $hashAnnotationList{'Eigen'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/eigen";
    $hashAnnotationList{'Eigen'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/eigen";
    $hashAnnotationList{'Eigen'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/eigen";
    $hashAnnotationList{'Eigen'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_eigen.txt";
    $hashAnnotationList{'Eigen'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'Eigen'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_eigen_dropped";
    $hashAnnotationList{'Eigen'}{'Output'}{1}{'Need_Info'}   = "Eigen=1=2,3,4,5,6";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列

    $hashAnnotationList{'dbscSNV'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'dbscSNV'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'dbscSNV'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype dbscsnv11 -otherinfo';
    $hashAnnotationList{'dbscSNV'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/dbscsnv";
    $hashAnnotationList{'dbscSNV'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/dbscsnv";
    $hashAnnotationList{'dbscSNV'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/dbscsnv";
    $hashAnnotationList{'dbscSNV'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_dbscsnv11.txt";
    $hashAnnotationList{'dbscSNV'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'dbscSNV'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_dbscsnv11_dropped";
    $hashAnnotationList{'dbscSNV'}{'Output'}{1}{'Need_Info'}   = "dbscSNV_ADA_SCORE=1/,/0=2,3,4,5,6\tdbscSNV_RF_SCORE=1/,/1=2,3,4,5,6"; 

    ###
    # 突变位点疾病数据库
    ### 
    $hashAnnotationList{'HGMD'}{'Species'}                   = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'HGMD'}{'From'}                      = 'Annovar';
    $hashAnnotationList{'HGMD'}{'Para'}                      = 'annotate_variation.pl -filter -dbtype HGMD_Site_VariantClass -otherinfo';
    $hashAnnotationList{'HGMD'}{'DBDir'}{'Human_hg19'}       = "/home/genesky/database/annovar/hg19/hgmd";
    $hashAnnotationList{'HGMD'}{'DBDir'}{'Human_hg19_HLA'}   = "/home/genesky/database/annovar/hg19/hgmd";
    $hashAnnotationList{'HGMD'}{'DBDir'}{'Human_hg38'}       = "/home/genesky/database/annovar/hg38/hgmd";
    $hashAnnotationList{'HGMD'}{'NecessaryFile'}{1}          = "prefix=AnnovarBuild\tsuffix=_HGMD_Site_VariantClass.txt";
    $hashAnnotationList{'HGMD'}{'GetResult'}                 = 'General SNP';
    $hashAnnotationList{'HGMD'}{'Output'}{1}{'File'}         = "prefix=outfile_library.AnnovarBuild\tsuffix=_HGMD_Site_VariantClass_dropped";
    $hashAnnotationList{'HGMD'}{'Output'}{1}{'Need_Info'}    = "HGMD_site=1/,/1,2=2,3,4,5,6\tHGMD_site_class=1/,/0=2,3,4,5,6";# 注释名称=注释内容在第几列（0开头）=SNP信息在第几列

    $hashAnnotationList{'ClinVar'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'ClinVar'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'ClinVar'}{'Para'}                    = 'annotate_variation.pl -filter -dbtype clinvar_20180603 -otherinfo';
    $hashAnnotationList{'ClinVar'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/clinvar";
    $hashAnnotationList{'ClinVar'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/annovar/hg19/clinvar";
    $hashAnnotationList{'ClinVar'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/annovar/hg38/clinvar";
    $hashAnnotationList{'ClinVar'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_clinvar_20180603.txt";
    $hashAnnotationList{'ClinVar'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'ClinVar'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_clinvar_20180603_dropped";
    $hashAnnotationList{'ClinVar'}{'Output'}{1}{'Need_Info'}  = "ClinVar=1=2,3,4,5,6";

    $hashAnnotationList{'COSMIC'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'COSMIC'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'COSMIC'}{'Para'}                    = 'annotate_variation.pl -filter -dbtype cosmic70';
    $hashAnnotationList{'COSMIC'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/cosmic";
    $hashAnnotationList{'COSMIC'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/annovar/hg19/cosmic";
    $hashAnnotationList{'COSMIC'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/annovar/hg38/cosmic";
    $hashAnnotationList{'COSMIC'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_cosmic70.txt";
    $hashAnnotationList{'COSMIC'}{'GetResult'}               = 'General SNP';
    $hashAnnotationList{'COSMIC'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_cosmic70_dropped";
    $hashAnnotationList{'COSMIC'}{'Output'}{1}{'Need_Info'}  = "COSMIC=1=2,3,4,5,6";
    
    $hashAnnotationList{'ICGC'}{'Species'}                   = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'ICGC'}{'From'}                      = 'Annovar';
    $hashAnnotationList{'ICGC'}{'Para'}                      = 'annotate_variation.pl -filter -dbtype icgc21 -otherinfo';
    $hashAnnotationList{'ICGC'}{'DBDir'}{'Human_hg19'}       = "/home/genesky/database/annovar/hg19/icgc";
    $hashAnnotationList{'ICGC'}{'DBDir'}{'Human_hg19_HLA'}   = "/home/genesky/database/annovar/hg19/icgc";
    $hashAnnotationList{'ICGC'}{'NecessaryFile'}{1}          = "prefix=AnnovarBuild\tsuffix=_icgc21.txt";
    $hashAnnotationList{'ICGC'}{'GetResult'}                 = 'General SNP';
    $hashAnnotationList{'ICGC'}{'Output'}{1}{'File'}         = "prefix=outfile_library.AnnovarBuild\tsuffix=_icgc21_dropped";
    $hashAnnotationList{'ICGC'}{'Output'}{1}{'Need_Info'}    = "ICGC=1=2,3,4,5,6";
    
    ###
    # 区域注释
    ###
    $hashAnnotationList{'tfbsConsSites'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'tfbsConsSites'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'tfbsConsSites'}{'Para'}                     = 'annotate_variation.pl -regionanno -dbtype tfbsConsSites';
    $hashAnnotationList{'tfbsConsSites'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/tfbsconssites";
    $hashAnnotationList{'tfbsConsSites'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/tfbsconssites";
    $hashAnnotationList{'tfbsConsSites'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_tfbsConsSites.txt";
    $hashAnnotationList{'tfbsConsSites'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'tfbsConsSites'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_tfbsConsSites";
    $hashAnnotationList{'tfbsConsSites'}{'Output'}{1}{'Need_Info'}   = "tfbsConsSites Score=1=2,3,4,5,6";

    ###
    # 突变蛋白注释
    ###
    $hashAnnotationList{'Interpro'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'Interpro'}{'From'}                     = 'Annovar';
    $hashAnnotationList{'Interpro'}{'Para'}                     = 'annotate_variation.pl -filter -dbtype dbnsfp31a_interpro';
    $hashAnnotationList{'Interpro'}{'DBDir'}{'Human_hg19'}      = "/home/genesky/database/annovar/hg19/dbnsfp_interpro";
    $hashAnnotationList{'Interpro'}{'DBDir'}{'Human_hg19_HLA'}  = "/home/genesky/database/annovar/hg19/dbnsfp_interpro";
    $hashAnnotationList{'Interpro'}{'DBDir'}{'Human_hg38'}      = "/home/genesky/database/annovar/hg38/dbnsfp_interpro";
    $hashAnnotationList{'Interpro'}{'NecessaryFile'}{1}         = "prefix=AnnovarBuild\tsuffix=_dbnsfp31a_interpro.txt";
    $hashAnnotationList{'Interpro'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'Interpro'}{'Output'}{1}{'File'}        = "prefix=outfile_library.AnnovarBuild\tsuffix=_dbnsfp31a_interpro_dropped";
    $hashAnnotationList{'Interpro'}{'Output'}{1}{'Need_Info'}   = "Interpro_domain=1=2,3,4,5,6";
    
    ###
    # 基因注释
    ###
    $hashAnnotationList{'GeneStrand'}{'Species'}                 = '[all]';
    $hashAnnotationList{'GeneStrand'}{'From'}                    = 'File';
    $hashAnnotationList{'GeneStrand'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/ucsc/hg19/gene";
    $hashAnnotationList{'GeneStrand'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/ucsc/hg19/gene";
    $hashAnnotationList{'GeneStrand'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/ucsc/hg38/gene";
    $hashAnnotationList{'GeneStrand'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_refGene.txt";
    $hashAnnotationList{'GeneStrand'}{'GetResult'}               = 'General Gene';
    $hashAnnotationList{'GeneStrand'}{'Output'}{1}{'Need_Info'}  = "Gene Strand Orientation=3=12";
    
    $hashAnnotationList{'HGMD_gene'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'HGMD_gene'}{'From'}                    = 'File';
    $hashAnnotationList{'HGMD_gene'}{'NecessaryFile'}{1}        = "prefix=\tsuffix=/home/genesky/database/annovar/hg19/hgmd/hg19_HGMD_Gene.txt";
    $hashAnnotationList{'HGMD_gene'}{'GetResult'}               = 'General Gene';
    $hashAnnotationList{'HGMD_gene'}{'Output'}{1}{'Need_Info'}  = "HGMD_gene=1=0";

    $hashAnnotationList{'Gene_Info'}{'Species'}                 = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'Gene_Info'}{'From'}                    = 'File';
    $hashAnnotationList{'Gene_Info'}{'NecessaryFile'}{1}        = "prefix=\tsuffix=/home/genesky/database/ncbi/gene_info/Homo_sapiens.gene_info";
    $hashAnnotationList{'Gene_Info'}{'GetResult'}               = 'General Gene';
    $hashAnnotationList{'Gene_Info'}{'Output'}{1}{'Need_Info'}  = "Gene_Info=12=1";

    $hashAnnotationList{'Gene_Info_Mouse'}{'Species'}                = '[Mouse_10]';
    $hashAnnotationList{'Gene_Info_Mouse'}{'From'}                   = 'File';
    $hashAnnotationList{'Gene_Info_Mouse'}{'NecessaryFile'}{1}       = "prefix=\tsuffix=/home/genesky/database/ncbi/gene_info/mm10.gene_info";
    $hashAnnotationList{'Gene_Info_Mouse'}{'GetResult'}              = 'General Gene';
    $hashAnnotationList{'Gene_Info_Mouse'}{'Output'}{1}{'Need_Info'} = "Gene_Info=12=1";

    $hashAnnotationList{'OMIM'}{'Species'}                     = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'OMIM'}{'From'}                        = 'File';
    $hashAnnotationList{'OMIM'}{'NecessaryFile'}{1}            = "prefix=\tsuffix=/home/genesky/database/self_build_database/human_gene/omim/gene2omim_20180921.txt";
    $hashAnnotationList{'OMIM'}{'GetResult'}                   = 'General Gene';
    $hashAnnotationList{'OMIM'}{'Output'}{1}{'Need_Info'}      = "OMIM=1=0";

    $hashAnnotationList{'MalaCards'}{'Species'}                = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'MalaCards'}{'From'}                   = 'File';
    $hashAnnotationList{'MalaCards'}{'NecessaryFile'}{1}       = "prefix=\tsuffix=/home/genesky/database/self_build_database/human_gene/malacards/MalaCards_20180921.txt";
    $hashAnnotationList{'MalaCards'}{'GetResult'}              = 'General Gene';
    $hashAnnotationList{'MalaCards'}{'Output'}{1}{'Need_Info'} = "MalaCards=1=0";

    $hashAnnotationList{'MGI'}{'Species'}                      = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'MGI'}{'From'}                         = 'File';
    $hashAnnotationList{'MGI'}{'NecessaryFile'}{1}             = "prefix=\tsuffix=/home/genesky/database/self_build_database/human_gene/mgi/MGI.xls";
    $hashAnnotationList{'MGI'}{'GetResult'}                    = 'General Gene';
    $hashAnnotationList{'MGI'}{'Output'}{1}{'Need_Info'}       = "MGI=1=0";

    $hashAnnotationList{'MGI_Mouse'}{'Species'}                = '[Mouse_10]';
    $hashAnnotationList{'MGI_Mouse'}{'From'}                   = 'File';
    $hashAnnotationList{'MGI_Mouse'}{'NecessaryFile'}{1}       = "prefix=\tsuffix=/home/genesky/database/self_build_database/mouse_gene/mgi/MGI.txt";
    $hashAnnotationList{'MGI_Mouse'}{'GetResult'}              = 'General Gene';
    $hashAnnotationList{'MGI_Mouse'}{'Output'}{1}{'Need_Info'} = "MGI=1=0";

    $hashAnnotationList{'HPO'}{'Species'}                      = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'HPO'}{'From'}                         = 'File';
    $hashAnnotationList{'HPO'}{'NecessaryFile'}{1}             = "prefix=\tsuffix=/home/genesky/database/self_build_database/human_gene/hpo/HPO.txt";
    $hashAnnotationList{'HPO'}{'GetResult'}                    = 'General Gene';
    $hashAnnotationList{'HPO'}{'Output'}{1}{'Need_Info'}       = "HPO=1=0";

    $hashAnnotationList{'GO'}{'Species'}                       = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'GO'}{'From'}                          = 'File';
    $hashAnnotationList{'GO'}{'NecessaryFile'}{1}              = "prefix=\tsuffix=/home/genesky/database/self_build_database/human_gene/go_kegg/go_info";
    $hashAnnotationList{'GO'}{'GetResult'}                     = 'General Gene';
    $hashAnnotationList{'GO'}{'Output'}{1}{'Need_Info'}        = "GO_BP=1=0\tGO_MF=2=0\tGO_CC=3=0\tKEGG_Pathway=4=0";

    $hashAnnotationList{'GO_Mouse'}{'Species'}                 = '[Mouse_10]';
    $hashAnnotationList{'GO_Mouse'}{'From'}                    = 'File';
    $hashAnnotationList{'GO_Mouse'}{'NecessaryFile'}{1}        = "prefix=\tsuffix=/home/genesky/database/self_build_database/mouse_gene/go_kegg/mm10.go_info";
    $hashAnnotationList{'GO_Mouse'}{'GetResult'}               = 'General Gene';
    $hashAnnotationList{'GO_Mouse'}{'Output'}{1}{'Need_Info'}  = "GO_BP=1=0\tGO_MF=2=0\tGO_CC=3=0";

    $hashAnnotationList{'GeneskyExonDB_lowfreq_geneCount'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'GeneskyExonDB_lowfreq_geneCount'}{'From'}                     = 'File';
    $hashAnnotationList{'GeneskyExonDB_lowfreq_geneCount'}{'NecessaryFile'}{1}         = "prefix=\tsuffix=/home/genesky/database/annovar/hg19/genesky_exon/hg19_geneskyExon_300_GATK4_lowfreq_geneCount_20190419.txt";
    $hashAnnotationList{'GeneskyExonDB_lowfreq_geneCount'}{'GetResult'}                = 'General Gene';
    $hashAnnotationList{'GeneskyExonDB_lowfreq_geneCount'}{'Output'}{1}{'Need_Info'}   = "GeneskyExonDB SNV Count(normal)=1/,/0=0\tGeneskyExonDB Mutation(0|1|2)(normal)=1/,/1=0\tGeneskyExonDB SNV Count(normal_pass)=2/,/0=0\tGeneskyExonDB Mutation(0|1|2)(normal_pass)=2/,/1=0\tGeneskyExonDB SNV Count(strict)=3/,/0=0\tGeneskyExonDB Mutation(0|1|2)(strict)=3/,/1=0\tGeneskyExonDB SNV Count(strict_pass)=4/,/0=0\tGeneskyExonDB Mutation(0|1|2)(strict_pass)=4/,/1=0";

    $hashAnnotationList{'GeneskyGenomeDB_lowfreq_geneCount'}{'Species'}                = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'GeneskyGenomeDB_lowfreq_geneCount'}{'From'}                   = 'File';
    $hashAnnotationList{'GeneskyGenomeDB_lowfreq_geneCount'}{'NecessaryFile'}{1}       = "prefix=\tsuffix=/home/genesky/database/annovar/hg19/genesky_genome/hg19_geneskyGenome_128_GATK_lowfreq_geneCount_20190423.txt";
    $hashAnnotationList{'GeneskyGenomeDB_lowfreq_geneCount'}{'GetResult'}              = 'General Gene';
    $hashAnnotationList{'GeneskyGenomeDB_lowfreq_geneCount'}{'Output'}{1}{'Need_Info'} = "GeneskyGenomeDB SNV Count(normal)=1/,/0=0\tGeneskyGenomeDB Mutation(0|1|2)(normal)=1/,/1=0\tGeneskyGenomeDB SNV Count(normal_pass)=2/,/0=0\tGeneskyGenomeDB Mutation(0|1|2)(normal_pass)=2/,/1=0\tGeneskyGenomeDB SNV Count(strict)=3/,/0=0\tGeneskyGenomeDB Mutation(0|1|2)(strict)=3/,/1=0\tGeneskyGenomeDB SNV Count(strict_pass)=4/,/0=0\tGeneskyGenomeDB Mutation(0|1|2)(strict_pass)=4/,/1=0";

    ###
    # 其他注释
    ###
    $hashAnnotationList{'snpEFF'}{'Species'}                  = '[Human_hg19][Human_hg38][Human_hg19_HLA]';
    $hashAnnotationList{'snpEFF'}{'From'}                     = 'snpEFF';
    $hashAnnotationList{'snpEFF'}{'GetResult'}                = 'General SNP';
    $hashAnnotationList{'snpEFF'}{'Output'}{1}{'File'}        = "prefix=outfile_library\tsuffix=.vcf.snpEFF.anno.final";
    $hashAnnotationList{'snpEFF'}{'Output'}{1}{'Need_Info'}   = "Predicted Protein Variants(snpEFF)=5=0,1,2,3,4"; 
 
    return %hashAnnotationList;    
}

1