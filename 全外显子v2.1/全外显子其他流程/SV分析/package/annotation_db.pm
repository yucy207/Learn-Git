package package::annotation_db;
use strict;
use warnings;

# 获取注释关键字
sub get_annotation_key_word{
    my $hashPara     = shift @_;
    my $hashConfig   = shift @_;
    my $species      = package::utils::get_species($hashConfig);
    my $AnnovarBuild = $hashPara->{$species}{'AnnovarBuild'};
    my $AnnovarDB    = $hashPara->{$species}{'AnnovarDB'};    
    my $report_dir   = $hashConfig->{'Report'};

    my %hashKeyWord;
    $hashKeyWord{'AnnovarBuild'}                   = "$AnnovarBuild";
    $hashKeyWord{'outfile_library'}                = "$report_dir/sv/library/library";
    $hashKeyWord{'outfile_library.AnnovarBuild'}   = "$report_dir/sv/library/library.$AnnovarBuild";
    $hashKeyWord{'GenomeFasta'}                    = $hashPara->{$species}{'Genome'};
    return %hashKeyWord;
}

sub get_annotation{
    my %hashAnnotationList;
    # 关键字说明
    # （1）Species       : 该注释对应的物种，用"[物种名]"表示,如果对应多个物种，用","分割
    # （2）From          : 当前注释的来源方式，Annovar表示需要用annovar进行注释（后面的Para必须提供）;File表示需要写代码从文件提取
    # （3）Para          : Annovar注释时，该参数提供注释参数
    # （4）DBDir         : 当前注释文件所在路径
    # （5）NecessaryFile : 注释时必须存在的文件，由"prefix=...;suffix=..."构成,关键字有"AnnovarDB/AnnovarBuild",代码会识别该关键字进行相应路径替换(部分注释会需要多个依赖文件，用不同的数字编号key区分)
    # （6）GetResult     : 结果获取方式，"General Pos"表示通用SNP处理代码即可处理
    # （7）Output:File   : 结果输出文件，由"prefix=...;suffix=..."构成,关键字有"outfile.AnnovarBuild"
    # （8）Output:Need_Info: 该注释需要的信息，（1）General Pos 处理模式：名称=信息位置
 
    ###
    # 基础注释
    ###    
    $hashAnnotationList{'refGene'}{'Species'}                 = '[all]';
    $hashAnnotationList{'refGene'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'refGene'}{'Para'}                    = 'annotate_variation.pl --hgvs --splicing_threshold 15 ';
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/ucsc/hg19/gene";
    # $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/annovar/hg19/refgene";
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19_HLA'} = "/home/genesky/database/ucsc/hg19/gene";
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg38'}     = "/home/genesky/database/ucsc/hg38/gene";
    $hashAnnotationList{'refGene'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_refGene.txt";
    $hashAnnotationList{'refGene'}{'NecessaryFile'}{2}        = "prefix=AnnovarBuild\tsuffix=_refGeneMrna.fa";
    $hashAnnotationList{'refGene'}{'GetResult'}               = 'General Pos';
    $hashAnnotationList{'refGene'}{'Output'}{1}{'File'}       = "prefix=outfile_library\tsuffix=.variant_function";# 该注释要包含的注释内容
    $hashAnnotationList{'refGene'}{'Output'}{1}{'Need_Info'}  = "Gene Region=0=7,\tGene=1=7,";# 该注释要包含的注释内容

    ###
    # 区域注释
    ###
    $hashAnnotationList{'cytoBand'}{'Species'}                = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'cytoBand'}{'From'}                   = 'Annovar';
    $hashAnnotationList{'cytoBand'}{'Para'}                   = 'annotate_variation.pl -regionanno -dbtype cytoBand';
    $hashAnnotationList{'cytoBand'}{'DBDir'}{'Human_hg19'}    = "/home/genesky/database/self_build_database/human_annovar/hg19/cytoband";
    $hashAnnotationList{'cytoBand'}{'NecessaryFile'}{1}       = "prefix=AnnovarBuild\tsuffix=_cytoBand.txt";
    $hashAnnotationList{'cytoBand'}{'GetResult'}              = 'General Pos';
    $hashAnnotationList{'cytoBand'}{'Output'}{1}{'File'}      = "prefix=outfile_library.AnnovarBuild\tsuffix=_cytoBand";
    $hashAnnotationList{'cytoBand'}{'Output'}{1}{'Need_Info'} = "cytoBand=1=7,";


    $hashAnnotationList{'dgvFreq'}{'Species'}                 = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'dgvFreq'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'dgvFreq'}{'Para'}                    = 'annotate_variation.pl -regionanno -dbtype dgvFreq';
    $hashAnnotationList{'dgvFreq'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/self_build_database/human_annovar/hg19/dgvfreq";
    $hashAnnotationList{'dgvFreq'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_dgvFreq.txt";
    $hashAnnotationList{'dgvFreq'}{'GetResult'}               = 'General Pos';
    $hashAnnotationList{'dgvFreq'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_dgvFreq";
    $hashAnnotationList{'dgvFreq'}{'Output'}{1}{'Need_Info'}  = "dgvFreq=1=0,";

    $hashAnnotationList{'iscaPathGainCum'}{'Species'}                      = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'iscaPathGainCum'}{'From'}                         = 'Annovar';
    $hashAnnotationList{'iscaPathGainCum'}{'Para'}                         = 'annotate_variation.pl -regionanno -dbtype iscaPathGainCum';
    $hashAnnotationList{'iscaPathGainCum'}{'DBDir'}{'Human_hg19'}          = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_path_gain_cum";
    $hashAnnotationList{'iscaPathGainCum'}{'NecessaryFile'}{1}             = "prefix=AnnovarBuild\tsuffix=_iscaPathGainCum.txt";
    $hashAnnotationList{'iscaPathGainCum'}{'GetResult'}                    = 'General Pos';
    $hashAnnotationList{'iscaPathGainCum'}{'Output'}{1}{'File'}            = "prefix=outfile_library.AnnovarBuild\tsuffix=_iscaPathGainCum";
    $hashAnnotationList{'iscaPathGainCum'}{'Output'}{1}{'Need_Info'}       = "iscaPathGainCum=1=7,";
                                                                           
    $hashAnnotationList{'iscaPathLossCum'}{'Species'}                      = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'iscaPathLossCum'}{'From'}                         = 'Annovar';
    $hashAnnotationList{'iscaPathLossCum'}{'Para'}                         = 'annotate_variation.pl -regionanno -dbtype iscaPathLossCum';
    $hashAnnotationList{'iscaPathLossCum'}{'DBDir'}{'Human_hg19'}          = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_path_loss_cum";
    $hashAnnotationList{'iscaPathLossCum'}{'NecessaryFile'}{1}             = "prefix=AnnovarBuild\tsuffix=_iscaPathLossCum.txt";
    $hashAnnotationList{'iscaPathLossCum'}{'GetResult'}                    = 'General Pos';
    $hashAnnotationList{'iscaPathLossCum'}{'Output'}{1}{'File'}            = "prefix=outfile_library.AnnovarBuild\tsuffix=_iscaPathLossCum";
    $hashAnnotationList{'iscaPathLossCum'}{'Output'}{1}{'Need_Info'}       = "iscaPathLossCum=1=7,";  

    $hashAnnotationList{'iscaLikelyPathogenic'}{'Species'}                 = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'iscaLikelyPathogenic'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'iscaLikelyPathogenic'}{'Para'}                    = 'annotate_variation.pl -regionanno -dbtype iscaLikelyPathogenic';
    $hashAnnotationList{'iscaLikelyPathogenic'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_likely_pathogenic";
    $hashAnnotationList{'iscaLikelyPathogenic'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_iscaLikelyPathogenic.txt";
    $hashAnnotationList{'iscaLikelyPathogenic'}{'GetResult'}               = 'General Pos';
    $hashAnnotationList{'iscaLikelyPathogenic'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_iscaLikelyPathogenic";
    $hashAnnotationList{'iscaLikelyPathogenic'}{'Output'}{1}{'Need_Info'}  = "iscaLikelyPathogenic=1=7,";      
    
    $hashAnnotationList{'iscaPathogenic'}{'Species'}                       = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'iscaPathogenic'}{'From'}                          = 'Annovar';
    $hashAnnotationList{'iscaPathogenic'}{'Para'}                          = 'annotate_variation.pl -regionanno -dbtype iscaPathogenic';
    $hashAnnotationList{'iscaPathogenic'}{'DBDir'}{'Human_hg19'}           = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_pathogenic";
    $hashAnnotationList{'iscaPathogenic'}{'NecessaryFile'}{1}              = "prefix=AnnovarBuild\tsuffix=_iscaPathogenic.txt";
    $hashAnnotationList{'iscaPathogenic'}{'GetResult'}                     = 'General Pos';
    $hashAnnotationList{'iscaPathogenic'}{'Output'}{1}{'File'}             = "prefix=outfile_library.AnnovarBuild\tsuffix=_iscaPathogenic";
    $hashAnnotationList{'iscaPathogenic'}{'Output'}{1}{'Need_Info'}        = "iscaPathogenic=1=7,";      
       
    $hashAnnotationList{'iscaCuratedPathogenic'}{'Species'}                = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'iscaCuratedPathogenic'}{'From'}                   = 'Annovar';
    $hashAnnotationList{'iscaCuratedPathogenic'}{'Para'}                   = 'annotate_variation.pl -regionanno -dbtype iscaCuratedPathogenic';
    $hashAnnotationList{'iscaCuratedPathogenic'}{'DBDir'}{'Human_hg19'}    = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_curated_pathogenic";
    $hashAnnotationList{'iscaCuratedPathogenic'}{'NecessaryFile'}{1}       = "prefix=AnnovarBuild\tsuffix=_iscaCuratedPathogenic.txt";
    $hashAnnotationList{'iscaCuratedPathogenic'}{'GetResult'}              = 'General Pos';
    $hashAnnotationList{'iscaCuratedPathogenic'}{'Output'}{1}{'File'}      = "prefix=outfile_library.AnnovarBuild\tsuffix=_iscaCuratedPathogenic";
    $hashAnnotationList{'iscaCuratedPathogenic'}{'Output'}{1}{'Need_Info'} = "iscaCuratedPathogenic=1=7,";      
         
    $hashAnnotationList{'CNVD'}{'Species'}                     = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'CNVD'}{'From'}                        = 'Annovar';
    $hashAnnotationList{'CNVD'}{'Para'}                        = 'annotate_variation.pl -regionanno -dbtype CNVD';
    $hashAnnotationList{'CNVD'}{'DBDir'}{'Human_hg19'}         = "/home/genesky/database/self_build_database/human_annovar/hg19/cnvd";
    $hashAnnotationList{'CNVD'}{'NecessaryFile'}{1}            = "prefix=AnnovarBuild\tsuffix=_CNVD.txt";
    $hashAnnotationList{'CNVD'}{'GetResult'}                   = 'General Pos';
    $hashAnnotationList{'CNVD'}{'Output'}{1}{'File'}           = "prefix=outfile_library.AnnovarBuild\tsuffix=_CNVD";
    $hashAnnotationList{'CNVD'}{'Output'}{1}{'Need_Info'}      = "CNVD=1=7,";      

    $hashAnnotationList{'DECIPHER'}{'Species'}                 = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'DECIPHER'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'DECIPHER'}{'Para'}                    = 'annotate_variation.pl -regionanno -dbtype DECIPHER';
    $hashAnnotationList{'DECIPHER'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/self_build_database/human_annovar/hg19/decipher";
    $hashAnnotationList{'DECIPHER'}{'NecessaryFile'}{1}        = "prefix=AnnovarBuild\tsuffix=_DECIPHER.txt";
    $hashAnnotationList{'DECIPHER'}{'GetResult'}               = 'General Pos';
    $hashAnnotationList{'DECIPHER'}{'Output'}{1}{'File'}       = "prefix=outfile_library.AnnovarBuild\tsuffix=_DECIPHER";
    $hashAnnotationList{'DECIPHER'}{'Output'}{1}{'Need_Info'}  = "DECIPHER=1=7,";    
                                                                        
    $hashAnnotationList{'dbVar'}{'Species'}                    = '[Human_hg19][Human_hg19_HLA]';
    $hashAnnotationList{'dbVar'}{'From'}                       = 'Annovar';
    $hashAnnotationList{'dbVar'}{'Para'}                       = 'annotate_variation.pl -regionanno -dbtype dbVar';
    $hashAnnotationList{'dbVar'}{'DBDir'}{'Human_hg19'}        = "/home/genesky/database/self_build_database/human_annovar/hg19/dbvar";
    $hashAnnotationList{'dbVar'}{'NecessaryFile'}{1}           = "prefix=AnnovarBuild\tsuffix=_dbVar.txt";
    $hashAnnotationList{'dbVar'}{'GetResult'}                  = 'General Pos';
    $hashAnnotationList{'dbVar'}{'Output'}{1}{'File'}          = "prefix=outfile_library.AnnovarBuild\tsuffix=_dbVar";
    $hashAnnotationList{'dbVar'}{'Output'}{1}{'Need_Info'}     = "dbVar=1=0,";      

    # gene 注释
    $hashAnnotationList{'OMIM'}{'Species'}                     = '[Human_hg19][Human_hg19_rmHAP][Human_hg38][Human_hg19_HLA]';                       
    $hashAnnotationList{'OMIM'}{'From'}                        = 'File';                       
    $hashAnnotationList{'OMIM'}{'NecessaryFile'}{1}            = "prefix=\tsuffix=/home/genesky/database/self_build_database/human_gene/omim/gene2omim_20180921.txt";  
    $hashAnnotationList{'OMIM'}{'GetResult'}                   = 'General Gene';                                             
    $hashAnnotationList{'OMIM'}{'Output'}{1}{'Need_Info'}      = "OMIM=1=0"; 

    $hashAnnotationList{'MalaCards'}{'Species'}                = '[Human_hg19][Human_hg19_rmHAP][Human_hg38][Human_hg19_HLA]';                       
    $hashAnnotationList{'MalaCards'}{'From'}                   = 'File';                       
    $hashAnnotationList{'MalaCards'}{'NecessaryFile'}{1}       = "prefix=\tsuffix=/home/genesky/database/self_build_database/human_gene/malacards/MalaCards_20180921.txt";
    $hashAnnotationList{'MalaCards'}{'GetResult'}              = 'General Gene';                                             
    $hashAnnotationList{'MalaCards'}{'Output'}{1}{'Need_Info'} = "MalaCards=1=0";     
                                                               
    $hashAnnotationList{'HGMD'}{'Species'}                     = '[Human_hg19][Human_hg19_rmHAP][Human_hg38][Human_hg19_HLA]';                       
    $hashAnnotationList{'HGMD'}{'From'}                        = 'File';                       
    $hashAnnotationList{'HGMD'}{'NecessaryFile'}{1}            = "prefix=\tsuffix=/home/genesky/database/annovar/hg19/hgmd/hg19_HGMD_Gene.txt";  
    $hashAnnotationList{'HGMD'}{'GetResult'}                   = 'General Gene';                                             
    $hashAnnotationList{'HGMD'}{'Output'}{1}{'Need_Info'}      = "HGMD=1=0";     
    
    return %hashAnnotationList;    
}                                                                          

1