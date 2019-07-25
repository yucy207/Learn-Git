package package::annotation_db;
use strict;
use warnings;

# 获取注释列表
sub get_annotation {
    my %hashAnnotationList;
    # 关键字说明
    # （1）Species       : 该注释对应的物种，用"[物种名]"表示,如果对应多个物种，用","分割
    # （2）From          : 当前注释的来源方式，Annovar表示需要用annovar进行注释（后面的Para必须提供）;File表示需要写代码从文件提取
    # （3）Para          : Annovar注释时，该参数提供注释参数
    # （4）DBDir         : 当前注释文件所在路径
    ###
    # 区域注释
    ###
    $hashAnnotationList{'refGene'}{'Species'}              = '[Human_hg19]';
    $hashAnnotationList{'refGene'}{'From'}                 = 'Annovar';
    $hashAnnotationList{'refGene'}{'Para'}                 = 'annotate_variation.pl --separate --hgvs --splicing_threshold 15';
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19'}  = "/home/genesky/database/ucsc/hg19/gene/";
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg38'}  = "/home/genesky/database/ucsc/hg38/gene";

    $hashAnnotationList{'cytoBand'}{'Species'}             = '[Human_hg19]';
    $hashAnnotationList{'cytoBand'}{'From'}                = 'Annovar';
    $hashAnnotationList{'cytoBand'}{'Para'}                = 'annotate_variation.pl -regionanno -dbtype cytoBand';
    $hashAnnotationList{'cytoBand'}{'DBDir'}{'Human_hg19'} = "/home/genesky/database/self_build_database/human_annovar/hg19/cytoband";

	$hashAnnotationList{'dgvFreq'}{'Species'}              = '[Human_hg19]';
    $hashAnnotationList{'dgvFreq'}{'From'}                 = 'Annovar';
    $hashAnnotationList{'dgvFreq'}{'Para'}                 = 'annotate_variation.pl -regionanno -dbtype dgvFreq';
    $hashAnnotationList{'dgvFreq'}{'DBDir'}{'Human_hg19'}  = "/home/genesky/database/self_build_database/human_annovar/hg19/dgvfreq";
	
	$hashAnnotationList{'CNVD'}{'Species'}                 = '[Human_hg19]';
    $hashAnnotationList{'CNVD'}{'From'}                    = 'Annovar';
    $hashAnnotationList{'CNVD'}{'Para'}                    = 'annotate_variation.pl -regionanno -dbtype CNVD';
    $hashAnnotationList{'CNVD'}{'DBDir'}{'Human_hg19'}     = "/home/genesky/database/self_build_database/human_annovar/hg19/cnvd";
	
	$hashAnnotationList{'DECIPHER'}{'Species'}             = '[Human_hg19]';
    $hashAnnotationList{'DECIPHER'}{'From'}                = 'Annovar';
    $hashAnnotationList{'DECIPHER'}{'Para'}                = 'annotate_variation.pl -regionanno -dbtype DECIPHER';
    $hashAnnotationList{'DECIPHER'}{'DBDir'}{'Human_hg19'} = "/home/genesky/database/self_build_database/human_annovar/hg19/decipher";
	
	$hashAnnotationList{'iscaPathGainCum'}{'Species'}                   = '[Human_hg19]';
    $hashAnnotationList{'iscaPathGainCum'}{'From'}                      = 'Annovar';
    $hashAnnotationList{'iscaPathGainCum'}{'Para'}                      = 'annotate_variation.pl -regionanno -dbtype iscaPathGainCum';
    $hashAnnotationList{'iscaPathGainCum'}{'DBDir'}{'Human_hg19'}       = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_path_gain_cum";
	                                                                    
	$hashAnnotationList{'iscaPathLossCum'}{'Species'}                   = '[Human_hg19]';
    $hashAnnotationList{'iscaPathLossCum'}{'From'}                      = 'Annovar';
    $hashAnnotationList{'iscaPathLossCum'}{'Para'}                      = 'annotate_variation.pl -regionanno -dbtype iscaPathLossCum';
    $hashAnnotationList{'iscaPathLossCum'}{'DBDir'}{'Human_hg19'}       = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_path_loss_cum";
	                                                                    
	$hashAnnotationList{'iscaLikelyPathogenic'}{'Species'}              = '[Human_hg19]';
    $hashAnnotationList{'iscaLikelyPathogenic'}{'From'}                 = 'Annovar';
    $hashAnnotationList{'iscaLikelyPathogenic'}{'Para'}                 = 'annotate_variation.pl -regionanno -dbtype iscaLikelyPathogenic';
    $hashAnnotationList{'iscaLikelyPathogenic'}{'DBDir'}{'Human_hg19'}  = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_likely_pathogenic";
	
	$hashAnnotationList{'iscaPathogenic'}{'Species'}                    = '[Human_hg19]';
    $hashAnnotationList{'iscaPathogenic'}{'From'}                       = 'Annovar';
    $hashAnnotationList{'iscaPathogenic'}{'Para'}                       = 'annotate_variation.pl -regionanno -dbtype iscaPathogenic';
    $hashAnnotationList{'iscaPathogenic'}{'DBDir'}{'Human_hg19'}        = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_pathogenic";
	
	$hashAnnotationList{'iscaCuratedPathogenic'}{'Species'}             = '[Human_hg19]';
    $hashAnnotationList{'iscaCuratedPathogenic'}{'From'}                = 'Annovar';
    $hashAnnotationList{'iscaCuratedPathogenic'}{'Para'}                = 'annotate_variation.pl -regionanno -dbtype iscaCuratedPathogenic';
    $hashAnnotationList{'iscaCuratedPathogenic'}{'DBDir'}{'Human_hg19'} = "/home/genesky/database/self_build_database/human_annovar/hg19/isca_curated_pathogenic";
	
    return %hashAnnotationList;    
}

1