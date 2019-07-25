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
    $hashAnnotationList{'refGene'}{'Species'}              = '[Human_hg19][Human_hg38]';
    $hashAnnotationList{'refGene'}{'From'}                 = 'Annovar';
    $hashAnnotationList{'refGene'}{'Para'}                 = 'annotate_variation.pl --hgvs --splicing_threshold 8';
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg19'}  = "/home/genesky/database/ucsc/hg19/gene/";
    $hashAnnotationList{'refGene'}{'DBDir'}{'Human_hg38'}  = "/home/genesky/database/ucsc/hg38/gene";
	
    return %hashAnnotationList;    
}

1