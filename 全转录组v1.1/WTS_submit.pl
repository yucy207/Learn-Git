#!/usr/bin/env perl

BEGIN{
	push @INC, qq{/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/};
}

use Getopt::Long;
use package::config;
use package::get_time;
use package::parallel;

######## 常规分析 ########
# QC
use package::data_eval;
use package::quality_control;
use package::reads_evaluate;
# mapping
use package::mapping;
use package::rseqc;
use package::RNA_evaluate;
# assemble
use package::assemble_isoform;
use package::merge_transcript;
# mRNA
use package::mRNA::mRNA_stringTie;
use package::mRNA::mRNA_deseq;
use package::mRNA::mRNA_enrichment;
use package::mRNA::mRNA_enrichment_GSEA;
use package::mRNA::mRNA_rmats;
use package::mRNA::mRNA_new_isoform;
use package::mRNA::mRNA_fusion;
use package::mRNA::mRNA_ppi;
use package::mRNA::mRNA_wgcna;
# lncRNA
use package::lncRNA::lncRNA_stringTie;
use package::lncRNA::lncRNA_deseq;
use package::lncRNA::lncRNA_target_deseq;
use package::lncRNA::lncRNA_enrichment;
use package::lncRNA::novel_lncRNA_identity;
# circRNA
use package::circRNA::circRNA_predict;
use package::circRNA::circRNA_profile;
use package::circRNA::circRNA_diff;
use package::circRNA::circRNA_enrich;
use package::circRNA::circRNA_target;
# association
use package::association::mRNA_miRNA;
use package::association::lncRNA_miRNA;
use package::association::circRNA_miRNA;
use package::association::lncRNA_cis;
use package::association::lncRNA_trans;
# ceRNA
use package::ceRNA::mRNA_miRNA_lncRNA;
use package::ceRNA::mRNA_miRNA_circRNA;

####### 个性化分析 #######
# overlap
use package::mRNA::mRNA_overlap_enrichment;
use package::lncRNA::lncRNA_overlap_target_deseq;
use package::lncRNA::lncRNA_overlap_enrichment;
use package::circRNA::circRNA_overlap_enrich;
# anova
use package::mRNA::mRNA_anova;
use package::mRNA::mRNA_anova_enrichment;
use package::lncRNA::lncRNA_anova;
use package::lncRNA::lncRNA_target_anova;
use package::lncRNA::lncRNA_anova_enrichment;
use package::circRNA::circRNA_anova;
use package::circRNA::circRNA_anova_enrich;
# time
use package::mRNA::mRNA_maSigPro;
use package::mRNA::mRNA_maSigPro_enrichment;
use package::lncRNA::lncRNA_maSigPro;
use package::lncRNA::lncRNA_maSigPro_target;
use package::lncRNA::lncRNA_maSigPro_enrichment;
use package::circRNA::circRNA_maSigPro;
use package::circRNA::circRNA_maSigPro_enrich;

my ($skip, 
	$only, 
	$force_sample, 
	$force_step,   
	$thread_qc, 
	$thread_map, 
	$adapter,
	$log2fc,
	$log2fc_circ,
	$help, 
	$list);

GetOptions(
	'skip=s'          => \$skip,
	'only=s'          => \$only,
	'force_sample=s'  => \$force_sample,
	'force_step=s'    => \$force_step,
	'thread_qc=i'     => \$thread_qc,
	'thread_map=i'    => \$thread_map,
	'adapter=s'       => \$adapter,
	'log2fc=f'        => \$log2fc,
	'log2fc_circ=f'   => \$log2fc_circ,
 	'help|h!'         => \$help,
 	'list!'           => \$list
	
);

if ($list) {
my $help =<<EOF;
1 质量评估
2 质量控制
3 数据量统计
4 比对参考基因组
5 文库QC
6 转录本组装
7 mRNA分析
8 lncRNA分析
9 circRNA分析
10 联合分析
EOF
print qq{$help\n};
exit;	
}

if ($help or $#ARGV == -1) {
	usage();
	exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl $0  config.txt
    [ Options ]
    --------------- steps  ----------------
    -skip          skip steps, eg [3,4,5]
    -only          only steps, eg [1,2,3]

    --------------- force -----------------
    -force_sample  force run sample, eg [sampleA, sampleB]          
    -force_step    force run steps,  eg [1,2,3]

    -------------- threads ----------------
    -thread_qc     thread of qc  [option]
    -thread_map    thread of map [option]

    -------------- others -----------------
    -adapter       adapter sequence [option:truseq, nextera]
    -log2fc        log2fold_change in mRNA/lncRNA different analysis, eg [1]
    -log2fc_circ   log2fold_change in circRNA different analysis, eg [0]

    -help|-h       print help message
    -list          list steps 

    [ Arguments ]
    config.txt
EOF
	print $help;
	exit;

}

if ($only and $skip) {
	die "you can just pick only one of options [only / skip]\n";
}

my ($config) = $ARGV[$#ARGV];

my $metadata   = package::config::read_config($config);
my $base       = package::config::base_config();
my @samples    = split/\,/,$metadata->{'samples'};
my $sample_num = @samples;

##### parse command options ########################
$base->{'thread_qc'}    = $thread_qc    if $thread_qc;
$base->{'thread_map'}   = $thread_map   if $thread_map;

$base->{'force_sample'} = $force_sample if $force_sample;
$base->{'force_step'}   = $force_step   if $force_step;

$base->{'log2fc'}       = $log2fc       if $log2fc;
$base->{'log2fc_circ'}  = $log2fc_circ  if $log2fc_circ;

my @run_steps = ();
if ($only) {
	foreach my $x (split /,/, $only) {
		push @run_steps, $x;
	}
} else {
	my @skips = split /,/, $skip;
	foreach my $x (1..10) {
		next if $x ~~ @skips;
		push @run_steps, $x;
	}
}

$SIG{INT} = sub {
    print "you has enter the CTRL+C keys,now exit\n";
    exit 0;
};

my $suffix = "="  x 30;

if (1 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行质量评估${suffix}\n};
	package::data_eval::run($metadata, $base);
}


if (2 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行质量控制${suffix}\n};
	package::quality_control::run($metadata, $base, $adapter);
}


if (3 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行数据量统计${suffix}\n};
	package::reads_evaluate::run($metadata, $base);
}

if (4 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行比对参考基因组${suffix}\n};
	package::mapping::run($metadata, $base);
}


if (5 ~~ @run_steps) {
	
	print package::get_time::run().qq{${suffix}运行文库QC${suffix}\n};
	my $rseqc_proc = sub {package::rseqc::run($metadata, $base);};
	my $rna_proc   = sub {package::RNA_evaluate::run($metadata, $base);};
	package::parallel::run_sub($rseqc_proc, $rna_proc);

}

if (6 ~~ @run_steps) {

	print package::get_time::run().qq{${suffix}运行转录本组装${suffix}\n};
	package::assemble_isoform::run($metadata, $base);
	print package::get_time::run().qq{${suffix}运行合并转录本${suffix}\n};
	package::merge_transcript::run($metadata, $base);

}

if (7 ~~ @run_steps) {
	
	print package::get_time::run().qq{${suffix}运行mRNA 定量分析${suffix}\n};
	package::mRNA::mRNA_stringTie::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行mRNA 差异分析${suffix}\n};
	package::mRNA::mRNA_deseq::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行mRNA 富集分析${suffix}\n};
	package::mRNA::mRNA_enrichment::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行mRNA GSEA富集分析${suffix}\n};
	package::mRNA::mRNA_enrichment_GSEA::run($metadata, $base);

	if (exists $metadata->{'overlap'}){

		print package::get_time::run().qq{${suffix}运行mRNA 差异的overlap基因富集分析${suffix}\n};
		package::mRNA::mRNA_overlap_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置overlap的组别,因此跳过overlap基因富集分析!\n";
	}

	if (exists $metadata->{'anova'}){

		print package::get_time::run().qq{${suffix}运行mRNA anova多组差异分析${suffix}\n};
		package::mRNA::mRNA_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行mRNA anova后续富集分析${suffix}\n};
		package::mRNA::mRNA_anova_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置anova分析的组别,因此跳过anova相关分析!\n";
	}

	if (exists $metadata->{'time'}){

		print package::get_time::run().qq{${suffix}运行mRNA 时间序列分析${suffix}\n};
		package::mRNA::mRNA_maSigPro::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行mRNA 时间序列后续富集分析${suffix}\n};
		package::mRNA::mRNA_maSigPro_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置design.txt的路径,因此跳过时间序列相关分析!\n";
	}

	print package::get_time::run().qq{${suffix}运行mRNA 可变剪切分析${suffix}\n};
	package::mRNA::mRNA_rmats::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行mRNA 转录本分析${suffix}\n};
	package::mRNA::mRNA_new_isoform::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行mRNA 融合基因分析${suffix}\n};
	package::mRNA::mRNA_fusion::run($metadata, $base);
	
	print package::get_time::run().qq{${suffix}运行mRNA ppi分析&&富集分析${suffix}\n};
	package::mRNA::mRNA_ppi::run($metadata, $base);

	if ($sample_num > 15){

		print package::get_time::run().qq{${suffix}运行mRNA WGCNA分析&&富集分析${suffix}\n};
		package::mRNA::mRNA_wgcna::run($metadata, $base);
		
	}else{
		print "样本个数小于15个,无法进行WGCNA分析!\n";
	}

}

if (8 ~~ @run_steps) {

	print package::get_time::run().qq{${suffix}运行lncRNA表达定量${suffix}\n};
	package::lncRNA::lncRNA_stringTie::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行lncRNA差异分析${suffix}\n};
	package::lncRNA::lncRNA_deseq::run($metadata, $base);

	print package::get_time::run().qq{${suffix}lncRNA靶基因预测分析${suffix}\n};
	package::lncRNA::lncRNA_target_deseq::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行lncRNA富集分析${suffix}\n};
	package::lncRNA::lncRNA_enrichment::run($metadata, $base);

	print package::get_time::run().qq{${suffix}lncRNA新转录本分析${suffix}\n};
	package::lncRNA::novel_lncRNA_identity::run($metadata, $base);

	if (exists $metadata->{'overlap'}){

		print package::get_time::run().qq{${suffix}运行lncRNA overlap靶基因预测分析${suffix}\n};
		package::lncRNA::lncRNA_overlap_target_deseq::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行lncRNA overlap富集分析${suffix}\n};
		package::lncRNA::lncRNA_overlap_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置overlap的组别,因此跳过lncRNA overlap相关分析!\n";
	}

	if (exists $metadata->{'anova'}){

		print package::get_time::run().qq{${suffix}运行lncRNA anova差异分析${suffix}\n};
		package::lncRNA::lncRNA_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行anova差异lncRNA靶基因预测分析${suffix}\n};
		package::lncRNA::lncRNA_target_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行lncRNA anova后续富集分析${suffix}\n};
		package::lncRNA::lncRNA_anova_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置anova的组别,因此跳过lncRNA anova相关分析!\n";
	}

	if (exists $metadata->{'time'}){

		print package::get_time::run().qq{${suffix}运行lncRNA 时间序列分析${suffix}\n};
		package::lncRNA::lncRNA_maSigPro::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行lncRNA 时间序列靶基因预测${suffix}\n};
		package::lncRNA::lncRNA_maSigPro_target::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行lncRNA 时间序列靶基因富集分析${suffix}\n};
		package::lncRNA::lncRNA_maSigPro_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置design.txt的路径,因此跳过lncRNA 时间序列相关分析!\n";
	}

}

if (9 ~~ @run_steps) {

	print package::get_time::run().qq{${suffix}运行circRNA预测分析${suffix}\n};
	package::circRNA::circRNA_predict::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行circRNA定量分析${suffix}\n};
	package::circRNA::circRNA_profile::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行circRNA差异分析${suffix}\n};
	package::circRNA::circRNA_diff::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行circRNA来源基因富集分析${suffix}\n};
	package::circRNA::circRNA_enrich::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行circRNA靶标预测分析${suffix}\n};
	package::circRNA::circRNA_target::run($metadata, $base);

	if (exists $metadata->{'overlap'}){

		print package::get_time::run().qq{${suffix}运行circRNA overlap基因富集分析${suffix}\n};
		package::circRNA::circRNA_overlap_enrich::run($metadata, $base);

	}else{
		print "config里没有设置overlap的组别,因此跳过circRNA overlap基因富集分析!\n";
	}

	if (exists $metadata->{'anova'}){

		print package::get_time::run().qq{${suffix}运行circRNA anova多组差异分析${suffix}\n};
		package::circRNA::circRNA_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行circRNA anova后续富集分析${suffix}\n};
		package::circRNA::circRNA_anova_enrich::run($metadata, $base);

	}else{
		print "config里没有设置anova的组别,因此跳过circRNA anova相关分析!\n";
	}

	if (exists $metadata->{'time'}){

		print package::get_time::run().qq{${suffix}运行circRNA 时间序列分析${suffix}\n};
		package::circRNA::circRNA_maSigPro::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行circRNA 时间序列富集分析${suffix}\n};
		package::circRNA::circRNA_maSigPro_enrich::run($metadata, $base);

	}else{
		print "config里没有设置design.txt的路径,因此跳过circRNA 时间序列相关分析!\n";
	}

}

if (10 ~~ @run_steps) {

	print package::get_time::run().qq{${suffix}运行mRNA和miRNA相互作用分析${suffix}\n};
	package::association::mRNA_miRNA::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行lncRNA和miRNA相互作用分析${suffix}\n};
	package::association::lncRNA_miRNA::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行circRNA和miRNA相互作用分析${suffix}\n};
	package::association::circRNA_miRNA::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行lncRNA Cis分析${suffix}\n};
	package::association::lncRNA_cis::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行lncRNA Trans分析${suffix}\n};
	package::association::lncRNA_trans::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行mRNA_miRNA_lncRNA ceRNA分析${suffix}\n};
	package::ceRNA::mRNA_miRNA_lncRNA::run($metadata, $base);

	print package::get_time::run().qq{${suffix}运行mRNA_miRNA_circRNA ceRNA分析${suffix}\n};
	package::ceRNA::mRNA_miRNA_circRNA::run($metadata, $base);

}



