#!/usr/bin/env perl


BEGIN{
	push @INC,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/";
}

use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Statistics::R;
use package::config;
use package::get_time;
use package::parallel;

# 常规
use package::mRNA::mRNA_stringTie;
use package::mRNA::mRNA_deseq;
use package::mRNA::mRNA_enrichment;
use package::mRNA::mRNA_enrichment_GSEA;
use package::mRNA::mRNA_rmats;
use package::mRNA::mRNA_new_isoform;
use package::mRNA::mRNA_fusion;
use package::mRNA::mRNA_ppi;
use package::mRNA::mRNA_wgcna;
# 个性
use package::mRNA::mRNA_overlap_enrichment;
use package::mRNA::mRNA_anova;
use package::mRNA::mRNA_anova_enrichment;
use package::mRNA::mRNA_maSigPro;
use package::mRNA::mRNA_maSigPro_enrichment;

my ($skip, 
	$only,
	$log2fc,      
	$help, 
	$list);

GetOptions(
	'skip=s'          => \$skip,
	'only=s'          => \$only,
	'log2fc=f'        => \$log2fc,
 	'help|h!'         => \$help,
 	'list!'           => \$list
);

if ($list) {
my $help =<<EOF;
1 stringtie定量分析
2 deseq差异分析
3 富集分析
4 GSEA富集分析
5 可变剪切分析
6 转录本分析
7 融合基因分析
8 ppi分析和富集分析
9 WGCNA分析和富集分析
10 差异overlap基因富集分析
11 anova多组差异分析及后续富集分析
12 时间序列分析及后续富集分析
EOF
print qq{$help\n};
exit;	
}

if ($help) {
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
    -log2fc        log2fold_change in mRNA/lncRNA different analysis, eg [1]
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

my $metadata         = package::config::read_config($config);
my $base             = package::config::base_config();
   $base->{'log2fc'} = $log2fc if $log2fc;
my @samples          = split/\,/,$metadata->{'samples'};
my $sample_num       = @samples;



my @run_steps = ();
if ($only) {
	foreach my $x (split /,/, $only) {
		push @run_steps, $x;
	}
} else {
	my @skips = split /,/, $skip;
	foreach my $x (1..12) {
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
	print package::get_time::run().qq{${suffix}运行mRNA 定量分析${suffix}\n};
	package::mRNA::mRNA_stringTie::run($metadata, $base);
}

if (2 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA 差异分析${suffix}\n};
	package::mRNA::mRNA_deseq::run($metadata, $base);
}

if (3 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA 富集分析${suffix}\n};
	package::mRNA::mRNA_enrichment::run($metadata, $base);
}

if (4 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA GSEA富集分析${suffix}\n};
	package::mRNA::mRNA_enrichment_GSEA::run($metadata, $base);
}

if (5 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA 可变剪切分析${suffix}\n};
	package::mRNA::mRNA_rmats::run($metadata, $base);
}

if (6 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA 转录本分析${suffix}\n};
	package::mRNA::mRNA_new_isoform::run($metadata, $base);
}

if (7 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA 融合基因分析${suffix}\n};
	package::mRNA::mRNA_fusion::run($metadata, $base);
}

if (8 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA ppi分析&&富集分析${suffix}\n};
	package::mRNA::mRNA_ppi::run($metadata, $base);
}

if (9 ~~ @run_steps) {
	if ($sample_num >= 15){

		print package::get_time::run().qq{${suffix}运行mRNA WGCNA分析&&富集分析${suffix}\n};
		package::mRNA::mRNA_wgcna::run($metadata, $base);
		
	}else{
		print "样本个数小于15个,无法进行WGCNA分析!\n";
	}
}

if (10 ~~ @run_steps) {
	if (exists $metadata->{'overlap'}){

		print package::get_time::run().qq{${suffix}运行mRNA 差异的overlap基因富集分析${suffix}\n};
		package::mRNA::mRNA_overlap_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置overlap的组别,因此跳过overlap基因富集分析!\n";
	}
}

if (11 ~~ @run_steps) {
	if (exists $metadata->{'anova'}){

		print package::get_time::run().qq{${suffix}运行mRNA anova多组差异分析${suffix}\n};
		package::mRNA::mRNA_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行mRNA anova后续富集分析${suffix}\n};
		package::mRNA::mRNA_anova_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置anova分析的组别,因此跳过anova相关分析!\n";
	}
}

if (12 ~~ @run_steps) {
	if (exists $metadata->{'time'}){

		print package::get_time::run().qq{${suffix}运行mRNA 时间序列分析${suffix}\n};
		package::mRNA::mRNA_maSigPro::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行mRNA 时间序列后续富集分析${suffix}\n};
		package::mRNA::mRNA_maSigPro_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置design.txt的路径,因此跳过时间序列相关分析!\n";
	}
}








	


	














