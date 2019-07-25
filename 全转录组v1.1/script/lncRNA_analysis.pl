#!/usr/bin/env perl


BEGIN{
	push @INC,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/";
}

use Getopt::Long;
use package::config;
use package::get_time;
use package::parallel;

# 常规
use package::lncRNA::lncRNA_stringTie;
use package::lncRNA::lncRNA_deseq;
use package::lncRNA::lncRNA_target_deseq;
use package::lncRNA::lncRNA_enrichment;
use package::lncRNA::novel_lncRNA_identity;
# 个性
use package::lncRNA::lncRNA_overlap_target_deseq;
use package::lncRNA::lncRNA_overlap_enrichment;
use package::lncRNA::lncRNA_anova;
use package::lncRNA::lncRNA_target_anova;
use package::lncRNA::lncRNA_anova_enrichment;
use package::lncRNA::lncRNA_maSigPro;
use package::lncRNA::lncRNA_maSigPro_target;
use package::lncRNA::lncRNA_maSigPro_enrichment;

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
3 靶基因预测分析
4 富集分析
5 新转录本分析
6 lncRNA overlap相关分析
7 lncRNA anova相关分析
8 lncRNA 时间序列相关分析
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

my $metadata = package::config::read_config($config);
my $base     = package::config::base_config();
$base->{'log2fc'} = $log2fc if $log2fc;

my @run_steps = ();
if ($only) {
	foreach my $x (split /,/, $only) {
		push @run_steps, $x;
	}
} else {
	my @skips = split /,/, $skip;
	foreach my $x (1..8) {
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
	print package::get_time::run().qq{${suffix}运行lncRNA表达定量${suffix}\n};
	package::lncRNA::lncRNA_stringTie::run($metadata, $base);
}

if (2 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行lncRNA差异分析${suffix}\n};
	package::lncRNA::lncRNA_deseq::run($metadata, $base);
}

if (3 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}lncRNA靶基因预测分析${suffix}\n};
	package::lncRNA::lncRNA_target_deseq::run($metadata, $base);
}

if (4 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行lncRNA富集分析${suffix}\n};
	package::lncRNA::lncRNA_enrichment::run($metadata, $base);
}

if (5 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}lncRNA新转录本分析${suffix}\n};
	package::lncRNA::novel_lncRNA_identity::run($metadata, $base);
}

if (6 ~~ @run_steps) {
	if (exists $metadata->{'overlap'}){

		print package::get_time::run().qq{${suffix}运行lncRNA overlap靶基因预测分析${suffix}\n};
		package::lncRNA::lncRNA_overlap_target_deseq::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行lncRNA overlap富集分析${suffix}\n};
		package::lncRNA::lncRNA_overlap_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置overlap的组别,因此跳过lncRNA overlap相关分析!\n";
	}
}

if (7 ~~ @run_steps) {
	if (exists $metadata->{'anova'}){

		print package::get_time::run().qq{${suffix}运行lncRNA anova多组差异分析${suffix}\n};
		package::lncRNA::lncRNA_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行anova差异lncRNA靶基因预测分析${suffix}\n};
		package::lncRNA::lncRNA_target_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行lncRNA anova后续富集分析${suffix}\n};
		package::lncRNA::lncRNA_anova_enrichment::run($metadata, $base);

	}else{
		print "config里没有设置anova的组别,因此跳过lncRNA anova相关分析!\n";
	}
}

if (8 ~~ @run_steps) {
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




	


	














