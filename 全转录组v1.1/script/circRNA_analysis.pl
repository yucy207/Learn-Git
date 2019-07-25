#!/usr/bin/env perl


BEGIN{
	push @INC,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/";
}

use Getopt::Long;
use package::config;
use package::get_time;
use package::parallel;

# 常规
use package::circRNA::circRNA_predict;
use package::circRNA::circRNA_profile;
use package::circRNA::circRNA_diff;
use package::circRNA::circRNA_enrich;
use package::circRNA::circRNA_target;
# 个性
use package::circRNA::circRNA_overlap_enrich;
use package::circRNA::circRNA_anova;
use package::circRNA::circRNA_anova_enrich;
use package::circRNA::circRNA_maSigPro;
use package::circRNA::circRNA_maSigPro_enrich;

my ($skip, 
	$only, 
	$log2fc_circ,     
	$help, 
	$list);

GetOptions(
	'skip=s'          => \$skip,
	'only=s'          => \$only,
	'log2fc_circ=f'   => \$log2fc_circ,
 	'help|h!'         => \$help,
 	'list!'           => \$list
);

if ($list) {
my $help =<<EOF;
1 circRNA预测分析
2 circRNA定量分析
3 circRNA差异分析
4 circRNA来源基因富集分析
5 circRNA靶标预测分析
6 circRNA overlap基因富集分析
7 circRNA anova多组差异分析及富集分析
8 circRNA 时间序列分析及富集分析
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
$base->{'log2fc_circ'}  = $log2fc_circ  if $log2fc_circ;

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
	print package::get_time::run().qq{${suffix}运行circRNA预测分析${suffix}\n};
	package::circRNA::circRNA_predict::run($metadata, $base);
}

if (2 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行circRNA定量分析${suffix}\n};
	package::circRNA::circRNA_profile::run($metadata, $base);
}

if (3 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行circRNA差异分析${suffix}\n};
	package::circRNA::circRNA_diff::run($metadata, $base);
}

if (4 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行circRNA来源基因富集分析${suffix}\n};
	package::circRNA::circRNA_enrich::run($metadata, $base);
}

if (5 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行circRNA靶标预测分析${suffix}\n};
	package::circRNA::circRNA_target::run($metadata, $base);
}
if (6 ~~ @run_steps) {
	if (exists $metadata->{'overlap'}){

		print package::get_time::run().qq{${suffix}运行circRNA overlap基因富集分析${suffix}\n};
		package::circRNA::circRNA_overlap_enrich::run($metadata, $base);

	}else{
		print "config里没有设置overlap的组别,因此跳过circRNA overlap基因富集分析!\n";
	}
}

if (7 ~~ @run_steps) {
	if (exists $metadata->{'anova'}){

		print package::get_time::run().qq{${suffix}运行circRNA anova差异分析${suffix}\n};
		package::circRNA::circRNA_anova::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行circRNA anova后续富集分析${suffix}\n};
		package::circRNA::circRNA_anova_enrich::run($metadata, $base);

	}else{
		print "config里没有设置anova的组别,因此跳过circRNA anova相关分析!\n";
	}
}

if (8 ~~ @run_steps) {
	if (exists $metadata->{'time'}){

		print package::get_time::run().qq{${suffix}运行circRNA 时间序列分析${suffix}\n};
		package::circRNA::circRNA_maSigPro::run($metadata, $base);

		print package::get_time::run().qq{${suffix}运行circRNA 时间序列富集分析${suffix}\n};
		package::circRNA::circRNA_maSigPro_enrich::run($metadata, $base);

	}else{
		print "config里没有设置design.txt的路径,因此跳过circRNA 时间序列相关分析!\n";
	}
}







	


	














