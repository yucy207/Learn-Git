#!/usr/bin/perl

BEGIN{
	push @INC,"/home/genesky/pipeline/absolute_quantification_16s/v2.1.1/packages";
}

$|=1;
use strict;
use warnings;
use Getopt::Long;
use config;
use quality_control;
use otu_noresample;
use quantitation;
use extract_and_resample;
use community;
use alphaDiversity;
use betaDiversity;
use combine;
use picrust;
use remove_rawdataIR;
use positive_control;
use environmentFactors;
use html_report;


my ($skip,
	$only,
	$help,
	$quan);

GetOptions(
	'skip=s'          => \$skip,
	'only=s'          => \$only,
	'help|h!'         => \$help,
    'quan|q=s'        => \$quan
);

Usage() if ( $help or (scalar @ARGV != 1) );

sub Usage{
my $help ="
Usage: perl $0  config.txt

    --------- [ Options ] ------------
    -skip          skip steps, eg [4,5,6,7,8]
    -only          only steps, eg [1,2,3]
    -help|-h       print help message
    -quan|-q       quantitative type of data, eg[A]

    --------- [ Arguments ] ----------
    config.txt

    --------- [ steps list ] ---------
    1 质量控制
    2 OTU聚类注释
    3 绝对拷贝数计算
    4 提取样本和抽平分析
    5 Community群落分析
    6 Alpha多样性分析
    7 Beta多样性分析
    8 功能预测分析
    9 环境因子分析
    10 联合分析
    11 html报告生成
    12 PC阳性对照样本分析
    13 删除原始fastq.gz中的spike-in
    14 删除中间结果

    --------- [quan list] -------------
    A|a 只分析绝对定量(Absolute Quantitation)
    R|r 只分析相对定量(Relative Quantitation)
    B|b 两种定量数据都分析(both AQ and RQ)，默认
\n";
print $help;
exit;
}

if ($only and $skip) {
	die "you can just pick only one of options [only / skip]\n";
}

my @run_steps = ();
if ($only) {
	foreach my $x (split /,/, $only) {
		push @run_steps, $x;
	}
} else {
	my @skips = split /,/, $skip;
	foreach my $x (1..14) {
		next if $x ~~ @skips;
		push @run_steps, $x;
	}
}

$SIG{INT} = sub {
    print "you has enter the CTRL+C keys,now exit\n";
    exit 0;
};

####分析数据的定量类型
my $quantitation = (defined $quan)? uc($quan) : "B";


my ($config) = @ARGV;
die "perl $0 config.txt \n" if scalar @ARGV != 1;

my $metadata = config::read_config($config);
my $base     = config::base_config();
			 
my $suffix = "="  x 15;
sub get_time{
	my $time = `date "+%Y-%m-%d %T"`;
	chomp $time;
	return "[$time]";
}



if (exists $metadata->{'raw_data'} and (1 ~~ @run_steps)){
	print get_time().qq{$suffix运行质量控制$suffix\n};
	quality_control::run($metadata, $base);
}

if (2 ~~ @run_steps){
	print get_time().qq{$suffix运行OTU聚类注释分析$suffix\n};
	otu_noresample::run($metadata, $base);
}


if (exists $metadata->{'group'}){

	if (3 ~~ @run_steps){
		print get_time().qq{$suffix运行绝对拷贝数计算$suffix\n};
		quantitation::run($metadata, $base);		
	}

	if (4 ~~ @run_steps){
		print get_time().qq{$suffix运行Extract_and_Resample样本提取和抽平分析$suffix\n};
		extract_and_resample::run($metadata, $base, $quantitation);		
	}

	if (5 ~~ @run_steps){
		print get_time().qq{$suffix运行Community群落组成分析$suffix\n};
		community::run($metadata, $base, $quantitation);		
	}

	if (6 ~~ @run_steps) {
		print get_time().qq{$suffix运行Alpha多样性分析$suffix\n};
		alphaDiversity::run($metadata, $base, $quantitation);	
	}

	if (7 ~~ @run_steps) {
		print get_time().qq{$suffix运行Beta多样性分析$suffix\n};
		betaDiversity::run($metadata, $base, $quantitation);		
	}

	if (8 ~~ @run_steps) {
		print get_time().qq{$suffix运行PICRUSt分析$suffix\n};
		picrust::run($metadata, $base, $quantitation);		
	}

	if (exists $metadata->{'env'} and 9 ~~ @run_steps){
		print get_time().qq{$suffix运行ENV分析$suffix\n};
		environmentFactors::run($metadata, $base, $quantitation);
	}

	if (10 ~~ @run_steps) {
		print get_time().qq{$suffix运行Combine联合分析$suffix\n};
		combine::run($metadata, $base);		
	}

}

if (11 ~~ @run_steps) {
	print get_time().qq{$suffix运行html报告制作分析$suffix\n};
	html_report::run($metadata, $base, $config);	
}

if (12 ~~ @run_steps) {
	print get_time().qq{$suffix运行PC样本分析$suffix\n};
	positive_control::run($metadata, $base);
}

if (13 ~~ @run_steps) {
	print get_time().qq{$suffix删除原始数据中spike-in序列$suffix\n};
	remove_rawdataIR::run($metadata, $base);	
}


if (14 ~~ @run_steps) {
	print get_time().qq{$suffix所有分析已完成$suffix\n是否删除中间结果？[y/n]};
	my $option = <STDIN>;
	$option =~ s/[\r\n]//g;
	if ($option eq "y"){
		
		system qq{cat $metadata->{intm_result}/qc/*/*.clean.fasta > $metadata->{result}/after.QC.reads.clean.fasta};
		system qq{gzip $metadata->{result}/after.QC.reads.clean.fasta};
		system qq{rm -rf $metadata->{intm_result}};

		print qq{中间结果已删除!\n};
	}	
}

