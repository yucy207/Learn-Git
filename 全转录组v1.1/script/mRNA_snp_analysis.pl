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


#SNP
use package::mRNA_SNP::utils;
use package::mRNA_SNP::annotation_db;
use package::mRNA_SNP::gatk;
use package::mRNA_SNP::gvcf;
use package::mRNA_SNP::call;
use package::mRNA_SNP::library;
use package::mRNA_SNP::gender;
use package::mRNA_SNP::database;
use package::mRNA_SNP::annotation;

my ($skip, 
	$only,      
	$help, 
	$list);

GetOptions(
	'skip=s'          => \$skip,
	'only=s'          => \$only,
 	'help|h!'         => \$help,
 	'list!'           => \$list
);

if ($list) {
my $help =<<EOF;
1 生成final bam
2 gvcf分析
3 call分析
4 library分析
5 gender分析
6 database分析
7 annotation分析
8 allSNV统计
9 低频功能突变输出
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

my @run_steps = ();
if ($only) {
	foreach my $x (split /,/, $only) {
		push @run_steps, $x;
	}
} else {
	my @skips = split /,/, $skip;
	foreach my $x (1..9) {
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
	print package::get_time::run().qq{${suffix}运行突变分析：生成final bam${suffix}\n};
	package::mRNA_SNP::gatk::run($metadata, $base);
}

if (2 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：gvcf分析${suffix}\n};
	package::mRNA_SNP::gvcf::run($metadata, $base);
}

if (3 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：call分析${suffix}\n};
	package::mRNA_SNP::call::run($metadata, $base);
}

if (4 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：library分析${suffix}\n};
	package::mRNA_SNP::library::run($metadata, $base);
}

if (5 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：gender分析${suffix}\n};
	package::mRNA_SNP::gender::run($metadata, $base);
}

if (6 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：database分析${suffix}\n};
	package::mRNA_SNP::database::run($metadata, $base);
}

if (7 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：annotation分析${suffix}\n};
	package::mRNA_SNP::annotation::run($metadata, $base);
}

if (8 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：allSNV统计${suffix}\n};
	system qq{perl /home/genesky/pipeline/whole_transriptome_sequencing/v1.1/package/mRNA_SNP/tools/allsnv/snvall.pl $config allSNV};
}

if (9 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行突变分析：低频功能突变输出${suffix}\n};
	system qq{perl /home/genesky/pipeline/whole_transriptome_sequencing/v1.1/package/mRNA_SNP/tools/report/report.pl $config};
}

	


	














