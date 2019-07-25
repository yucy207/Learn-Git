BEGIN{
	push @INC,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/package/mRNA_SNP/tools/report/";
}

$|=1;
use strict;
use warnings;
use Excel::Writer::XLSX;

use package::config;
use package::main;
use package::time;
use package::format;
use package::filter;
use package::priority;
use package::gene;
use package::sample;
use package::snv;
use package::genotype;
use package::readme;
use package::gender;
use package::p1p2p3p4;
my $tab="##########";

# 读取配置文件
die "Usage perl $0 ConfigFile\n" if(@ARGV!=1);
my $configFile = shift @ARGV;
my %config     = package::config::read_config($configFile);
my @sample     = split /,/, $config{'samples'};
my $result     = "$config{'project'}/mRNA/snp/result";

# 生成run.ini
my $runConfig = "$result/run.ini";
open(RUN,">$runConfig") or die "cannot make run.ini!\n";
print RUN "::RUN LowFreq_Func_Mutations\n";
if (exists $config{'group'}) {
	my @groups = @{$config{'group'}};
	my $control;
	my $case;
	foreach my $x(@groups){
		$control= (split /;/, $x->[2])[0];
    	$case   = (split /;/, $x->[2])[1];
	}
	print RUN "Case:$case\n";
	print RUN "Control:$control\n";
}else {
	print RUN "Case:$config{'samples'}\n";
	print RUN "Control:\n";
}
print RUN "Freq_Alt (1000g):0.01\n";
print RUN "ExAC03_EAS:0.01\n";
print RUN "gnomAD_exome_EAS:0.01\n";
print RUN "GeneskyExonDB_Freq:0.05\n";
print RUN "Model:Dominance\n";
print RUN "GeneControlModel:normal\n";
print RUN "Rm_Function:unknown\n";
print RUN "Report:All\n";
print RUN "MOD_SET:ONE\n";
print RUN "OR_Site:";
foreach my $x(@sample){
	print RUN "$x,1,";
}
print RUN "\n";
print RUN "OR_Site:";
foreach my $x(@sample){
	print RUN "$x,2,";
}
close RUN;

my %run = package::main::readRunConfig($runConfig);

# 开始分析
die "[ERR] Loss Library\n" if(not -e "$result/database/database.snv");
die "[ERR] Loss Annotation\n" if(not -e "$result/annotation.txt");
die "[ERR] Loss Analysis\n" if(not -e "$result/database/analysis.txt");
die "[ERR] Loss Gender\n" if(not -e "$result/gender/gender.log");
die "[ERR] Loss ReadMe File\n" if(not -e "/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/package/mRNA_SNP/tools/report/readme.txt");

my %library    = package::main::readLibrary("$result/database/database.snv");
my %annotation = package::main::readAnnotation("$result/annotation.txt");
my %analysis   = package::main::readAnalysis("$result/database/analysis.txt");
my %gender     = package::gender::read("$result/gender/gender.log");
my %hashassigngender = package::gender::read2(\%config);

my $report = "$config{'report'}/03_mRNA_Analysis/09_mRNA_SNP_Analysis";
mkdir $report if(not -e $report);

foreach my $Excel(sort keys %run){
	# Excel
	my $workbook = Excel::Writer::XLSX->new("$report/$Excel.xlsx");
	my %format   = package::format::run($workbook);
	my %filter   = package::filter::run(\%library,\%annotation,\%analysis,$run{$Excel});
	# Sheet
	my $GENE     = $workbook->add_worksheet("Gene");
	my $SAMPLE   = $workbook->add_worksheet("Sample");
	my $SNV      = $workbook->add_worksheet("SNV Information");
	my $Genotype = $workbook->add_worksheet("GENOTYPES");
	my $Call     = $workbook->add_worksheet("GENOTYPE CALLING");
	my $Readme   = $workbook->add_worksheet("READ ME");

	package::gene::run($GENE,$report,$Excel,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$Excel});
	package::sample::run($SAMPLE,\%format,\%library,\%annotation,\%analysis,\%filter,\%gender,$run{$Excel},\%hashassigngender);
	package::snv::run($SNV,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$Excel});
	package::genotype::run($Genotype,$Call,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$Excel});
	package::readme::run($Readme,\%format,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/package/mRNA_SNP/tools/report/readme.txt");
}

print "低频功能突变输出完成!\n";