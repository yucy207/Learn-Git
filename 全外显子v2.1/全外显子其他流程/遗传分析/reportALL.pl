$|=1;
use strict;
use warnings;
use Excel::Writer::XLSX;
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

# 读取配置
die "Usage perl $0 ConfigFile RunConfig\n" if(@ARGV!=2);
my $configFile=shift @ARGV;
my $runConfig=shift @ARGV;
my %config=package::main::readConfig($configFile);
my %run=package::main::readRunConfig($runConfig);

# 开始流程
system "clear";
print "$tab Annotation ".package::time::ymd." ".package::time::hms." $tab\n";
die "[ERR] Loss Report DIR\n" if(not exists $config{"Report"} or not -e $config{"Report"});
die "[ERR] Loss Library\n" if(not -e $config{"Report"}."/database.snv");
die "[ERR] Loss Annotation\n" if(not -e $config{"Report"}."/annotation.txt");
die "[ERR] Loss Analysis\n" if(not -e $config{"Report"}."/analysis_all.txt");
die "[ERR] Loss Gender\n" if(not -e $config{"Report"}."/gender.log");
die "[ERR] Loss ReadMe File\n" if(not -e "readme.txt");
my %library=package::main::readLibrary($config{"Report"}."/database.snv");
my %annotation=package::main::readAnnotation($config{"Report"}."/annotation.txt");
my %analysis=package::main::readAnalysis($config{"Report"}."/analysis_all.txt");
my %gender=package::gender::read($config{"Report"}."/gender.log");
my %hashassigngender=package::gender::read2(\%config);#自定义性别或完全可以通过tag判定性别
mkdir $config{"Report"}."/document" if(not -e $config{"Report"}."/document");
mkdir $config{"Report"}."/document/2_SNV" if(not -e $config{"Report"}."/document/2_SNV");
foreach my $Excel(keys %run){
	# Excel
	print $config{"Report"}."/document/2_SNV/$Excel.xlsx\n";
	my $workbook=Excel::Writer::XLSX->new($config{"Report"}."/document/2_SNV/$Excel.xlsx");
	my %format=package::format::run($workbook);
	my %filter=package::filter::run(\%library,\%annotation,\%analysis,$run{$Excel});
	# Sheet
	my $GENE=$workbook->add_worksheet("Gene");
	my $SAMPLE=$workbook->add_worksheet("Sample");
	my $SNV=$workbook->add_worksheet("SNV Information");
	my $Genotype=$workbook->add_worksheet("GENOTYPES");
	my $Call=$workbook->add_worksheet("GENOTYPE CALLING");
	my $Readme=$workbook->add_worksheet("READ ME");
	package::gene::run($GENE,$config{"Report"},$Excel,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$Excel});
	package::sample::run($SAMPLE,\%format,\%library,\%annotation,\%analysis,\%filter,\%gender,$run{$Excel},\%hashassigngender);
	package::snv::run($SNV,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$Excel});
	package::genotype::run($Genotype,$Call,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$Excel});
	package::readme::run($Readme,\%format,"readme.txt");
}

