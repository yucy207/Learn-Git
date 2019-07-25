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

# ��ȡ����
die "Usage perl $0 ConfigFile RunConfig OutFileName\n" if(@ARGV!=3);
my $configFile=shift @ARGV;
my $runConfig=shift @ARGV;
my $outfile=shift @ARGV;
my %config=package::main::readConfig($configFile);
my %run=package::main::readRunConfig($runConfig);

# ��ʼ����
system "clear";
print "$tab Annotation ".package::time::ymd." ".package::time::hms." $tab\n";
die "[ERR] Loss Report DIR\n" if(not exists $config{"Report"} or not -e $config{"Report"});
die "[ERR] Loss Library\n" if(not -e $config{"Report"}."/database.snv");
die "[ERR] Loss Annotation\n" if(not -e $config{"Report"}."/annotation.txt");
die "[ERR] Loss Analysis\n" if(not -e $config{"Report"}."/analysis.txt");
die "[ERR] Loss Gender\n" if(not -e $config{"Report"}."/gender.log");
die "[ERR] Loss ReadMe File\n" if(not -e "readme.txt");
my %library=package::main::readLibrary($config{"Report"}."/database.snv");
my %annotation=package::main::readAnnotation($config{"Report"}."/annotation.txt");
my %analysis=package::main::readAnalysis($config{"Report"}."/analysis.txt");
my %gender=package::gender::read($config{"Report"}."/gender.log");
my %hashassigngender=package::gender::read2(\%config);#�Զ����Ա����ȫ����ͨ��tag�ж��Ա�
mkdir $config{"Report"}."/document" if(not -e $config{"Report"}."/document");
mkdir $config{"Report"}."/document/2_SNV" if(not -e $config{"Report"}."/document/2_SNV");
# Excel
print $config{"Report"}."/document/2_SNV/$outfile.xlsx\n";
my $workbook=Excel::Writer::XLSX->new($config{"Report"}."/document/2_SNV/$outfile.xlsx");
my %format=package::format::run($workbook);
my %tmp=();
my %genesnv=();
my $OneExcel="";
foreach my $Excel(sort keys %run){
	print "Process $Excel \n";
    $OneExcel=$Excel;
	my %filterResult=package::filter::run(\%library,\%annotation,\%analysis,$run{$Excel});
	my @ids=keys %filterResult;
	foreach my $id(@ids){
	    my @titles=keys %{$filterResult{$id}};
		foreach my $title(@titles){
		    my @alts=keys %{$filterResult{$id}{$title}};
			foreach my $alt(@alts){
			    my @userinfos=keys %{$filterResult{$id}{$title}{$alt}};
				foreach my $userinfo(@userinfos){
				    next if($userinfo eq 'SNV NO.');
				    if($userinfo eq 'gene'){
					   my $gene=$filterResult{$id}{$title}{$alt}{$userinfo};
					   $genesnv{$gene}{$title}{$alt}++;
					}else{
					    $tmp{$title}{$alt}{$userinfo}=$filterResult{$id}{$title}{$alt}{$userinfo};
					}
				}
			}
		}
	}	
}
my %filter=getfilter(\%tmp,\%genesnv);
my $GENE=$workbook->add_worksheet("Gene");
my $SAMPLE=$workbook->add_worksheet("Sample");
my $SNV=$workbook->add_worksheet("SNV Information");
my $Genotype=$workbook->add_worksheet("GENOTYPES");
my $Call=$workbook->add_worksheet("GENOTYPE CALLING");
my $Readme=$workbook->add_worksheet("READ ME");
package::gene::run($GENE,$config{"Report"},$outfile,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$OneExcel});
package::sample::run($SAMPLE,\%format,\%library,\%annotation,\%analysis,\%filter,\%gender,$run{$OneExcel},\%hashassigngender);
package::snv::run($SNV,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$OneExcel});
package::genotype::run($Genotype,$Call,\%format,\%library,\%annotation,\%analysis,\%filter,$run{$OneExcel});
package::readme::run($Readme,\%format,"readme.txt");

sub getfilter{
    my $tmp=shift @_;
	my $genesnv=shift @_;
	my @genes=keys %$genesnv;
	my $id=1;
	my %filter=();
	foreach my $gene(@genes){
	    my @titles=keys %{$genesnv->{$gene}};
		foreach my $title(@titles){
		    my @alts=keys %{$genesnv->{$gene}{$title}};
			foreach my $alt(@alts){
			    $filter{$id}{$title}{$alt}{"SNV NO."}="SNV"."0"x(5-length($id)).$id;
				my @userinfos=keys %{$tmp->{$title}{$alt}};
				foreach my $userinfo(@userinfos){
				    $filter{$id}{$title}{$alt}{$userinfo}=$tmp->{$title}{$alt}{$userinfo};
				}
				$id++;
			}
		}
	}
    return %filter;
}
