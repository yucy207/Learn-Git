BEGIN{
	# 加载主流程目录
	use File::Spec;
	my $path_curf = File::Spec->rel2abs(__FILE__);
	my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
	push @INC,"$dirs";
}
use strict;
use warnings;
use File::Spec;
use Bio::SeqIO;
use Bio::SearchIO;
use Statistics::R;
use Excel::Writer::XLSX;
use Parallel::ForkManager;
use package::utils;
use package::check;
use package::hwe;
use package::qc;
use package::check_methylation;
use package::mapping;
use package::target_typing;
use package::methyl_typing;
use package::report_methyl;
use package::mapping_bwa;
use package::gatk_VCF;
use package::report_snv;

my $tab="##########";

# # 读取配置
die "Usage: perl $0 ConfigFile\n" if(@ARGV!=1);
my $config_file = shift @ARGV;
my %hashPara    = package::utils::get_para();
my %hashConfig  = package::utils::read_config($config_file);

# # 开始流程
system "clear";
print "$tab Start NGS MethylTarget Pipline ".package::utils::get_time." $tab\n";
package::utils::show_config(\%hashConfig);
package::utils::is_continue();

my $input = "";
my $start = 1;
my $end   = "Finish"; # 设定流程的起始与终止，默认是完整流程
print "\n[Option] WhereStart?\nDefault:1\n1:check_methylation\n2:QC\n3:Mapping\n4:TargetTyping\n5:MethylTarget\n6:Report_Methyl\n7:BWA Mapping\n8:GATK VCF\n9:Report SNV\nInput:";$input=package::utils::get_input();$start=$input if($input=~/\d/);
print "\n[Option] WhereEnd?\nDefault:Finish\n1:check_methylation\n2:QC\n3:Mapping\n4:TargetTyping\n5:MethylTarget\n6:Report_Methyl\n7:BWA Mapping\n8:GATK VCF\n9:Report SNV\nInput:";$input=package::utils::get_input();$end=$input if($input=~/\d/ and $input<=9);

# # 1 check methylation info
package::check_methylation::run(\%hashPara, \%hashConfig) if($start<=1 and ($end eq "Finish" or $end>=1));

# # 2 QC
package::qc::run(\%hashPara, \%hashConfig) if($start<=2 and ($end eq "Finish" or $end>=2));

# # 3 Mapping
package::mapping::run(\%hashPara, \%hashConfig) if($start<=3 and ($end eq "Finish" or $end>=3));

# # 4 TargetTyping
package::target_typing::run(\%hashPara, \%hashConfig) if($start<=4 and ($end eq "Finish" or $end>=4));

# # 5 MethylTyping
package::methyl_typing::run(\%hashPara, \%hashConfig) if($start<=5 and ($end eq "Finish" or $end>=5));

# # 6 ReportMethylation
package::report_methyl::run(\%hashPara, \%hashConfig) if($start<=6 and ($end eq "Finish" or $end>=6));

# # 7 BWA Mapping
package::mapping_bwa::run(\%hashPara, \%hashConfig) if($start<=7 and ($end eq "Finish" or $end>=7));

# # 8 GATK VCF
package::gatk_VCF::run(\%hashPara, \%hashConfig) if($start<=8 and ($end eq "Finish" or $end>=8));

# # 9 Report SNV
package::report_snv::run(\%hashPara, \%hashConfig) if($start<=9 and ($end eq "Finish" or $end>=9));