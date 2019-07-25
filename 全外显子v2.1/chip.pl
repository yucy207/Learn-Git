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
use Statistics::R;
use Parallel::ForkManager;
use package::utils;
use package::check;
use package::qc;
use package::mapping;
use package::gatk_GVCF;
use package::gatk_VCF;
use package::annotation_db;
use package::library;
use package::database;
use package::annotation;
use package::gender;
use package::snv_check;
use package::format;

my $tab="##########";

# # 读取配置
die "Usage: perl $0 ConfigFile\n" if(@ARGV!=1);
my $config_file = shift @ARGV;
my %hashPara    = package::utils::get_para();
my %hashConfig  = package::utils::read_config($config_file);

# # 开始流程
system "clear";
print "$tab Start NGS Chip Pipline ".package::utils::get_time." $tab\n";
package::utils::show_config(\%hashConfig);
package::utils::is_continue();

my $input = "";
my $start = 1;
my $end   = "Finish"; # 设定流程的起始与终止，默认是完整流程
print "\n[Option] WhereStart?\nDefault:1\n1:QC\n2:Mapping\n3:gatk_GVCF\n4:gatk_VCF\n5:library\n6:gender\n7:database\n8:annotation\n9:snv_check\nInput:";$input=package::utils::get_input();$start=$input if($input=~/\d/);
print "\n[Option] WhereEnd?\nDefault:Finish\n1:QC\n2:Mapping\n3:gatk_GVCF\n4:gatk_VCF\n5:library\n6:gender\n7:database\n8:annotation\n9:snv_check\nInput:";$input=package::utils::get_input();$end=$input if($input=~/\d/ and $input<=8);

# # 1 QC
package::qc::run(\%hashPara, \%hashConfig) if($start<=1 and ($end eq "Finish" or $end>=1));

# # 2 Mapping
package::mapping::run(\%hashPara, \%hashConfig) if($start<=2 and ($end eq "Finish" or $end>=2));

# # 3 GATK GVCF [GVCF+VCF 4全外，13G = 4小时]
package::gatk_GVCF::run(\%hashPara, \%hashConfig) if($start<=3 and ($end eq "Finish" or $end>=3));

# # 4 GATK VCF
package::gatk_VCF::run(\%hashPara, \%hashConfig) if($start<=4 and ($end eq "Finish" or $end>=4));

# # 5 Library
package::library::run(\%hashPara, \%hashConfig) if($start<=5 and ($end eq "Finish" or $end>=5));

# # 6 Gender
package::gender::run(\%hashConfig) if($start<=6 and ($end eq "Finish" or $end>=6));

# # 7 Database
package::database::run(\%hashConfig) if($start<=7 and ($end eq "Finish" or $end>=7));

# # 8 Annotation
package::annotation::run(\%hashPara, \%hashConfig) if($start<=8 and ($end eq "Finish" or $end>=8));

# # 9 SNV check
package::snv_check::run(\%hashPara, \%hashConfig) if($start<=9 and ($end eq "Finish" or $end>=9));

