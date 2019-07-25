$|=1;
use strict;
use warnings;
use Statistics::R;
use Excel::Writer::XLSX;
use Parallel::ForkManager;
use package::utils;
use package::sv;
use package::merge;
use package::annotation;
use package::annotation_db;
use package::homology;
use package::filter;
use package::circos;
use package::region_plot;
use package::evidence;


my $tab="##########";


# # 读取配置
die "Usage: perl $0 ConfigFile\n" if(@ARGV!=1);
my $config_file = shift @ARGV;
my %hashPara    = package::utils::get_para();
my %hashConfig  = package::utils::read_config($config_file);

# # 开始流程
system "clear";
print "$tab Start SV Pipline ".package::utils::get_time." $tab\n";
package::utils::show_config(\%hashConfig);
package::utils::is_continue();

my $input = "";
my $start = 1;
my $end   = "Finish"; # 设定流程的起始与终止，默认是完整流程
print "\n[Option] WhereStart?\nDefault:1\n1:sv\n2:merge\n3:annotation\n4:Homology\n5:Excel and Ciocos Output\n6:Region Plot\n7:bam evidence\nInput:";$input=package::utils::get_input();$start=$input if($input=~/\d/);
print "\n[Option] WhereEnd?\nDefault:Finish\n1:sv\n2:merge\n3:annotation\n4:Homology\n5:Excel and Ciocos Output\n6:Region Plot\n7:bam evidence\nInput:";$input=package::utils::get_input();$end=$input if($input=~/\d/ and $input<=8);


# # 1 sv
package::sv::run(\%hashPara, \%hashConfig) if($start<=1 and ($end eq "Finish" or $end>=1));

# # 2 merge
package::merge::run(\%hashPara, \%hashConfig) if($start<=2 and ($end eq "Finish" or $end>=2));

# # 3  注释
package::annotation::run(\%hashPara, \%hashConfig) if($start<=3 and ($end eq "Finish" or $end>=3));

# # 4 断点同源性检测
package::homology::run(\%hashPara, \%hashConfig) if($start<=4 and ($end eq "Finish" or $end>=4));

# # 5 filter and Result output (excel and circos)
package::filter::run(\%hashPara, \%hashConfig) if($start<=5 and ($end eq "Finish" or $end>=5));

# # 6 Region plot 需要用到cnvseq的结果
package::region_plot::run(\%hashPara, \%hashConfig) if($start<=6 and ($end eq "Finish" or $end>=6));

# # 7 bam evidence
package::evidence::run(\%hashPara, \%hashConfig) if($start<=7 and ($end eq "Finish" or $end>=7));