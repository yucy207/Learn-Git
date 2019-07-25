$|=1;
use strict;
use warnings;
use Encode;
use Excel::Writer::XLSX;
use Parallel::ForkManager;
use package::main;
use package::check;
use package::annotation_db;
use package::excavator2;
use package::split2target;
use package::single_sample_process;
use package::target_filter_combine;
use package::output;
use package::prepare_for_plot_n_snv;
use package::plot;

# # 读取配置
die "Usage perl $0 Config\n" if(@ARGV<1);
my $config_file = shift @ARGV;
my %hashPara    = package::main::get_para(); 
my %hashConfig  = package::main::read_config($config_file); 

# # R环境配置
my $RPATH      = $hashPara{"Soft"}{"R_dir"};
$ENV{'PATH'}   = "$RPATH:".$ENV{'PATH'};
my $RLIB       = $hashPara{"Soft"}{"RLib"};
$ENV{'R_LIBS'} = $RLIB;

my $input = "";
my $start = 1;
my $end   = "Finish"; # 设定流程的起始与终止，默认是完整流程
print "\n[Option] WhereStart?\nDefault:1\n1:excavator2\n2:split2target\n3:single_sample_process\n4:target_filter_combine\n5:output\n6:prepare_for_plot_n_snv\n7:plot\nInput:";$input=package::main::get_input();$start=$input if($input=~/\d/);
print "\n[Option] WhereEnd?\nDefault:Finish\n1:excavator2\n2:split2target\n3:single_sample_process\n4:target_filter_combine\n5:output\n6:prepare_for_plot_n_snv\n7:plot\nInput:";$input=package::main::get_input();$end=$input if($input=~/\d/ and $input<=6);

# # 1 excavator2运行 及 性别校正
package::excavator2::run(\%hashPara, \%hashConfig) if($start<=1 and ($end eq "Finish" or $end>=1));

# # 2 按target文件拆分区域
package::split2target::run(\%hashPara, \%hashConfig) if($start<=2 and ($end eq "Finish" or $end>=2));

# # 3 单样本区域注释 及 统计
package::single_sample_process::run(\%hashPara, \%hashConfig) if($start<=3 and ($end eq "Finish" or $end>=3));

# # 4 拆分区域的过滤 及 区域合并、注释
package::target_filter_combine::run(\%hashPara, \%hashConfig) if($start<=4 and ($end eq "Finish" or $end>=4));

# # 5 结果输出
package::output::run(\%hashPara, \%hashConfig) if($start<=5 and ($end eq "Finish" or $end>=5));

# # 6 为绘图和突变注释准备文件
package::prepare_for_plot_n_snv::run(\%hashPara, \%hashConfig) if($start<=6 and ($end eq "Finish" or $end>=6));

# # 7 绘图
package::plot::run(\%hashPara, \%hashConfig) if($start<=7 and ($end eq "Finish" or $end>=7));