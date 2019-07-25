package package::check;
use strict;
use warnings;

# check_methylation运行条件检测
sub check_methylation{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck;
    my %hashConfigCheck;
    my %hashParaCheck_warning;
    my $species = package::utils::get_species($hashConfig);
    $hashParaCheck{'Soft'}{'BlastPlus'}  = '-dir';

    $hashParaCheck_warning{$species}{'Genome'}      = '-file'; 
    $hashParaCheck_warning{$species}{'RefGene'}     = '-file'; 
    $hashParaCheck_warning{$species}{'SpecialMark'} = '-file'; 
 
    $hashConfigCheck{'Report'}       = '-dir';
    $hashConfigCheck{'TargetFasta'}  = '-file';
    $hashConfigCheck{'Primer'}       = '-file';
    $hashConfigCheck{'AnalysisType'} = '\w';
    $hashConfigCheck{'Species'}      = '\w';
    my $isOK_para   = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题
    my $isOK_para_w = check_para($hashPara, \%hashParaCheck_warning, 'warning');# 检查配置是否有问题，只给出警告
    my $isOK_config = check_config($hashConfig, \%hashConfigCheck); # 检查路径是否有问题
    my $isOK = ($isOK_para==1 and $isOK_config==1);
    return $isOK;
}

# QC运行条件检测
sub qc{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck;
    my %hashConfigCheck;
    $hashParaCheck{'Process'}{'QC'}     = '\d'; # 数值型
    $hashParaCheck{'Soft'}{'FastQC'}    = '-file';# 文件型
    $hashParaCheck{'Soft'}{'FastqStat'} = '-file';
    $hashParaCheck{'Soft'}{'R'}         = '-file';
    $hashParaCheck{'Soft'}{'RLib'}      = '-dir';# 路径
    $hashConfigCheck{'Fastq'}  = '-dir';
    $hashConfigCheck{'Report'} = '-dir';
    my $isOK_para   = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题
    my $isOK_config = check_config($hashConfig, \%hashConfigCheck); # 检查路径是否有问题
    my $isOK = ($isOK_para==1 and $isOK_config==1);
    return $isOK;
}

# mapping运行条件检测
sub mapping{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck; # 必须存在的参数
    my %hashParaCheck_warning; # 非必须
    my %hashConfigCheck;# 必须
    my $species = package::utils::get_species($hashConfig);
    $hashParaCheck{'Process'}{'Mapping'} = '\d'; # 数值型
    $hashParaCheck{'Soft'}{'FastQC'}     = '-file';# 文件型
    $hashParaCheck{'Soft'}{'FastqStat'}  = '-file';    
    $hashParaCheck{'Soft'}{'BlastPlus'}  = '-dir';
    $hashParaCheck{'Soft'}{'FastX'}      = '-dir'; 
    $hashParaCheck{'Soft'}{'FLASH'}      = '-file';

    $hashConfigCheck{'Fastq'}        = '-dir';
    $hashConfigCheck{'Output'}       = '-dir';
    $hashConfigCheck{'Report'}       = '-dir';
    $hashConfigCheck{'TargetFasta'}  = '-file';
    $hashConfigCheck{'AnalysisType'} = '\w';
    $hashConfigCheck{'Species'}      = '\w';
    my $isOK_para    = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题,这些参数必须要给出的
    my $isOK_config  = check_config($hashConfig, \%hashConfigCheck); # 检查是否有问题
    my $isOK = ($isOK_para==1 and $isOK_config==1);
    return $isOK;
}

# target_typing运行条件检测
sub target_typing{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck; # 必须存在的参数
    my %hashParaCheck_warning; # 非必须
    my %hashConfigCheck;# 必须
    my $species = package::utils::get_species($hashConfig);
    $hashParaCheck{'Process'}{'TargetTyping'} = '\d'; # 数值型

    $hashConfigCheck{'Output'}       = '-dir';
    $hashConfigCheck{'Report'}       = '-dir';
    my $isOK_para    = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题,这些参数必须要给出的
    my $isOK_config  = check_config($hashConfig, \%hashConfigCheck); # 检查是否有问题
    my $isOK = ($isOK_para==1 and $isOK_config==1);
    return $isOK;
}

# methyl_typing运行条件检测
sub methyl_typing{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck; # 必须存在的参数
    my %hashParaCheck_warning; # 非必须
    my %hashConfigCheck;# 必须
    my $species = package::utils::get_species($hashConfig);
    $hashParaCheck{'Process'}{'MethylTyping'} = '\d'; # 数值型

    $hashConfigCheck{'Output'}       = '-dir';
    $hashConfigCheck{'Report'}       = '-dir';
    $hashConfigCheck{'TargetFasta'}  = '-file';
    my $isOK_para    = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题,这些参数必须要给出的
    my $isOK_config  = check_config($hashConfig, \%hashConfigCheck); # 检查是否有问题
    my $isOK = ($isOK_para==1 and $isOK_config==1);
    return $isOK;
}

sub report_methyl{
    my $hashConfig = shift @_;
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control', 'Ycontrol');
    my %hashConfigCheck;# 必须
    $hashConfigCheck{'Report'} = '-dir';
    my $isOK = check_config($hashConfig, \%hashConfigCheck); # 检查是否有问题
    return $isOK;
}

sub mapping_bwa{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck; # 必须存在的参数
    my %hashParaCheck_warning; # 非必须
    my %hashConfigCheck;# 必须
    my $species = package::utils::get_species($hashConfig);
    $hashParaCheck{'Process'}{'MappingBWA'} = '\d'; # 数值型
    $hashParaCheck{'Soft'}{'Java'}          = '-file';# 文件型
    $hashParaCheck{'Soft'}{'BWA'}           = '-file';    
    $hashParaCheck{'Soft'}{'SamTools'}      = '-file';
    $hashParaCheck{'Soft'}{'Picard'}        = '-file'; 
    $hashParaCheck{'Soft'}{'Tmp'}           = '-dir';

    $hashConfigCheck{'Output'}       = '-dir';
    $hashConfigCheck{'Report'}       = '-dir';
    my $isOK_para    = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题,这些参数必须要给出的
    my $isOK_config  = check_config($hashConfig, \%hashConfigCheck); # 检查是否有问题
    my $isOK = ($isOK_para==1 and $isOK_config==1);
    return $isOK;
}

# gatk_GVCF运行条件检测
sub gatk_VCF{
    my $hashPara = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck; # 必须存在的参数
    my %hashParaCheck_warning; # 非必须
    my %hashConfigCheck;# 必须
    my $species = package::utils::get_species($hashConfig);
    $hashParaCheck{'Process'}{'GATK_VCF'} = '\d'; # 数值型
    $hashParaCheck{'Soft'}{'Tmp'}         = '-dir'; 
    $hashParaCheck{'Soft'}{'Java'}        = '-file'; 
    $hashParaCheck{'Soft'}{'GATK4_Loc'}   = '-file'; 
    $hashParaCheck{'Soft'}{'bcftools'}    = '-file';
    $hashParaCheck{'Soft'}{'Bgzip'}       = '-file';

    $hashConfigCheck{'Report'} = '-dir';
    my $isOK_para    = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题,这些参数必须要给出的
    my $isOK_para_w  = check_para($hashPara, \%hashParaCheck_warning, 'warning');# 检查配置是否有问题，只给出警告
    my $isOK_config  = check_config($hashConfig, \%hashConfigCheck); # 检查是否有问题
    my $isOK = ($isOK_para==1 and $isOK_config==1);
    return $isOK;
}
 

# 检查配置文件中的路径是否有问题
sub check_config{
    my $hashConfig = shift @_;
    my $hashConfigCheck = shift @_;
    my $isOK = 1;
    foreach my $content(sort keys %$hashConfigCheck)
    {
        
        my $condition = '[ERR]';
        $condition = 'OK' if(exists $hashConfig->{$content} and $hashConfigCheck->{$content} eq '\w' and $hashConfig->{$content}=~/\w/);
        $condition = 'OK' if(exists $hashConfig->{$content} and $hashConfigCheck->{$content} eq '-file' and package::utils::is_file_ok($hashConfig->{$content})); # 文件型配置
        $condition = 'OK' if(exists $hashConfig->{$content} and $hashConfigCheck->{$content} eq '-dir' and package::utils::is_dir_ok($hashConfig->{$content}));# 路径型配置
        system "echo 'Check Config $content \t\t\t$condition'" if($condition eq 'OK');
        system "echo -e '\\033[41;37mCheck Config $content \t\t\t$condition \\033[0m'" if($condition eq '[ERR]'); # 错误，红底白字
        $isOK = 0 if($condition ne 'OK');       
    }
    return $isOK;
}

# 检查Para是否有问题
sub check_para{
    my $hashPara      = shift @_;
    my $hashParaCheck = shift @_;
    my $need_level    = shift @_;
    my $isOK = 1;
    foreach my $content(sort keys %$hashParaCheck)
    {
        foreach my $value(sort keys %{$hashParaCheck->{$content}})
        {
            my $condition = ($need_level eq 'warning') ? '[Warning]' : '[ERR]';
            $condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '\w' and $hashPara->{$content}{$value}=~/\w/); # 字符型配置
            $condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '\d' and $hashPara->{$content}{$value}=~/^\d+$/); # 数值型配置
            $condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '-file' and package::utils::is_file_ok($hashPara->{$content}{$value})); # 文件型配置
            $condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '-dir' and package::utils::is_dir_ok($hashPara->{$content}{$value}));# 路径型配置
            system "echo 'Check Para $content -> $value\t\t\t$condition'"  if($condition eq 'OK');
            system "echo -e '\\033[41;37mCheck Para $content -> $value\t\t\t$condition \\033[0m'"  if($condition eq '[ERR]'); # 错误，红底白字
            system "echo -e '\\033[43;37mCheck Para $content -> $value\t\t\t$condition \\033[0m'"  if($condition eq '[Warning]'); # 警告，黄底白字
            $isOK = 0 if($condition ne 'OK');
        }
    }
    return $isOK;
}

1