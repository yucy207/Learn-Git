package package::check;
use strict;
use warnings;

# excavator2运行条件检测
sub repeat_LOH{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    # 检验线程与软件配置是否有问题
    my %hashParaCheck;
    my %hashConfigCheck;
	my $species = (exists $hashConfig->{'Species'}) ? $hashConfig->{'Species'} : 'Lost Species';
    $hashParaCheck{'Process'}{'repeat_LOH'}  = '\d';   # 数值型
    $hashParaCheck{'Soft'}{'AnnovarDIR'}     = '-dir'; # 路径型
    $hashParaCheck{'Soft'}{'bedtools'}       = '-file'; # 文件型
	
	$hashParaCheck{$species}{'Genome'}       = '-file';
	$hashParaCheck{$species}{'Dict'}         = '-file';
	$hashParaCheck{$species}{'AnnovarBuild'} = '\w';
	
    $hashConfigCheck{'Report'}     = '-dir';
    $hashConfigCheck{'Species'}    = '\w';
	
    my $isOK_para   = check_para($hashPara, \%hashParaCheck, 'necessary'); # 检查配置是否有问题
    my $isOK_config = check_config($hashConfig, \%hashConfigCheck); # 检查路径是否有问题
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
        $condition = 'OK' if(exists $hashConfig->{$content} and $hashConfigCheck->{$content} eq '\w' and $hashConfig->{$content} =~ /\w/);
        $condition = 'OK' if(exists $hashConfig->{$content} and $hashConfigCheck->{$content} eq '-file' and package::main::is_file_ok($hashConfig->{$content})); # 文件型配置
        $condition = 'OK' if(exists $hashConfig->{$content} and $hashConfigCheck->{$content} eq '-dir' and package::main::is_dir_ok($hashConfig->{$content}));# 路径型配置
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
            $condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '\w' and $hashPara->{$content}{$value} =~ /\w/);
			$condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '\d' and $hashPara->{$content}{$value}=~/^\d+$/); # 数值型配置
            $condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '-file' and package::main::is_file_ok($hashPara->{$content}{$value})); # 文件型配置
            $condition = 'OK' if(exists $hashPara->{$content}{$value} and $hashParaCheck->{$content}{$value} eq '-dir' and package::main::is_dir_ok($hashPara->{$content}{$value}));# 路径型配置
            system "echo 'Check Para $content -> $value\t\t\t$condition'"  if($condition eq 'OK');
            system "echo -e '\\033[41;37mCheck Para $content -> $value\t\t\t$condition \\033[0m'"  if($condition eq '[ERR]'); # 错误，红底白字
            system "echo -e '\\033[43;37mCheck Para $content -> $value\t\t\t$condition \\033[0m'"  if($condition eq '[Warning]'); # 警告，黄底白字
            $isOK = 0 if($condition ne 'OK');
        }
    }
    return $isOK;
}

1