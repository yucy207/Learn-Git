package package::gender;
use strict;
use warnings;

# 计算样本性别
sub run{
    my $hashConfig = shift @_;
    print "########## Start Gender ".package::utils::get_time()." ##########\n";
    my $report_dir  = $hashConfig->{'Report'}; 
    my $gender_file = exists $hashConfig->{"GenderFile"} ? $hashConfig->{"GenderFile"} : '';
    my $gender_log  = "$report_dir/gender.log";
    my @samples     = package::utils::get_sample($hashConfig, 'case', 'control');
    #####
    # 写日志
    #####
    package::utils::write_log_simple("$report_dir/run.log", "Gender");

    # 性别判定
    # 1.自定义性别
    my %hashGender_Defined = read_gender_file($gender_file) if($gender_file=~/\w/ and package::utils::is_file_ok($gender_file) == 1);
    # 2.log已存在
    my %hashGender_Log     = read_gender_file($gender_log) if(package::utils::is_file_ok($gender_log) == 1);
    # 3. log是否已经包含了所有样本
    my $is_sample_has_gender = 1;
    foreach my $sample(@samples)
    {
        $is_sample_has_gender = 0 if(!exists $hashGender_Log{$sample});
    }
    if($is_sample_has_gender == 1 and package::utils::is_force_sample($hashConfig) == 0)
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }
    # 4. 计算样本性别
    my %hashGender;
    open GENDER, ">$gender_log";
    foreach my $sample(@samples)
    {   
        print "Process $sample\n";
        if(exists $hashGender_Defined{$sample}) # 自定义
        {   
            print GENDER "# $sample gender defined in GenderFile\n"; 
            $hashGender{$sample} = $hashGender_Defined{$sample};
        }
        else
        {
            my $tag = "$report_dir/status/$sample.tag"; 
            my ($y_cov, $x_cov, $ydividx) = calculate_y_cov($tag);
            my $gender = 'female'; # 默认女性
            if($y_cov=~/\d/ and $x_cov=~/\d/ and $ydividx=~/\d/)
            {
                $gender = 'male'   if($ydividx > 0.4);
                $gender = 'female' if($ydividx < 0.1);
                print GENDER "$sample\t$y_cov/$x_cov=$ydividx\n"; # x/y覆盖情况判定样本性别
            }
            else
            {
                print GENDER "$sample NO X/Y chromosome in bed! default defined as female\n" if($x_cov eq 'X_NO_COV' or $y_cov eq 'Y_NO_COV'); # X/Y没有数据
            }
            $hashGender{$sample} = $gender;                
        }
    }
    map{print GENDER "$_\t$hashGender{$_}\n";}@samples;
    close GENDER;
 

}

sub calculate_y_cov{
    my $tag = shift @_;
    my %hashData;
    my ($y_cov, $x_cov, $ydividx) = (0, 0, 0.5);
    return ($y_cov, $x_cov, 'LostTag') if(package::utils::is_file_ok($tag) == 0);
    open TAG, $tag;
    while(<TAG>)
    {
        $_=~s/[\r\n]//g;
        my ($chr, $start, $end, $length, $name, $gc, $mean_cov, $norm_cov, $tmp) = split /\t/, $_, 9;
        $chr =~s/^chr//;
        next if($chr ne 'X' and $chr ne 'Y');
        $hashData{$chr}{'Sum'} += $norm_cov;
        $hashData{$chr}{'Count'}++;;
    }
    close TAG;
    $y_cov = $hashData{'Y'}{'Sum'}/$hashData{'Y'}{'Count'} if(exists $hashData{'Y'});
    $x_cov = $hashData{'X'}{'Sum'}/$hashData{'X'}{'Count'} if(exists $hashData{'X'});
    $ydividx = $y_cov/$x_cov if($x_cov > 0);
    $y_cov = 'Y_NO_COV' if(!exists $hashData{'Y'});
    $x_cov = 'X_NO_COV' if(!exists $hashData{'X'});
    return ($y_cov, $x_cov, $ydividx);
    
}
# 读性别文件
sub read_gender_file{
    my $gender_file = shift @_;

    my %hashGender;
    open GENDER, $gender_file;
    while(<GENDER>)
    {
        $_=~s/[\r\n]//g;
        next if($_!~/\w/ or $_=~/^#/);
        my ($sample, $gender) = split /\t/, $_;
        next if(!defined $gender or ($gender ne 'female' and $gender ne 'male'));
        $hashGender{$sample} = $gender;
    }
    close GENDER;
    return %hashGender;
}

1