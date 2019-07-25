$|=1;
use strict;
use warnings;
use File::Spec;
use Parallel::ForkManager;

die "Usage: perl $0 ConfigFile Suffix[final/sort]\n" if(@ARGV!=2);
my $configFile = shift @ARGV;
my $suffix     = shift @ARGV;
my %hashPara   = read_para("para.ini");
my %hashConfig = read_config($configFile);
# 分析

# 样本、目录准备
my @samples = get_sample(\%hashConfig, 'case', 'control');
my $report_dir = $hashConfig{'Report'};
my $output_dir = $hashConfig{'Output'};
my $status_dir = "$report_dir/status";
mkdir $status_dir if(not -e $status_dir);

# 分析文件准备检测
my $Threshold = $hashPara{'Process'}{'Metrics'};# 线程数
my @run_array=();
foreach my $sample(@samples){
	next if(checkFileSizeOK("$status_dir/$sample.$suffix.RmXY.metrics") and checkFileSizeOK("$status_dir/$sample.$suffix.RmXY.tag"));
	next if(not -e "$output_dir/$sample/$sample\_$suffix.bam");			
	push @run_array, $sample;
}

# 并行处理
my $target_bed      = $hashConfig{"TargetBed"};
my $target_bed_name = get_basename($target_bed);
my $target_bed_rmXY = "$report_dir/$target_bed_name.RmXY";
bed_remove_xy($target_bed, $target_bed_rmXY) if(scalar(@run_array) > 0 ); # bed去掉XY染色体

my $pm = Parallel::ForkManager->new($Threshold);
   $pm->run_on_start(sub{my ($pid, $sample) = @_; process_bar_array($sample, \@run_array)});# 进度条
foreach my $sample(@run_array)
{
    $pm->start($sample) and next;
    subRun(\%hashPara, \%hashConfig, $sample, $suffix, $target_bed_rmXY);
    $pm->finish;    
}
$pm->wait_all_children;

# 结果汇总显示
my %hashGender = read_gender_file("$report_dir/gender.log");
my @tag_files;
my @gender_values;
foreach my $sample(@samples)
{
	push @tag_files, "Sample=$sample=tag=$status_dir/$sample.$suffix.RmXY.tag";
	push @gender_values, "Sample=$sample=sex=$hashGender{$sample}";
}
print "" . (join "\n", @tag_files) . "\n";
print "" . (join "\n", @gender_values) . "\n";


sub subRun{
	my $hashPara        = shift @_;
	my $hashConfig      = shift @_;
	my $sample          = shift @_;
	my $suffix          = shift @_;
	my $target_bed_rmXY = shift @_;
 
	# 路径准备
	my $report_dir = $hashConfig->{'Report'};
	my $output_dir = $hashConfig->{'Output'};
	my $status_dir = "$report_dir/status";

	# 软件准备
	my $tmp_dir = $hashPara->{"Soft"}{"Tmp"};
	my $Java    = $hashPara->{"Soft"}{"Java"};
    my $Picard  = $hashPara->{"Soft"}{"Picard"};
	
	# 输入文件准备
	my $Species = $hashConfig->{"Species"};
 	my $genome  = $hashPara->{$Species}{"Genome"};
 	my $bam     = "$output_dir/$sample/$sample\_$suffix.bam";

	# 输出文件
	my $metrics = "$status_dir/$sample.$suffix.RmXY.metrics";
	my $tag     = "$status_dir/$sample.$suffix.RmXY.tag";  
	system "$Java -jar -Xmx3g $Picard CollectHsMetrics INPUT=$bam OUTPUT=$metrics R=$genome BI=$target_bed_rmXY TI=$target_bed_rmXY PER_TARGET_COVERAGE=$tag MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp_dir";
	
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
sub checkFileSizeOK{
	my $file=shift @_;
	if(-e $file and -s $file>0){
		return 1;
	}else{
		return 0;
	}
}
# 进度条， 根据输入向量计算
sub process_bar_array{
    my $process_name   = shift @_; # 需要分析的对象
    my $process_arrays = shift @_; # 总列表
    my $process_all_count = @$process_arrays;# 需要分析的总列表数量
    my ($pos) = grep $$process_arrays[$_] eq $process_name, 0..$process_all_count-1;
    my $process_count = $pos+ 1; # 当前对象的位置
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count; # 进度
    my $windows = 100; # 窗口宽度
    my $finished = $windows * int($process_perc) / 100; # 完成
    my $unfinished = $windows - $finished; # 未完成
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\n" if($process_count == $process_all_count);   
}

# 读参数
sub read_para{
	my $file = shift @_;
	my %hashPara = ();
	open FILE, $file or die "[ERR] Loss Para File,$file\n";
	while(my $line = <FILE>){
		$line =~ s/[\s\t\r\n]//g;
		my @split_line = split /\=/, $line;
		next if(@split_line != 3);
		$hashPara{$split_line[0]}{$split_line[1]} = $split_line[2];
	}
	close FILE;
	return %hashPara;
}

sub read_config{
	my $file = shift @_;
	my %hashConfig = ();
	open FILE, $file or die "[ERR] Loss Config File,$file\n";
	while(my $line = <FILE>){
		$line =~ s/[\s\t\r\n]//g;
		my @split_line = split /\=/,$line;
		next if(@split_line!=2);
		my ($name, $value) = ($split_line[0], $split_line[1]);

		if(exists $hashConfig{$name} and ($name eq "case" or $name eq "control") ){
			$hashConfig{$name} .= ",$value";
		}else{
			$hashConfig{$name} = $value;
		}
	}
	close FILE;
	return %hashConfig;
}
# 获取当前项目分析的样本
sub get_sample{
    my $hashConfig = shift @_;
    my @need_types = @_;
    my @samples;

    foreach my $need_type(@need_types)
    {
        next if(!exists($hashConfig->{$need_type}));
        foreach my $sample(split /,/, $hashConfig->{$need_type})
        {
            push @samples, $sample if($sample=~/\w/);
        }
    }

    return @samples;
}
# 获取路径文件名
sub get_basename{
    my $path = shift @_;
    my $path_curf = File::Spec->rel2abs($path);
    my ($vol, $dirs, $file) = File::Spec->splitpath($path_curf);
    return $file;
}
# 去除性染色体bed区域
sub bed_remove_xy{
	my $target_bed = shift @_;
	my $target_bed_rmXY = shift @_;
	open IN, $target_bed;
	open OUT,">$target_bed_rmXY";
	while(my $line = <IN>){
		my ($chr, $tmp) = split /\t/, $line;
		next if($chr eq 'X' or $chr eq 'Y');
		print OUT $line;
	}
	close IN;
	close OUT;
}