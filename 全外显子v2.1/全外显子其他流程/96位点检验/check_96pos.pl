# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Parallel::ForkManager;
use Excel::Writer::XLSX;
use Encode;

# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
my $gatk4  = "/home/genesky/software/gatk/4.0.12.0/gatk-package-4.0.12.0-local.jar";
my $java   = "/home/genesky/software/java/1.8.0_181/bin/java";
my $genome = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.fa"; #基因组固化
my %hashDB = database_list(); # SNP数据库。

# 检测 -> 脚本输入
my ($config_file, $snp_file, $output_name, $database, $database_list, $threshold, $if_help);
GetOptions(
    "config_file|c=s"  => \$config_file,
    "snp_file|s=s"     => \$snp_file,
    "output_name|p=s"  => \$output_name,
    "database|d=s"     => \$database,
    "threshold|t=s"    => \$threshold,
    "database_list|l!" => \$database_list,
    "help|h"           => \$if_help,
);
 
die show_database_list(\%hashDB) if(defined $database_list); # 查看数据库
die help() if(defined $if_help or (not defined $config_file or not defined $snp_file));

$threshold   = 10 if(not defined $threshold);
$output_name = '96SNP_QC' if(not defined $output_name);
$database    = '20190515' if(not defined $database);
die "not find database\n" if(not exists $hashDB{$database}{'bed_96'} or not exists $hashDB{$database}{'vcf_in'} or not -e $hashDB{$database}{'bed_96'} or not -e $hashDB{$database}{'vcf_in'}); # 数据库不存在

my $bed_96 = $hashDB{$database}{'bed_96'}; # 96位点区域文件
my $vcf_in = $hashDB{$database}{'vcf_in'}; # 96位点vcf文件

###################################################################### 主程序
 
# (1) 数据读入
my @samples;
my %config   = read_config($config_file); # 读取 配置文件
my %hash_vcf = read_snpscan($snp_file,\@samples); # 读取 SNPSCAN 分型结果

# (2) 二代96位点分型
if(@samples > 0)
{
    my $pm = Parallel::ForkManager->new($threshold);
    foreach my $sample(@samples)
    {
        $pm->start and next;
        vcf_to_genotype(\%config, $sample);
        $pm->finish;    
    }
    $pm->wait_all_children;
}
else
{
    print "[Note] Process None, for reProcess Delete Result\n";
}

read_vcf(\%config, \%hash_vcf, \@samples); # 读取 二代分型结果

# (3) 二者分型结果统计
my %hash_tongji = check_geno(\%hash_vcf, \@samples);

# (4) 结果输出
output(\%config, \%hash_vcf, \%hash_tongji, $output_name);

###################################################################### 子函数
sub database_list {
    my %hashDB;

    $hashDB{'20190413'}{'bed_96'} = SCRIPTDIR."database/20190413/96_pos_region_20190413.bed"; # 96位点区域文件固化
    $hashDB{'20190413'}{'vcf_in'} = SCRIPTDIR."database/20190413/96_pos_region_20190413.vcf"; # #96位点vcf文件固化
	
	$hashDB{'20190515'}{'bed_96'} = SCRIPTDIR."database/20190515/96_pos_region_20190515.bed"; # 96位点区域文件固化
    $hashDB{'20190515'}{'vcf_in'} = SCRIPTDIR."database/20190515/96_pos_region_20190515.vcf"; # #96位点vcf文件固化

    return %hashDB;
}

# # 展示数据库列表
sub show_database_list {
    my $hashDB = shift @_;
    print "\nSupport DB Version (Named after date) : \n\n";
    foreach my $version(sort keys %$hashDB)
    {   
	    print "\t$version\n";
	    foreach my $file(sort keys %{$hashDB->{$version}})
		{
            print "\t\t$file:$hashDB->{$version}{$file}\n";
		}
		print "\n";
    }
    return "\n";
}

sub write_readme{
    my $workbook = shift @_;
    my $format   = shift @_;
    my $file     = shift @_;
    open IN, $file;
    my $readme   = $workbook->add_worksheet("Read Me");
    my $row = 0;
    $readme->set_column(0,0,15);
    $readme->set_column(1,1,20);
    $readme->set_column(2,2,60);
    while(<IN>){
        $_       =~s/[\r\n]//g;
        my @data = split /\t/, $_;
        my $type = $row == 0 ? "title" : "readme_me";
        $readme->write($row, 0, decode("gb2312", $data[0]),$format->{$type});
        $readme->write($row, 1, decode("gb2312", $data[1]),$format->{$type});
        $readme->write($row, 2, decode("gb2312", $data[2]),$format->{$type});
        $row++;        
    }
    close IN;
}

# # 结果输出
sub output{
    my $config      = shift @_;
    my $hash_vcf    = shift @_;
    my $hash_tongji = shift @_;
    my $output_name = shift @_;
    my $document    = $config->{'Report'}."/document/";
    my $qc_dir      = "$document/1_QC";
    mkdir $document if(not -e $document);
    mkdir $qc_dir if(not -e $qc_dir);
    my $workbook    = Excel::Writer::XLSX->new("$qc_dir/$output_name.xlsx");
    my %format      = format_run($workbook); 
    ######
    # 输出 sheet1：Genotyping Data
    ######
    my $geno_sheet  = $workbook->add_worksheet("Genotyping Data");
    $geno_sheet->set_row(0, 60);
    $geno_sheet->set_column(1,1,15);
    $geno_sheet->set_column(6,8,15);
    my $row = 0;
    my @titles_1    = ('Samples','rs_ID','Ref Allele','Alt Allele','Ref_dp','Alt_dp','Genotype Calling','Genotype Check','SNPscan Result');
    foreach my $col(0..@titles_1-1){
        $geno_sheet->write($row, $col, $titles_1[$col], $format{"title"});
    }
    foreach my $sample(sort {$a cmp $b} keys %{$hash_vcf}) {
        foreach my $rs(sort {$a cmp $b} keys %{$hash_vcf->{$sample}}) {
            $row ++;
            foreach my $col(0..@titles_1-1){
                $geno_sheet->write($row, $col, $hash_vcf->{$sample}{$rs}{$titles_1[$col]}, $format{"normal"});
            }
        }
    }    
    ######
    # 输出 sheet2：DeepSeq Quality Control
    ######
    my $qual_sheet  = $workbook->add_worksheet("DeepSeq Quality Control");
    $qual_sheet->set_row(0, 60);
    $qual_sheet->set_column(6,7,18);
    $row = 0;
    my @titles_2    = ('Samples','Total SNPs','Miss','Cons','SeqError','CallError','Genotyping Accuracy%','False Negative%');
    foreach my $col(0..@titles_2-1){
        $qual_sheet->write($row, $col, $titles_2[$col], $format{"title"});
    }
    foreach my $sample(sort {$a cmp $b} keys %{$hash_tongji}) {
        $row ++;
        foreach my $col(0..@titles_2-1){
            $qual_sheet->write($row, $col, $hash_tongji->{$sample}{$titles_2[$col]}, $format{"normal"});
        }
    }
    
    my $flag = 0;
    my $string = "Warnings : These samples' consistency < 90% : ";
    foreach my $sample(sort {$a cmp $b} keys %{$hash_tongji}){
        my $value = (exists $hash_tongji->{$sample}{"Genotyping Accuracy"}) ? $hash_tongji->{$sample}{"Genotyping Accuracy"} : 0;
        if($value < 0.9){
            $string.="$sample, ";
            $flag = 1;
        }
    }
    system "echo -e \"\\033[41;37m $string \\033[0m\"" if($flag == 1); # 警告 红底白字
    
    ######
    # 输出 sheet3：Read Me
    ######
    write_readme($workbook,\%format,"readme_check_96pos.txt");
}

# # 一代、二代分型结果统计
sub check_geno {
    my $hash_vcf = shift @_;
    my $sample   = shift @_;
    my %hash_tongji;
    foreach my $sample(keys %{$hash_vcf}) {
        my ($total,$miss,$cons,$seq_error,$call_error) = (0,0,0,0,0);    
        foreach my $rs(keys %{$hash_vcf->{$sample}}) {
            $total ++;
            if(not exists $hash_vcf->{$sample}{$rs}{'Genotype Calling'} or $hash_vcf->{$sample}{$rs}{'Genotype Calling'} !~/\w/ or not exists $hash_vcf->{$sample}{$rs}{'SNPscan Result'} or $hash_vcf->{$sample}{$rs}{'SNPscan Result'} !~/\w/) {
                $hash_vcf->{$sample}{$rs}{'Genotype Check'} = 'Miss';
                $miss ++;
            }
            else {
                my ($allele_1,$allele_2) = split /\//, $hash_vcf->{$sample}{$rs}{'Genotype Calling'};
                if("$allele_1/$allele_2" eq $hash_vcf->{$sample}{$rs}{'SNPscan Result'} or "$allele_2/$allele_1" eq $hash_vcf->{$sample}{$rs}{'SNPscan Result'}) {
                    $hash_vcf->{$sample}{$rs}{'Genotype Check'} = 'Cons';
                    $cons ++;
                }
                elsif(exists $hash_vcf->{$sample}{$rs}{'Genotype Check'} and $hash_vcf->{$sample}{$rs}{'Genotype Check'} eq 'CallError') {
                    $call_error ++;
                }
                else{
                    $hash_vcf->{$sample}{$rs}{'Genotype Check'} = 'SeqError';
                    $seq_error ++;
                }
            }
        }
        $hash_tongji{$sample}{'Samples'}              = $sample;
        $hash_tongji{$sample}{'Total SNPs'}           = $total;
        $hash_tongji{$sample}{'Miss'}                 = $miss;
        $hash_tongji{$sample}{'Cons'}                 = $cons;
        $hash_tongji{$sample}{'SeqError'}             = $seq_error;
        $hash_tongji{$sample}{'CallError'}            = $call_error;
        $hash_tongji{$sample}{'Genotyping Accuracy'}  = ($cons / ($total - $miss));
        $hash_tongji{$sample}{'Genotyping Accuracy%'} = sprintf "%.2f%%", ($cons / ($total - $miss)) * 100;
        $hash_tongji{$sample}{'False Negative%'}      = sprintf "%.2f%%", (($call_error / 2) / ($total - $miss)) * 100;
    }
    return %hash_tongji;
}

# # 读取二代测序结果
sub read_vcf { 
    my $config     = shift @_;
    my $hash_vcf   = shift @_;
    my $samples    = shift @_;
    my $report_dir = $config->{'Report'};
	my %hash_rs_trans = trans_rs();
    foreach my $sample(@{$samples}) {
        my $vcf_dir    = "$report_dir/check/$sample\_filtered.vcf";
        open VCF, $vcf_dir;
        while(<VCF>) {
            next if($_ =~/^#/);
            $_ =~s/[\r\n]//g;
            my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$format_info) = split /\t/, $_;
            my $mark = "$chr\_$pos\_$ref";
            my $rs   = $hash_rs_trans{$mark};
            ## test_filter($sample, $mark, $info); #验证gatk的filter结果
            my @formats      = split /:/, $format;
            my @format_infos = split /:/, $format_info;        
            my %tmp;
            foreach my $i(0..@formats-1) {
                $tmp{$formats[$i]} = $format_infos[$i];
            }
            $tmp{'GT'} =~s/0/$ref/g;
            $tmp{'GT'} =~s/1/$alt/g;
            my ($ref_dp,$alt_dp) = exists $tmp{'AD'} ? split /,/, $tmp{'AD'} : ("","");
            $hash_vcf->{$sample}{$rs}{'Samples'}          = $sample;
            $hash_vcf->{$sample}{$rs}{'rs_ID'}            = $rs;
            $hash_vcf->{$sample}{$rs}{'Ref Allele'}       = $ref;
            $hash_vcf->{$sample}{$rs}{'Alt Allele'}       = $alt;
            $hash_vcf->{$sample}{$rs}{'Ref_dp'}           = $ref_dp;
            $hash_vcf->{$sample}{$rs}{'Alt_dp'}           = $alt_dp;
            $hash_vcf->{$sample}{$rs}{'Genotype Quality'} = exists $tmp{'GQ'} ? $tmp{'GQ'} : "";
            $hash_vcf->{$sample}{$rs}{'Genotype Calling'} = $tmp{'GT'};    
            $hash_vcf->{$sample}{$rs}{'Genotype Check'}   = 'CallError' if($filter =~/CallError/); #若位点未通过筛选条件，标记为CallError
        }
        close VCF;
    }
}

# # 二代96位点分型
sub vcf_to_genotype{
    my $config = shift @_;
    my $sample = shift @_;
    
    my $output_dir = $config->{'Output'};
    my $report_dir = $config->{'Report'};
    my $check_dir  = "$report_dir/check";
    mkdir $check_dir if(not -e $check_dir);
    
    my $final_bam   = "$output_dir/$sample/$sample\_final.bam";
    my $raw_vcf     = "$check_dir/$sample.vcf";
    my $filter_vcf  = "$check_dir/$sample\_filtered.vcf";
    my $filter_mark = "CallError";
    system("$java -jar $gatk4 HaplotypeCaller -R $genome -L $bed_96 --alleles $vcf_in --min-base-quality-score 20 --genotyping-mode GENOTYPE_GIVEN_ALLELES --output-mode EMIT_ALL_SITES -I $final_bam -O $raw_vcf");
    system("$java -jar $gatk4 VariantFiltration -R $genome -V $raw_vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name $filter_mark -O $filter_vcf");

}

# # 读取SNPSCAN结果
sub read_snpscan { 
    my $snp_file = shift @_;
    my $samples  = shift @_;
    my %hash_vcf;
    open SNP, $snp_file;
    my $line1  = <SNP>;
    $line1     =~ s/[\r\n]//g;
    my @titles = split /\t/, $line1;
    while(<SNP>) {
        next if($_ !~/\w/);
        $_         =~ s/[\r\n]//g;
        my @infos  = split /\t/, $_;
        my $sample = $infos[0];
        push @{$samples}, $sample;
        foreach my $i(1..@titles-1) {
            $hash_vcf{$sample}{$titles[$i]}{'SNPscan Result'} = $infos[$i];
            $hash_vcf{$sample}{$titles[$i]}{'Samples'}        = $sample;
            $hash_vcf{$sample}{$titles[$i]}{'rs_ID'}          = $titles[$i];
        }
    }
    close SNP;
    return %hash_vcf;
}

# # 读取配置文件
sub read_config {
    my $config_dir = shift @_;
    my %config;
    open FILE, $config_dir;
    while(<FILE>) {
        $_ =~s/[\r\n\s]//g;
        next if($_ =~ /^#/);
        next if($_ !~ /\w/);
        my ($key,$value) = split /=/, $_;
        $config{$key}    = $value;
    }
    close FILE;
    return %config;
}

sub format_run{
    my ($workbook)=@_;
    my %format=();
    $format{'title'} = $workbook->add_format();
    $format{'title'} ->set_align('center');
    $format{'title'} ->set_align('vcenter');
    $format{'title'} ->set_size(12);
    $format{'title'} ->set_font("Times New Roman");
    $format{'title'} ->set_border();
    $format{'title'} ->set_bg_color("yellow");
    $format{'title'} ->set_color("black");
    
    $format{'normal'} = $workbook->add_format();
    $format{'normal'} ->set_align('center');
    $format{'normal'} ->set_align('vcenter');
    $format{'normal'} ->set_size(12);
    $format{'normal'} ->set_font("Times New Roman");
    $format{'normal'} ->set_border();
    
    $format{'readme_me'} = $workbook->add_format();
    $format{'readme_me'} ->set_align('vcenter');
    $format{'readme_me'} ->set_size(12);
    $format{'readme_me'} ->set_font("Times New Roman");
    $format{'readme_me'} ->set_border();
    
    $format{'small'} = $workbook->add_format();
    $format{'small'} ->set_align('vcenter');
    $format{'small'} ->set_size(10);
    $format{'small'} ->set_font("Times New Roman");
    $format{'small'} ->set_border();
    
    $format{'seq'} = $workbook->add_format();
    $format{'seq'} ->set_align('vcenter');
    $format{'seq'} ->set_size(11);
    $format{'seq'} ->set_font("Courier New");
    $format{'seq'} ->set_border();
    
    $format{'left'} = $workbook->add_format();
    $format{'left'} ->set_align('vcenter');
    $format{'left'} ->set_size(12);
    $format{'left'} ->set_font("Times New Roman");
    $format{'left'} ->set_border();
    
    $format{'orange'} = $workbook->add_format();
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();

    $format{'skyblue'} = $workbook->add_format();
    $format{'skyblue'} ->set_align('vcenter');
    $format{'skyblue'} ->set_size(12);
    $format{'skyblue'} ->set_font("Times New Roman");
    $format{'skyblue'} ->set_bg_color("#538ed5");
    $format{'skyblue'} ->set_border();

    $format{'bold'} = $workbook->add_format( bold => 1 );
    $format{'blue'} = $workbook->add_format( color => "#538ed5" );
    $format{'redbold'} = $workbook->add_format( color => "#ff0000", bold => 1, );
    $format{'italic'} = $workbook->add_format( italic => 1 );
    $format{'boldblue'} = $workbook->add_format( bold => 1, color => "#538ed5" );
    $format{'bolditalic'} = $workbook->add_format( bold => 1, italic => 1 );
    $format{'blueitalic'} = $workbook->add_format( color => "#538ed5", italic => 1 );
    $format{'boldblueitalic'} = $workbook->add_format( bold => 1, color => "#538ed5", italic => 1 );
    
    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    
    $format{'readme1'} = $workbook->add_format();
    $format{'readme1'}->set_align('center');
    $format{'readme1'}->set_align('vcenter');
    $format{'readme1'}->set_bold();
    $format{'readme1'}->set_size(14);
    $format{'readme1'}->set_font("Times New Roman");
    $format{'readme1'}->set_border();

    $format{'readme2'} = $workbook->add_format();
    $format{'readme2'}->set_align('vcenter');
    $format{'readme2'}->set_bold();
    $format{'readme2'}->set_size(14);
    $format{'readme2'}->set_font("Times New Roman");

    $format{'readme2tmp'} = $workbook->add_format();
    $format{'readme2tmp'}->set_right();

    $format{'readme3'} = $workbook->add_format();
    $format{'readme3'}->set_align('center');
    $format{'readme3'}->set_align('vcenter');
    $format{'readme3'}->set_bold();
    $format{'readme3'}->set_size(11);
    $format{'readme3'}->set_font("Times New Roman");
    $format{'readme3'}->set_border();

    $format{'readme4'} = $workbook->add_format();
    $format{'readme4'}->set_align('vcenter');
    $format{'readme4'}->set_bold();
    $format{'readme4'}->set_size(11);
    $format{'readme4'}->set_font("Times New Roman");
    $format{'readme4'}->set_border();

    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    return %format;
}

sub trans_rs{
    my %hash_rs_trans;
	open VCF, $vcf_in;
	while(<VCF>) {
	    $_ =~ s/[\r\n]//g;
		next if($_ =~ /^#/);
		my ($chr, $pos, $id, $ref, $alt, $others) = split /\t/, $_;
		$hash_rs_trans{"$chr\_$pos\_$ref"}    = $id;
	}
	close VCF;
    return %hash_rs_trans;
}

sub help {
    my $info = "
Program: check_96pos， 96位点验证
Version: 2019-05-15
Contact: 239 田也

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:

         --config_file/-c    [必填] 基础流程的配置文件。
         --snp_file/-s       [必填] 实验提供的一代测序SNP分型文件，97列，第一列样本名，其他列为需要验证的96个位点，必须有表头，表头为SNP的rs号。
		                           （一般直接复制实验给出表格，注意样本名需要与基础分析的样本名保持一致，不能有空格等标点存在）
         --output_name/-p    [选填] 输出文件名称。默认：96SNP_QC
         --database/-d       [选填] SNP数据库版本。默认：20190515
         --threshold/-t      [选填] 并行数。默认：10
         --database_list/-l         列出支持的SNP数据库版本及版本对应的文件。
         --help/-h                  查看帮助文档。
    \n";
    return $info;
}
