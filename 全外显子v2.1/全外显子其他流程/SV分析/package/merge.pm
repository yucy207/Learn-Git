package package::merge;
use strict;
use warnings;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start vcf merge ".package::utils::get_time()." ##########\n";
    my $report_dir  = $hashConfig->{'Report'}; 
    my $output_dir  = $hashConfig->{'Output'};
    my $SVdbtype    = exists $hashConfig->{'GenomeBed'} ? "SVdbWGS" : "SVdbWES"; #判断是WES 还是WGS    
	my $svReportDir = "$report_dir/sv";
	my $tmp_dir     = "$svReportDir/TMP";
	mkdir $svReportDir if(not -e $svReportDir);

    # 数据状态检测
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control');
    my %hashCondition; 
    
	# 已完成该步骤
    my $sv_vcf_merged = "$svReportDir/sample.all.final.sv.vcf";
	if(package::utils::is_file_ok($sv_vcf_merged)) {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }
	foreach my $sample(@samples) {
	    my $sampel_tmp_dir       = "$tmp_dir/$sample";
		my $sample_readdepth_txt = "$sampel_tmp_dir/$sample.sv.$sample\_final.bam.readdepth.txt";
		my $sample_readdepth_bed = "$sampel_tmp_dir/$sample.sv.$sample\_final.bam.readdepth.bed";
		my $sample_sv_vcf        = "$sampel_tmp_dir/$sample.sv.vcf";
		
        # 原始数据没问题
        ##if(package::utils::is_file_ok($sample_readdepth_txt, $sample_readdepth_bed, $sample_sv_vcf)) {
		if(package::utils::is_file_ok($sample_sv_vcf)) {
            ##$hashCondition{"Good2Run"}{$sample} = "$sample_readdepth_txt $sample_readdepth_bed $sample_sv_vcf";
			$hashCondition{"Good2Run"}{$sample} = "$sample_sv_vcf";
            next;
        }
        # 原始数据丢失
        $hashCondition{"Error"}{$sample} = "$sample_sv_vcf";
    }

    # 写日志
    package::utils::write_log("$report_dir/run.log", "SV_merge", \%hashCondition); 
	die "[ERR]LOSS sv.vcf, please check run.log\n" if(exists $hashCondition{"Error"}); # # 有样本没有结果即停止
    
	# # 数据库样本获取及去重
	my %sample_project;
	map{$sample_project{$_} ++} @samples;
	my $Species   = $hashConfig->{'Species'};
	my %hashSVDB  = package::utils::GetSVDBFiles($hashConfig, $hashPara, \%sample_project, $SVdbtype);
	my @db_sv_vcfs = split /,/, $hashSVDB{'sv_vcfs'};
	
	# # 生成vcf列表
    my @sv_vcfs = map{"$tmp_dir/$_/$_.sv.vcf"} @samples;
	push @sv_vcfs, @db_sv_vcfs;
	my $sv_vcf_list = join "\n", @sv_vcfs;
	my $list_file = "$svReportDir/sv_vcf_list.txt";
	open LIST, ">$list_file";
	print LIST "$sv_vcf_list\n";
	close LIST;
	
	# # 进行样本vcf合并
    my $svtools   = $hashPara->{'Soft'}{'svtools'};
	my $tmp_dir2  = $hashPara->{'Soft'}{'Tmp'};
	my $lsort_vcf = "$svReportDir/sample.lsort.sv.vcf";
	system("$svtools lsort -r -f $list_file -b 20 -t $tmp_dir2 > $lsort_vcf");
	system("$svtools lmerge -g -i $lsort_vcf > $sv_vcf_merged");
}

1