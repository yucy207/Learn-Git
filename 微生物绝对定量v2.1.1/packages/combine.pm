package combine;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata    = shift;
	my $base        = shift;

	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);

	for my $groupfile (@groups){

		my @groupname = split /\//,$groupfile;
		# my ($groupname) = $groupname[-1] =~ /(.*)\.[\w+]/;
		my ($groupname) = ($groupname[-1] =~ /::/)? $groupname[-1] =~ /(.*)\.\w\w\w::/ : $groupname[-1] =~ /(.*)\.[\w+]/;
		my ($groupfile, $select) = split /::/, $groupfile;
		my $intm_result = qq{$metadata->{intm_result}/$groupname};
		my $projectpath = qq{$metadata->{result}/Report/$groupname};

		die "ERROR: There is no $groupfile, Please Check It! $!\n" if not -e $groupfile;

		# Forks and returns the pid for the child:
		$manager->start and next;

		### Combine分析
		Combine($intm_result,$projectpath,$groupfile,$metadata,$base);

		# system qq{touch $intm_result/tmp/combine.finish};
		print qq{$groupname Combine_Analysis分析已经运行完成\n};

		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{Combine_Analysis分析已经运行完成\n};
}


sub Combine
{
	
	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $kronatools  = qq{$base->{kronatools_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $Rscript_xlsx  = qq{$base->{Rscript_xlsx}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $Combine = qq{$projectpath/Combine_Analysis};
	system qq{mkdir -p $Combine} if not -d $Combine;

	# return if -e qq{$tmpPath/combine.finish};


	my $RQ_dir = qq{$projectpath/Relative_Quantitation};
	my $AQ_dir = qq{$projectpath/Absolute_Quantitation};

	my $RQ_dir_int = qq{$intm_result/Relative_Quantitation};
	my $AQ_dir_int = qq{$intm_result/Absolute_Quantitation};

	die "ERROR: There is no $RQ_dir, Please Check It! $!\n" if not -d $RQ_dir;
	die "ERROR: There is no $AQ_dir, Please Check It! $!\n" if not -d $AQ_dir;
	die "ERROR: There is no $RQ_dir_int, Please Check It! $!\n" if not -d $RQ_dir_int;
	die "ERROR: There is no $AQ_dir_int, Please Check It! $!\n" if not -d $AQ_dir_int;
	# die "ERROR: There is no $groupfile, Please Check It! $!\n" if not -e $groupfile;

	# prepare output path
	my %path;
	foreach("Community", "AlphaDiversity", "BetaDiversity", "PICRUSt"){

		my $str = $_."Path"; 
		$path{$str} = $Combine."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}
	my $num = `less $groupfile | cut -f 2 | sort -u | wc -l`;
	my $samplenum = `less $groupfile | cut -f 1 | sort -u | wc -l`;
	
	###============= Community ==========
	# ###============================ Braplot ===============
	if ($samplenum == 1) {
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Combine.Community.Barplot_one.r $RQ_dir/Community $AQ_dir/Community $groupfile $path{"CommunityPath"} &> $tmpPath/combine.log};
	}else{
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Combine.Community.Barplot.r $RQ_dir/Community $AQ_dir/Community $groupfile $path{"CommunityPath"} &>> $tmpPath/combine.log};
	}

	# # # # ###============================ ANOVA ===============
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Combine.Community.ANOVA.r $RQ_dir/Community $AQ_dir/Community $groupfile $path{"CommunityPath"} $scriptdir &>> $tmpPath/combine.log};

	# # # # ###============================ Kruskal_Wallis ===============
	system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Combine.Community.Kruskal_Wallis.boxplot.r $RQ_dir/Community $AQ_dir/Community $groupfile $path{"CommunityPath"} $scriptdir &>> $tmpPath/combine.log};

	# # # # ###============================ Wilcoxon ===============
	system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Combine.Community.Wilcoxon.boxplot.r $RQ_dir/Community $AQ_dir/Community $groupfile $path{"CommunityPath"} $scriptdir &>> $tmpPath/combine.log};


	# # # ###============= AlphaDiversity ==========
	# # # ###============================ Index ===============
	system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Combine.Alpha.Index.boxplot.r $RQ_dir/AlphaDiversity $AQ_dir/AlphaDiversity $groupfile $path{"AlphaDiversityPath"} &>> $tmpPath/combine.log};	



	# # # ###============= BetaDiversity ==========
	# # # ###============================ ADONIS ===============
	system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Combine.Beta.ADONIS.Distance.r $RQ_dir_int $AQ_dir_int $groupfile $path{"BetaDiversityPath"} &>> $tmpPath/combine.log};

	# # # ###===========================Lefse========================
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Combine.Beta.lefse.r $RQ_dir/BetaDiversity $AQ_dir/BetaDiversity $groupfile $path{"BetaDiversityPath"} $scriptdir &>> $tmpPath/combine.log};



	# ###============= PICRUSt ==========
	# ###============================ t-test ===============

	my $RQ_check = `ls $RQ_dir/PICRUSt/pair_group_different_analysis/* | grep "_diff.P0.05.xls" |wc -l`;
	my $AQ_check = `ls $AQ_dir/PICRUSt/pair_group_different_analysis/* | grep "_diff.P0.05.xls" |wc -l`;
	if ($RQ_check > 0 and $AQ_check > 0) {
		system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Combine.PI_pair_diff_analysis.r $RQ_dir/PICRUSt $AQ_dir/PICRUSt $groupfile COG $path{"PICRUStPath"} $scriptdir &>> $tmpPath/combine.log};
		system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Combine.PI_pair_diff_analysis.r $RQ_dir/PICRUSt $AQ_dir/PICRUSt $groupfile KEGG $path{"PICRUStPath"} $scriptdir &>> $tmpPath/combine.log};		
	}

	# # # ###============================ ANOVA ===============
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Combine.PI_all_diff_analysis.r $RQ_dir/PICRUSt $AQ_dir/PICRUSt $groupfile COG $path{"PICRUStPath"} $scriptdir &>> $tmpPath/combine.log};
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.Combine.PI_all_diff_analysis.r $RQ_dir/PICRUSt $AQ_dir/PICRUSt $groupfile KEGG $path{"PICRUStPath"} $scriptdir &>> $tmpPath/combine.log};



	###检查结果
	my $Wilcoxonresult = `ls $path{"CommunityPath"}/Wilcoxon|wc -l`;
	if ($Wilcoxonresult == 0) {
		system qq{rm -rf $path{"CommunityPath"}/Wilcoxon};
	}

	my $Kruskal_Wallis_result = `ls $path{"CommunityPath"}/Kruskal_Wallis|wc -l`;
	if ($Kruskal_Wallis_result == 0) {
		system qq{rm -rf $path{"CommunityPath"}/Kruskal_Wallis};
	}

	if ($num < 3) {
		system qq{rm -rf $path{"CommunityPath"}/Kruskal_Wallis};
		system qq{rm -rf $path{"CommunityPath"}/ANOVA};
	}

	my $PI_result = `ls $path{"PICRUStPath"} |wc -l`;
	system qq{rm -rf $path{"PICRUStPath"}} if $PI_result == 0;

}

1;
