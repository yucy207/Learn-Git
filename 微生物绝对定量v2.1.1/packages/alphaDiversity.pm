package alphaDiversity;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata     = shift;
	my $base         = shift;
	my $quantitation = shift;

	my @groups       = @{$metadata->{'group'}};
	my $group_num    = @groups;
	my $manager      = new Parallel::ForkManager($group_num);

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
		my $tmpPath = qq{$intm_result/tmp};
		system qq{mkdir -p $tmpPath} if not -d $tmpPath;

		if (uc($quantitation) eq "A" || uc($quantitation) eq "B"){
			Absolute_Quantitation($projectpath,$intm_result,$groupfile,$select,$metadata,$base);
		}
		
		if (uc($quantitation) eq "R" || uc($quantitation) eq "B"){
			Relative_Quantitation($projectpath,$intm_result,$groupfile,$select,$metadata,$base);
		}
		# system qq{touch $tmpPath/alpha.finish};
		print qq{$groupname alpha分析完成\n};
		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{alpha分析已经运行完成\n};
}


sub Absolute_Quantitation{

	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $mothur      = qq{$base->{alpha_mothur_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	my $tmpPath = qq{$intm_result/Absolute_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $AlphaDiversity = qq{$projectpath/Absolute_Quantitation/AlphaDiversity};
	if (-d $AlphaDiversity) {
		system "rm -rf $AlphaDiversity";
		system qq{mkdir -p $AlphaDiversity};
	}else{
		system qq{mkdir -p $AlphaDiversity};
	}


	my $otufile = qq{$tmpPath/otu.tax.$dissi.xls};                    ####因为样本水平有可能出现数值太大，R语言无法计算，所以用基于DNA水平计算
	my $otufile_rare = qq{$intm_result/Relative_Quantitation/tmp/otu.tax.$dissi.xls};  ####用于绘制稀释曲线

	### prepare output path
	my %path;
	my @function = ("DiversityIndex", "Rarefaction", "RankAbundance", "Specaccum");
	foreach (@function){

		my $str = $_."Path"; 
		$path{$str} = $AlphaDiversity."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}

	### prepare
	my $id = `id`;
	my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;

	my $num = `less $groupfile | cut -f 2 | sort -u | wc -l`;
	my $num_sample =  `less $groupfile | cut -f 1 | wc -l`;

	system qq{$Rscript $scriptdir/otu2shared.r $otufile_rare $groupfile $dissi $AlphaDiversity &>> $tmpPath/IRalpha.log};   ###用于绘制稀释曲线
	my $sharedfile = qq{$AlphaDiversity/sample.shared};

	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);


	# print "$select\n";
	for my $f (@function){

		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;

		if ($f eq "DiversityIndex"){
			# =============== DiversityIndex ===============
			system qq{$Rscript $scriptdir/Absolute_Abundance/16S.AlphaDiversity.DiversityIndex.IR.r $otufile $otufile_rare $groupfile $scriptdir "$select" $path{"DiversityIndexPath"} &>> $tmpPath/IRalpha.log};
		}
		
		if ($f eq "Rarefaction"){
			# =============== Rarefaction ===============

			my $RarefactionPath = $path{"RarefactionPath"};
			system qq{rm $RarefactionPath/*.xls} if -e qq{$RarefactionPath/*.xls};
			system qq{cp $sharedfile $RarefactionPath};
			my $sharedfile2 = qq{$RarefactionPath/sample.shared};
			die qq{ERROR:***mistakes when generate $sharedfile, Please Check It!: $!\n} if not (-e $sharedfile2 and -s $sharedfile2 != 0);
			system qq{$mothur bash -c 'mothur "#rarefaction.single(shared=$sharedfile2,label=$dissi,groupmode=f,processors=20)";chown -R $uid:$gid $RarefactionPath' &>> $tmpPath/IRalpha.log};
			system qq{$Rscript $scriptdir/16S.AlphaDiversity.RarefactionCurves.r $RarefactionPath 40 $RarefactionPath &>> $tmpPath/IRalpha.log};
			
			my @files = `ls $path{"RarefactionPath"} | grep '.rabund\$\\|.rarefaction\$'`;
			foreach my $file (@files){

				chomp $file;
				system qq{mv $path{RarefactionPath}/$file $path{RarefactionPath}/$file.xls};
			
			}
			
		}
		
		if ($f eq "RankAbundance"){
		# 	# =============== RankAbundance ===============
			system qq{$Rscript $scriptdir/Absolute_Abundance/16S.AlphaDiversity.RankAbundance_IR.r $otufile $groupfile $path{"RankAbundancePath"} &>> $tmpPath/IRalpha.log};
		}
		
		if ($f eq "Specaccum"){
		# 	# =============== pecaccum ===============
			system qq{$Rscript $scriptdir/16S.AlphaDiversity.Specaccum.r $otufile $groupfile $path{"SpecaccumPath"} &>> $tmpPath/IRalpha.log};
		}
		
		$pm->finish; # Terminates the child process

	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

	system qq{rm -rf $sharedfile};

	if ($num_sample < 3) {
		system qq{rm -rf $AlphaDiversity/Specaccum};
		delete $path{"SpecaccumPath"};
	}


	### check: tell you Where is the error!
	foreach my $key (keys %path){

		my $errorpath = qq{$path{$key}};
		my $minsize = `ls -lR $errorpath|grep "pdf"|awk '{print \$5}'| sort -n |head -n 1`;
		$key =~ s/Path//;

		print qq{WARINING/ERROR:**** mistakes when generate pdf file into $errorpath, 
		Please check it using 16S.AlphaDiversity.$key.r *** $!\n} if $minsize < 4500;   # exists & size>5000
	}

}

sub Relative_Quantitation{

	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $mothur      = qq{$base->{alpha_mothur_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	my $tmpPath = qq{$intm_result/Relative_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $AlphaDiversity = qq{$projectpath/Relative_Quantitation/AlphaDiversity};
	if (-d $AlphaDiversity) {
		system "rm -rf $AlphaDiversity";
		system qq{mkdir -p $AlphaDiversity};
	}else{
		system qq{mkdir -p $AlphaDiversity};
	}

	my $otufile     = qq{$tmpPath/subsample_otu.tax.$dissi.xls};
	my $num = `less $groupfile | cut -f 2 | sort -u | wc -l`;
	my $num_sample =  `less $groupfile | cut -f 1 | wc -l`;
	### prepare output path
	my %path;
	my @function = ("DiversityIndex", "Shannon", "Rarefaction", "RankAbundance", "Specaccum");
	foreach (@function){

		my $str = $_."Path"; 
		$path{$str} = $AlphaDiversity."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}

	### prepare
	my $id = `id`;
	my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;

	system qq{$Rscript $scriptdir/otu2shared.r $otufile $groupfile $dissi $AlphaDiversity &>> $tmpPath/alpha.log};
	my $sharedfile = qq{$AlphaDiversity/sample.shared};

	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	for my $f (@function){

		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;

		if ($f eq "DiversityIndex"){
			# =============== DiversityIndex ===============
			system qq{$Rscript $scriptdir/16S.AlphaDiversity.DiversityIndex.r $otufile $groupfile "$select" $path{"DiversityIndexPath"} &>> $tmpPath/alpha.log};
		}
		
		if ($f eq "Shannon"){
			# =============== Shannon ===============
			my $ShannonPath = $path{"ShannonPath"};
			system qq{rm $ShannonPath/*.xls} if -e qq{$ShannonPath/*.xls};
			system qq{cp $sharedfile $ShannonPath};
			my $sharedfile1 = qq{$ShannonPath/sample.shared};
			die qq{ERROR:***mistakes when generate $sharedfile, Please Check It!: $!\n} if not (-e $sharedfile1 and -s $sharedfile1 != 0);
			system qq{$mothur bash -c 'mothur "#rarefaction.single(shared=$sharedfile1,label=$dissi,calc=shannon,groupmode=f,processors=20)";chown -R $uid:$gid $ShannonPath' &>> $tmpPath/alpha.log};
			system qq{$Rscript $scriptdir/16S.AlphaDiversity.ShannonCurves.r $ShannonPath 40 $ShannonPath &>> $tmpPath/alpha.log};

			my @files = `ls $path{"ShannonPath"} | grep '.rabund\$\\|.r_shannon\$'`;			
			foreach my $file (@files){
			
				chomp $file;
				system qq{mv $path{ShannonPath}/$file $path{ShannonPath}/$file.xls};
			
			}
			
		}
		
		if ($f eq "Rarefaction"){
			# =============== Rarefaction ===============

			my $RarefactionPath = $path{"RarefactionPath"};
			system qq{rm $RarefactionPath/*.xls} if -e qq{$RarefactionPath/*.xls};
			system qq{cp $sharedfile $RarefactionPath};
			my $sharedfile2 = qq{$RarefactionPath/sample.shared};
			die qq{ERROR:***mistakes when generate $sharedfile, Please Check It!: $!\n} if not (-e $sharedfile2 and -s $sharedfile2 != 0);
			system qq{$mothur bash -c 'mothur "#rarefaction.single(shared=$sharedfile2,label=$dissi,groupmode=f,processors=20)";chown -R $uid:$gid $RarefactionPath' &>> $tmpPath/alpha.log};
			system qq{$Rscript $scriptdir/16S.AlphaDiversity.RarefactionCurves.r $RarefactionPath 40 $RarefactionPath &>> $tmpPath/alpha.log};
			
			my @files = `ls $path{"RarefactionPath"} | grep '.rabund\$\\|.rarefaction\$'`;
			foreach my $file (@files){

				chomp $file;
				system qq{mv $path{RarefactionPath}/$file $path{RarefactionPath}/$file.xls};
			
			}
			
		}
		
		if ($f eq "RankAbundance"){
			# =============== RankAbundance ===============
			system qq{$Rscript $scriptdir/16S.AlphaDiversity.RankAbundance.r $otufile $groupfile $path{"RankAbundancePath"} &>> $tmpPath/alpha.log};
		}
		
		if ($f eq "Specaccum"){
			# =============== pecaccum ===============
			system qq{$Rscript $scriptdir/16S.AlphaDiversity.Specaccum.r $otufile $groupfile $path{"SpecaccumPath"} &>> $tmpPath/alpha.log};
		}
		
		$pm->finish; # Terminates the child process

	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

	system qq{rm -rf $sharedfile};

	if ($num_sample < 3) {
		system qq{rm -rf $AlphaDiversity/Specaccum};
		delete $path{"SpecaccumPath"};
	}
	### check: tell you Where is the error!
	foreach my $key (keys %path){

		my $errorpath = qq{$path{$key}};
		my $minsize = `ls -lR $errorpath|grep "pdf"|awk '{print \$5}'| sort -n |head -n 1`;
		$key =~ s/Path//;

		print qq{WARINING/ERROR:**** mistakes when generate pdf file into $errorpath, 
		Please check it using 16S.AlphaDiversity.$key.r *** $!\n} if $minsize < 4500;   # exists & size>5000
	}

}

1;
