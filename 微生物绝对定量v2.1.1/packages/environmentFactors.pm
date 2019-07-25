package environmentFactors;
use Parallel::ForkManager;
use strict;

sub run

{

	my $metadata      = shift;
	my $base          = shift;
	my $quantitation  = shift;
	# my $quantitation  = qq{$metadata->{quantitation}};

	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);

	for my $groupfile (@groups){
		my @groupname     = split /\//,$groupfile;
		my ($groupname) = ($groupname[-1] =~ /::/)? $groupname[-1] =~ /(.*)\.\w\w\w::/ : $groupname[-1] =~ /(.*)\.[\w+]/;
		my ($groupfile, $select) = split /::/, $groupfile;		
		# my ($groupname)   = $groupname[-1] =~ /(.*)\.[\w+]/;
		my $intm_result   = qq{$metadata->{intm_result}/$groupname};
		my $projectpath   = qq{$metadata->{result}/Report/$groupname};

		die "ERROR: *** There is no $groupfile, please check it: $!***\n" if not -e $groupfile;
	
		# Forks and returns the pid for the child:
		$manager->start and next;

		my $tmpPath = qq{$intm_result/tmp};
		system qq{mkdir -p $tmpPath} if not -d $tmpPath;

		if (uc($quantitation) eq "A" || uc($quantitation) eq "B"){
			Absolute_Quantitation($projectpath,$intm_result,$groupfile,$metadata,$base);
		}		

		if (uc($quantitation) eq "R" || uc($quantitation) eq "B"){
			Relative_Quantitation($projectpath,$intm_result,$groupfile,$metadata,$base);
		}

		# system qq{touch $tmpPath/env.finish};
		print qq{$groupname env分析完成\n};

		$manager->finish;

	}

	$manager->wait_all_children;
	print qq{env分析已经运行完成\n};

}



# sub rm_finish
# {

# 	my $intm_result = shift;
# 	my $groupname   = shift;

# 	if (-e qq{$intm_result/Absolute_Quantitation/tmp/env.finish}){

# 		print qq{$groupname Absolute_Quantitation IRenv分析已经运行完成!\n 是否重新运行? [y/n]};
# 		my $option = <STDIN>;
# 		$option =~ s/[\r\n]//g;

# 		if ($option eq "y"){
# 			system qq{rm $intm_result/Absolute_Quantitation/tmp/env.finish};
# 			system qq{rm $intm_result/Absolute_Quantitation/tmp/env.log} if -e qq{$intm_result/Absolute_Quantitation/tmp/env.log};
# 		}

# 	}



# 	if (-e qq{$intm_result/Relative_Quantitation/tmp/env.finish}){
# 		print qq{$groupname Relative_Quantitation env分析已经运行完成!\n 是否重新运行? [y/n]};
# 		my $option = <STDIN>;
# 		$option =~ s/[\r\n]//g;
# 		if ($option eq "y"){

# 			system qq{rm $intm_result/Relative_Quantitation/tmp/env.finish};
# 			system qq{rm $intm_result/Relative_Quantitation/tmp/env.log} if -e qq{$intm_result/Relative_Quantitation/tmp/env.log};

# 		}

# 	}

# }



sub Absolute_Quantitation{

	### prepare Absolute_Quantitation ###
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi     = qq{$metadata->{dissi}};
	my $envpath   = qq{$metadata->{env}};
	my $Rscript   = qq{$base->{Rscript_bin}};
	my $vegan     = qq{$base->{vegan_bin}};
	my $scriptdir = qq{$base->{scriptdir}/sub_script};


	my $EnvironmentFactorsPath = qq{$projectpath/Absolute_Quantitation/EnvironmentFactors};
	system qq{mkdir -p $EnvironmentFactorsPath} if not -d $EnvironmentFactorsPath;
	my $tmpPath = qq{$intm_result/Absolute_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $otufile     = qq{$tmpPath/otu.tax.$dissi.xls};
	my $taxonpath   = qq{$projectpath/Absolute_Quantitation/Community/Community_Structure};

	# return if -e qq{$tmpPath/env.finish};

	### Run 
	my $envfile     = qq{$projectpath/env.xls};
	system qq{$Rscript $scriptdir/16S.EnvironmentFactors.select.r $envpath $groupfile $envfile};
	# system qq{cp $envfile $projectpath/env.xls};

	my %path;
	my @function = ("BioEnv", "MantelTest", "RDA_CCA", "Correlation");
	foreach (@function){
	    my $str = $_."Path"; 
		$path{$str} = $EnvironmentFactorsPath."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}

	my @taxonfiles = `find $taxonpath -name *.Abundance.xls`;
	die "ERROR: *** There is no files named *.taxon.Abundance.xls in $taxonpath, please ckeck it (may rename files?) ***\n" if scalar @taxonfiles < 1;

	#### use functions
	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	for my $f (@function){
		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;
		if ($f eq "BioEnv"){
			# =============== BioEnv 筛选与群落相关性最大的环境因子组合 ===============
			system qq{$vegan Rscript $scriptdir/16S.EnvironmentFactors.BioEnv.r $envfile $otufile $path{"BioEnvPath"} &>> $tmpPath/env.log};
		}

		if ($f eq "MantelTest"){
			# =============== Mental Test 检验样本群落分布与环境因子整体之间的相关性 ===============
			system qq{$Rscript $scriptdir/16S.EnvironmentFactors.MantelTest.r $envfile $otufile $path{"MantelTestPath"} &>> $tmpPath/env.log};
		}	

		if ($f eq "RDA_CCA"){
			foreach my $taxonfile (@taxonfiles) {
				chomp $taxonfile;
				my ($taxon) = $taxonfile =~ /\/(\w+)\.[^\/]+$/;
				# ========== RDA_CCA ===============
				my $RDA_CCAPath = qq{$path{RDA_CCAPath}/$taxon};
				system qq{mkdir -p $RDA_CCAPath} if not -d $RDA_CCAPath;
				system qq{$vegan Rscript $scriptdir/16S.EnvironmentFactors.RDA_CCA.r $envfile $taxonfile $groupfile $RDA_CCAPath &>> $tmpPath/env.log};
				system qq{$vegan Rscript $scriptdir/RDA_degree.R $RDA_CCAPath/sample.position.csv $RDA_CCAPath/environmentfactors.position.csv $RDA_CCAPath/sample_with_env_relationship.xls &>> $tmpPath/env.log};
		    	system qq{$vegan Rscript $scriptdir/RDA_degree.R $RDA_CCAPath/species.position.csv $RDA_CCAPath/environmentfactors.position.csv $RDA_CCAPath/species_with_env_relationship.xls &>> $tmpPath/env.log};

			}

		}

		if ($f eq "Correlation"){
			foreach my $taxonfile (@taxonfiles) {
				chomp $taxonfile;
				my ($taxon) = $taxonfile =~ /\/(\w+)\.[^\/]+$/;
				# ========== Correlation ===============
				my $CorrelationPath = qq{$path{CorrelationPath}/$taxon};
				system qq{mkdir -p $CorrelationPath} if not -d $CorrelationPath;
				system qq{$vegan Rscript $scriptdir/16S.EnvironmentFactors.Correlation.r $envfile $taxonfile 100 $CorrelationPath &>> $tmpPath/env.log};
			}
		}

		$pm->finish; # Terminates the child process

	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

	### check: tell you Where is the error!
	foreach my $p (keys %path){

		my $errorpath = $path{$p};
		$p =~ s/Path//;

		if ($p eq "BioEnv" or $p eq "MantelTest"){
			my $file = `ls $errorpath`;
			chomp $file;
			my $filelinesnum = `less $errorpath/$file | wc -l`;
			print qq{ERROR:**** mistakes when generate result into $errorpath, 
			Please check it using 16S.EnvironmentFactors.$p.r *** $!\n} if $filelinesnum < 2;    # 结果行数小于2
		}

		
		if ($p eq "RDA_CCA" or $p eq "Correlation"){

			opendir DIR,$errorpath or die "Can not open $errorpath\n";
			my @files = grep {-d "$errorpath/$_" && ! /^\.{1,2}$/} readdir(DIR);      # list folder under DIR

			foreach my $file (@files){
				my $testpath = qq{$errorpath/$file};
				my @list = `ls -l $testpath|grep "pdf"`;

				if (scalar @list == 0){

					print qq{ERROR:*** mistake when generate pdf into $testpath,
					Please check it using 16S.EnvironmentFactors.$p.r *** $!\n};

				}else{

					my $minsize = `ls -l $testpath|grep "pdf"|awk '{print \$5}'|sort -n |head -n 1`;
					print qq{ERROR:
					*** mistake when generate pdf into $testpath,
					Please check it using 16S.EnvironmentFactors.$p.r *** $!\n} if $minsize < 4500;   # exists & size>4000

				}
			}

			closedir DIR;

		}
	}
	# system qq{touch $tmpPath/env.finish};

}


sub Relative_Quantitation{

	### prepare  Rbsolute_Quantitation ###
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi     = qq{$metadata->{dissi}};
	my $envpath   = qq{$metadata->{env}};
	my $Rscript   = qq{$base->{Rscript_bin}};
	my $vegan     = qq{$base->{vegan_bin}};
	my $scriptdir = qq{$base->{scriptdir}/sub_script};

	### prepare output path
	my $EnvironmentFactorsPath = qq{$projectpath/Relative_Quantitation/EnvironmentFactors};
	system qq{mkdir -p $EnvironmentFactorsPath} if not -d $EnvironmentFactorsPath;

	my $tmpPath = qq{$intm_result/Relative_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $otufile     = qq{$tmpPath/subsample_otu.tax.$dissi.xls};
	my $taxonpath   = qq{$projectpath/Relative_Quantitation/Community/Community_Structure};

	# return if -e qq{$tmpPath/env.finish};

	my $envfile     = qq{$projectpath/env.xls};
	system qq{$Rscript $scriptdir/16S.EnvironmentFactors.select.r $envpath $groupfile $envfile};


	my %path;
	my @function = ("BioEnv", "MantelTest", "RDA_CCA", "Correlation");
	foreach (@function){

	    my $str = $_."Path"; 
		$path{$str} = $EnvironmentFactorsPath."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};

	}

	my @taxonfiles = `find $taxonpath -name *.Abundance.xls`;
	die "ERROR: *** There is no files named *.taxon.Abundance.xls in $taxonpath, please ckeck it (may rename files?) ***\n" if scalar @taxonfiles < 1;

	#### use functions
	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	for my $f (@function){
		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;
		if ($f eq "BioEnv"){
			# =============== BioEnv 筛选与群落相关性最大的环境因子组合 ===============
			system qq{$vegan Rscript $scriptdir/16S.EnvironmentFactors.BioEnv.r $envfile $otufile $path{"BioEnvPath"} &>> $tmpPath/env.log};
		}

		if ($f eq "MantelTest"){
			# =============== Mental Test 检验样本群落分布与环境因子整体之间的相关性 ===============
			system qq{$vegan Rscript $scriptdir/16S.EnvironmentFactors.MantelTest.r $envfile $otufile $path{"MantelTestPath"} &>> $tmpPath/env.log};
		}	

		if ($f eq "RDA_CCA"){
			foreach my $taxonfile (@taxonfiles) {
				chomp $taxonfile;
				my ($taxon) = $taxonfile =~ /\/(\w+)\.[^\/]+$/;

				# ========== RDA_CCA ===============
				my $RDA_CCAPath = qq{$path{RDA_CCAPath}/$taxon};
				system qq{mkdir -p $RDA_CCAPath} if not -d $RDA_CCAPath;
				system qq{$vegan Rscript $scriptdir/16S.EnvironmentFactors.RDA_CCA.r $envfile $taxonfile $groupfile $RDA_CCAPath &>> $tmpPath/env.log};
				system qq{$vegan Rscript $scriptdir/RDA_degree.R $RDA_CCAPath/sample.position.csv $RDA_CCAPath/environmentfactors.position.csv $RDA_CCAPath/sample_with_env_relationship.xls &>> $tmpPath/env.log};
		    	system qq{$vegan Rscript $scriptdir/RDA_degree.R $RDA_CCAPath/species.position.csv $RDA_CCAPath/environmentfactors.position.csv $RDA_CCAPath/species_with_env_relationship.xls &>> $tmpPath/env.log};
			}
		}

		if ($f eq "Correlation"){
			foreach my $taxonfile (@taxonfiles) {
				chomp $taxonfile;
				my ($taxon) = $taxonfile =~ /\/(\w+)\.[^\/]+$/;
				# ========== Correlation ===============
				my $CorrelationPath = qq{$path{CorrelationPath}/$taxon};
				system qq{mkdir -p $CorrelationPath} if not -d $CorrelationPath;
				system qq{$vegan Rscript $scriptdir/16S.EnvironmentFactors.Correlation.r $envfile $taxonfile 100 $CorrelationPath &>> $tmpPath/env.log};
			}
		}

		$pm->finish; # Terminates the child process
	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

	### check: tell you Where is the error!

	foreach my $p (keys %path){

		my $errorpath = $path{$p};
		$p =~ s/Path//;

		if ($p eq "BioEnv" or $p eq "MantelTest"){	
			my $file = `ls $errorpath`;
			chomp $file;
			my $filelinesnum = `less $errorpath/$file | wc -l`;

			print qq{ERROR:**** mistakes when generate result into $errorpath, 
			Please check it using 16S.EnvironmentFactors.$p.r *** $!\n} if $filelinesnum < 2;    # 结果行数小于2
		}

		if ($p eq "RDA_CCA" or $p eq "Correlation"){

			opendir DIR,$errorpath or die "Can not open $errorpath\n";
			my @files = grep {-d "$errorpath/$_" && ! /^\.{1,2}$/} readdir(DIR);      # list folder under DIR

			foreach my $file (@files){	
				my $testpath = qq{$errorpath/$file};		
				my @list = `ls -l $testpath|grep "pdf"`;

				if (scalar @list == 0){
					print qq{ERROR:*** mistake when generate pdf into $testpath,
					Please check it using 16S.EnvironmentFactors.$p.r *** $!\n};
					
				}else{
						
					my $minsize = `ls -l $testpath|grep "pdf"|awk '{print \$5}'|sort -n |head -n 1`;
					print qq{ERROR:
					*** mistake when generate pdf into $testpath,
					Please check it using 16S.EnvironmentFactors.$p.r *** $!\n} if $minsize < 4500;   # exists & size>4000					

				}
			}
			closedir DIR;
		}
	}
	# system qq{touch $tmpPath/env.finish};

}

1;

