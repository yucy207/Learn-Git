package environmentFactors;

use Parallel::ForkManager;
use strict;

sub run
{
	my $metadata    = shift;
	my $base        = shift;

	my @envs      = @{$metadata->{'env'}};
	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);

	for my $groupfile (@groups){

		die "ERROR: *** There is no $groupfile, please check it: $!***\n" if not -e $groupfile;
		
		my @groupname   = split /\//,$groupfile;
		my ($groupname) = $groupname[-1] =~ /(.*).xls/;
		my $intm_result = qq{$metadata->{intm_result}/$groupname};
		my $projectpath = qq{$metadata->{result}/$groupname};
		
		my $envfile;
		for my $file (@envs){
        	$envfile = $file if $file  =~ qq{$groupname\_env.xls};
        }
		
		if (-e qq{$intm_result/tmp/env.finish}){
			print qq{$groupname ENV分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/env.finish};
			system qq{rm $intm_result/tmp/env.log};
		} 

		# Forks and returns the pid for the child:
		$manager->start and next;

		### ENV分析
		environmentFactors($intm_result,$projectpath,$groupfile,$envfile,$metadata,$base);
		
		system qq{touch $intm_result/tmp/env.finish};
		print qq{$groupname ENV分析已经运行完成\n};

		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{ENV分析已经运行完成\n};
}

sub environmentFactors
{
	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $envfile 	= shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi     = qq{$metadata->{dissi}};
	my $vegan     = qq{$base->{vegan_bin}};
	my $scriptdir = qq{$base->{scriptdir}};

	### prepare output path
	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $EnvironmentFactorsPath = qq{$projectpath/EnvironmentFactors};
	system qq{mkdir -p $EnvironmentFactorsPath} if not -d $EnvironmentFactorsPath;

	#system qq{touch $tmpPath/env.finish};
	return if -e qq{$tmpPath/env.finish};
	
	#my $otufile   = qq{$intm_result/tmp/subsample_otu.tax.$dissi.xls};
	my $otufile   = -e qq{$tmpPath/subsample_otu.tax.$dissi.xls} ? qq{$tmpPath/subsample_otu.tax.$dissi.xls} : qq{$tmpPath/otu.tax.$dissi.xls};
	my $taxonpath = qq{$projectpath/Community/Community_Structure};

	die "ERROR: *** There is no $envfile, please check it: $!***\n" if not -e $envfile;
	die "ERROR: *** There is no $otufile, please check it: $!***\n" if not -e $otufile;

	system qq{cp $envfile $projectpath/env.xls};

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

}
1;
