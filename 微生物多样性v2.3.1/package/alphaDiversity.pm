package alphaDiversity;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata    = shift;
	my $base        = shift;

	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);

	for my $group (@groups){
		
		my $groupfile = (split /:/,$group)[0];
		my $groupbase = (split /:/,$group)[1];

		die "ERROR: There is no $groupfile, Please Check It! $!\n" if not -e $groupfile;
		
		my @groupname = split /\//,$groupfile;
		my ($groupname) = $groupname[-1] =~ /(.*).xls/;
		my $intm_result = qq{$metadata->{intm_result}/$groupname};
		my $projectpath = qq{$metadata->{result}/$groupname};
		
		if (-e qq{$intm_result/tmp/alpha.finish}){
		
			print qq{$groupname Alpha分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/alpha.finish};
			system qq{rm $intm_result/tmp/alpha.log};
		}

		# Forks and returns the pid for the child:
		$manager->start and next;
		
		### Alpha分析	
		alphaDiversity($intm_result,$projectpath,$groupfile,$groupbase,$metadata,$base);

		system qq{touch $intm_result/tmp/alpha.finish};
		print qq{$groupname Alpha分析已经运行完成\n};

		$manager->finish;
	}
	
	$manager->wait_all_children;
	print qq{Alpha分析已经运行完成\n};
}


sub alphaDiversity
{

	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $groupbase   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $mothur      = qq{$base->{alpha_mothur_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}};

	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	return if -e qq{$tmpPath/alpha.finish};


	my $AlphaDiversity = qq{$projectpath/AlphaDiversity};
	system "rm -rf $AlphaDiversity" if -d $AlphaDiversity;
	system qq{mkdir -p $AlphaDiversity};

	#system qq{touch $tmpPath/alpha.finish};

	
	#my $otufile = qq{$tmpPath/subsample_otu.tax.$dissi.xls};
	my $otufile  = -e qq{$tmpPath/subsample_otu.tax.$dissi.xls} ? qq{$tmpPath/subsample_otu.tax.$dissi.xls} : qq{$tmpPath/otu.tax.$dissi.xls};
	die "ERROR: There is no $otufile, Please Check It! $!\n" if not -e $otufile;
	
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
	my $num = `less $groupfile | cut -f 2 | sort -u | wc -l`;
	my $num_sample =  `less $groupfile | cut -f 1 | wc -l`;

	system qq{$Rscript $scriptdir/otu2shared.r $otufile $groupfile $dissi $AlphaDiversity &>> $tmpPath/alpha.log};
	my $sharedfile = qq{$AlphaDiversity/sample.shared};

	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);
	
	###差异分析
	my @samplegroup = samplegroup($groupfile,$groupbase,$intm_result);
	
	for my $samplegroup (@samplegroup){
		
		my $path = "";
		
		($path) = (split /\//,$samplegroup)[-1] =~ /(.*).xls/ if ((split /\//,$samplegroup)[-1] =~ /_vs_/);
		system qq{mkdir -p $path{"DiversityIndexPath"}/$path} if not -d qq{$path{"DiversityIndexPath"}/$path};
		system qq{$Rscript $scriptdir/16S.AlphaDiversity.DiversityIndex.r $otufile $samplegroup $path{"DiversityIndexPath"}/$path &>> $tmpPath/alpha.log};
		
	}
	
	for my $f (@function){

		# Forks and returns the pid for the child:
		# my $pid = $pm->start and next;

		#if ($f eq "DiversityIndex"){
		#	# =============== DiversityIndex ===============
		#	system qq{$Rscript $scriptdir/16S.AlphaDiversity.DiversityIndex.r $otufile $groupfile $path{"DiversityIndexPath"} &>> $tmpPath/alpha.log};
		#}
		
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

			if ($num_sample > $num and $num >= 2){
				my $ShannonPath_groupPath = qq{$path{"ShannonPath"}/group};
				system qq{mkdir -p  $ShannonPath_groupPath} if not -d qq{$ShannonPath_groupPath};
				system qq{rm $ShannonPath_groupPath/*.xls} if -e qq{$ShannonPath_groupPath/*.xls};
				system qq{$Rscript $scriptdir/16S.otu.group.r $otufile $groupfile $ShannonPath_groupPath &>> $tmpPath/alpha.log};
				system qq{$Rscript $scriptdir/otu2shared.r "$ShannonPath_groupPath/otu.group.xls" "$ShannonPath_groupPath/sample.group" $dissi $ShannonPath_groupPath &>> $tmpPath/alpha.log};
				my $sharedfile2 = qq{$ShannonPath_groupPath/sample.shared};
				system qq{$mothur bash -c 'mothur "#rarefaction.single(shared=$sharedfile2,label=$dissi,calc=shannon,groupmode=f,processors=20)";chown -R $uid:$gid $ShannonPath_groupPath' &>> $tmpPath/alpha.log};
				system qq{$Rscript $scriptdir/16S.AlphaDiversity.ShannonCurves.r $ShannonPath_groupPath 40 $ShannonPath_groupPath &>> $tmpPath/alpha.log};
				system qq{mv $ShannonPath_groupPath/Samples_Shannon_Curves.pdf $ShannonPath_groupPath/Groups_Shannon_Curves.pdf};

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

			if ($num_sample > $num and $num >= 2){
				my $Rarefaction_groupPath = qq{$path{"RarefactionPath"}/group};
				system qq{mkdir -p  $Rarefaction_groupPath} if not -d qq{$Rarefaction_groupPath};
				system qq{rm $Rarefaction_groupPath/*.xls} if -e qq{$Rarefaction_groupPath/*.xls};
				system qq{$Rscript $scriptdir/16S.otu.group.r $otufile $groupfile $Rarefaction_groupPath &>> $tmpPath/alpha.log};
				system qq{$Rscript $scriptdir/otu2shared.r "$Rarefaction_groupPath/otu.group.xls" "$Rarefaction_groupPath/sample.group" $dissi $Rarefaction_groupPath &>> $tmpPath/alpha.log};
				my $sharedfile2 = qq{$Rarefaction_groupPath/sample.shared};
				system qq{$mothur bash -c 'mothur "#rarefaction.single(shared=$sharedfile2,label=$dissi,groupmode=f,processors=20)";chown -R $uid:$gid $Rarefaction_groupPath' &>> $tmpPath/alpha.log};
				system qq{$Rscript $scriptdir/16S.AlphaDiversity.RarefactionCurves.r $Rarefaction_groupPath 40 $Rarefaction_groupPath &>> $tmpPath/alpha.log};
				system qq{mv $Rarefaction_groupPath/Samples_Rarefaction_Curves.pdf $Rarefaction_groupPath/Groups_Rarefaction_Curves.pdf};

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
		
		# $pm->finish; # Terminates the child process

	}

	# $pm->wait_all_children;  # blocks until all forked processes have exited

	system qq{rm -rf $sharedfile};

	### check: tell you Where is the error!
	foreach my $key (keys %path){

		my $errorpath = qq{$path{$key}};
		my $minsize = `ls -lR $errorpath|grep "pdf"|awk '{print \$5}'| sort -n |head -n 1`;
		$key =~ s/Path//;

		print qq{WARINING/ERROR:**** mistakes when generate pdf file into $errorpath, 
		Please check it using 16S.AlphaDiversity.$key.r *** $!\n} if $minsize < 4500;   # exists & size>5000
	}

}

sub samplegroup	{
	
	my $groupfile = shift;
	my $groupbase = shift;
	my $grouppath = shift;
	
	my @samplegroup;
	
	push @samplegroup, $groupfile;
	
	if (defined $groupbase){
		
		my @groupname = split /;/,$groupbase;
		
		for my $groupname (@groupname){

			my @group = split /,/, $groupname;
			
			my $group = join "_vs_", @group; 
			
			my $samplegroup = qq{$grouppath/$group\.xls};
			
			open FILE ,$groupfile;
			open OUT,">$samplegroup";
			while (<FILE>) {
				
				s/[\r\n]//g;
				
				my @arr = split /\t/, $_;
				
				print OUT qq{$_\n} if ($arr[1] ~~ @group);

			}
			close FILE;
			close OUT;
			push @samplegroup, $samplegroup;
			
		}
		
	}
	# else{
		
		# push @samplegroup, $groupfile;
	# }
	
	return @samplegroup;
}

1;
