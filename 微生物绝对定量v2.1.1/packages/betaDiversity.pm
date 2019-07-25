package betaDiversity;

use strict;
use Parallel::ForkManager;
use List::MoreUtils qw/uniq/;

sub run
{
	my $metadata      = shift;
	my $base          = shift;
	my $quantitation  = shift;

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
			Absolute_Quantitation($projectpath,$intm_result,$groupfile,$select,$metadata,$base);
		}

		if (uc($quantitation) eq "R" || uc($quantitation) eq "B"){
			Relative_Quantitation($projectpath,$intm_result,$groupfile,$select,$metadata,$base);
		}

		# system qq{touch $tmpPath/beta.finish};
		print qq{$groupname beta分析完成\n};
		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{beta分析已经运行完成\n};
}


sub Absolute_Quantitation{
	### prepare Absolute_Quantitation ###
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi           = qq{$metadata->{dissi}};
	my $lefse           = qq{$base->{lefse_bin}};
	my $Rscript         = qq{$base->{Rscript_bin}};
	my $Rscript_3_5_3   = qq{$base->{Rscript_3_5_3}};
	my $scriptdir       = qq{$base->{scriptdir}/sub_script};

	### prepare Sample.groups
	precheck($groupfile);

	my $tmpPath = qq{$intm_result/Absolute_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $BetaDiversity = qq{$projectpath/Absolute_Quantitation/BetaDiversity};
	system qq{mkdir -p $BetaDiversity} if not -d $BetaDiversity;

	my $otufile       = qq{$tmpPath/otu.tax.$dissi.xls};
	# my $otufile1      = qq{$tmpPath/copies_DNA_level_otu.tax.$dissi.xls};
	my $otustatfile = qq{$tmpPath/otu.tax.0.03.stat.for.lefse.xls};
	my $trefile       = qq{$tmpPath/otu.repseq.fasta.tre};

	# return if -e qq{$tmpPath/IRbeta.finish};
	
	my %path;

	my @function = ("SamplesTree", "PCA", "PCoA", "NMDS", "PLS_DA", "ADONIS", "Venn", "Lefse", "BetaNTI");
	foreach (@function){

	    my $str = $_."Path"; 
		$path{$str} = $BetaDiversity."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}


	#### use functions
	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	for my $f (@function){

		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;

		if ($f eq "SamplesTree"){
			# =============== SamplesTree ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.SamplesTree.r $otufile $groupfile $path{"SamplesTreePath"} &>> $tmpPath/IRbeta.log};

		}
		
		if ($f eq "PCA"){
			# =============== PCA ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.PCA.r $otufile $groupfile "$select" $path{"PCAPath"} &>> $tmpPath/IRbeta.log};
		}
		
		if ($f eq "PCoA"){
			# =============== PCoA ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.PCoA.r $otufile $trefile $groupfile "$select" $path{"PCoAPath"} &>> $tmpPath/IRbeta.log};

		}
		
		if ($f eq "NMDS"){
			# =============== NMDS ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.NMDS.r $otufile $trefile $groupfile "$select" $path{"NMDSPath"} &>> $tmpPath/IRbeta.log};
		}

		if ($f eq "PLS_DA"){
			# =============== PLS_DA ===============
			system qq{$Rscript_3_5_3 $scriptdir/16S.BetaDiversity.PLS_DA.r $otufile $groupfile "$select" $path{"PLS_DAPath"} &>> $tmpPath/IRbeta.log};
		}


		
		if ($f eq "ADONIS"){
			# =============== PERMANOVA ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.ADONIS.new.r $otufile $trefile $groupfile "$select" $path{"ADONISPath"} &>> $tmpPath/IRbeta.log};
		}

		
		# OTU 水平
		if ($f eq "Venn"){
			# =============== Venn ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.Venn.r $otufile $groupfile "$select" $path{"VennPath"} $scriptdir &>> $tmpPath/IRbeta.log};

		}


		if ($f eq "Lefse"){
			# =============== Lefse ===============
			my $lefsepath = qq{$path{"LefsePath"}};

		
			statlefse($otustatfile, $groupfile, $lefsepath);

			system qq{$lefse format_input.py $lefsepath/taxon.stat.for.lefse.xls $lefsepath/taxon.stat.format.xls -c 1 -o 1000000 &>> $tmpPath/IRbeta.log};
			system qq{$lefse run_lefse.py $lefsepath/taxon.stat.format.xls $lefsepath/taxon.stat.res.xls &>> $tmpPath/IRbeta.log};

			statres($lefsepath);
			system qq{$Rscript $scriptdir/lefse4plot.r $lefsepath/taxon.stat.res1.xls 50 $lefsepath &>> $tmpPath/IRbeta.log};
			system qq{$lefse plot_cladogram.py $lefsepath/50.taxon.stat.res4plot.xls $lefsepath/taxon.stat.cladogram.pdf --format pdf --dpi 1000 --title_font_size 8 --label_font_size 3 --class_legend_font_size 5 --right_space_prop 0.25 --left_space_prop 0.1 &>> $tmpPath/IRbeta.log};
			
			statres2($lefsepath);
			system qq{$Rscript $scriptdir/lefse4plot.r $lefsepath/taxon.stat.res2.xls 10 $lefsepath &>> $tmpPath/IRbeta.log};
			system qq{$lefse plot_res.py $lefsepath/10.taxon.stat.res4plot.xls $lefsepath/taxon.stat.res.pdf --width 10 --left_space 0.35 --format pdf --dpi 1000 &>> $tmpPath/IRbeta.log};

			#system qq{rm -rf $lefsepath/taxon.stat.for.lefse.xls};
			system qq{rm -rf $lefsepath/50.taxon.stat.res4plot.xls};
			system qq{rm -rf $lefsepath/10.taxon.stat.res4plot.xls};
			system qq{rm -rf $lefsepath/taxon.stat.format.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res1.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res2.xls};


			# my $lefseresult = `cat $lefsepath/taxon.lefse.result.xls| wc -l`;
			# if ($lefseresult == 1) {
			# 	system qq{rm -rf $lefsepath};
			# 	print "$groupfile的Lefse分析无差异物种，结果已删除\n";
			# }

		}	


		if ($f eq "BetaNTI"){
			# =============== BetaNTI ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.betaNTI.r $otufile $groupfile $trefile $path{"BetaNTIPath"} &>> $tmpPath/IRbeta.log};

		}
		
			

		$pm->finish; # Terminates the child process
		
	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

	### check: tell you Where is the error!
	foreach my $p (keys %path){

		my $errorpath = $path{$p};
		$p =~ s/Path//;
		my $errorinfo = qq{ERROR:**** mistakes when generate pdf file into $errorpath, 
		Please check it using 16S.BetaDiversity.$p.r *** $!\n};

		if ($p eq "SamplesTree"){
			my @list = `ls -l $errorpath|grep "pdf"`;
			if (scalar @list == 0){
				print $errorinfo;
			}else{
				my $minsize = `ls -l $errorpath|grep "pdf"|awk '{print \$5}'| sort -n |head -n 1`;
				print $errorinfo if $minsize < 4000;   # exists & size>4000
			}
		}
		
		if ($p eq "PCA" or $p eq "PCoA" or $p eq "NMDS" or $p eq "Lefse"){
			my @list = `ls -lR $errorpath|grep "pdf"`;
			if (scalar @list == 0){
				print $errorinfo;
			}
		}
	}

	# system qq{touch $tmpPath/IRbeta.finish};
}


sub Relative_Quantitation{
	### prepare  Rbsolute_Quantitation ###
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi           = qq{$metadata->{dissi}};
	my $lefse           = qq{$base->{lefse_bin}};
	my $Rscript         = qq{$base->{Rscript_bin}};
	my $Rscript_3_5_3   = qq{$base->{Rscript_3_5_3}};
	my $scriptdir       = qq{$base->{scriptdir}/sub_script};

	### prepare Sample.groups
	precheck($groupfile);

	my $tmpPath = qq{$intm_result/Relative_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $BetaDiversity = qq{$projectpath/Relative_Quantitation/BetaDiversity};
	system qq{mkdir -p $BetaDiversity} if not -d $BetaDiversity;

	my $otufile       = qq{$tmpPath/subsample_otu.tax.$dissi.xls};
	my $otustatRAfile = qq{$tmpPath/subsample_otu.tax.$dissi.stat.RA.xls};
	my $trefile       = qq{$tmpPath/subsample_otu.repseq.fasta.tre};
	# my $taxonpath   = qq{$projectpath/Community/Community_Structure};


	my %path;
	my @function = ("SamplesTree", "PCA", "PCoA", "NMDS", "PLS_DA", "ADONIS", "Lefse", "Venn", "BetaNTI");
	foreach (@function){

	    my $str = $_."Path"; 
		$path{$str} = $BetaDiversity."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}

	#### use functions
	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	for my $f (@function){

		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;

		if ($f eq "SamplesTree"){
			# =============== SamplesTree ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.SamplesTree.r $otufile $groupfile $path{"SamplesTreePath"} &>> $tmpPath/beta.log};

		}
		
		if ($f eq "PCA"){
			# =============== PCA ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.PCA.r $otufile $groupfile "$select" $path{"PCAPath"} &>> $tmpPath/beta.log};
		}
		
		if ($f eq "PCoA"){
			# =============== PCoA ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.PCoA.r $otufile $trefile $groupfile "$select" $path{"PCoAPath"} &>> $tmpPath/beta.log};

		}
		
		if ($f eq "NMDS"){
			# =============== NMDS ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.NMDS.r $otufile $trefile $groupfile "$select" $path{"NMDSPath"} &>> $tmpPath/beta.log};
		}
		

		if ($f eq "PLS_DA"){
			# =============== PLS_DA ===============
			system qq{$Rscript_3_5_3 $scriptdir/16S.BetaDiversity.PLS_DA.r $otufile $groupfile "$select" $path{"PLS_DAPath"} &>> $tmpPath/beta.log};
		}


		if ($f eq "ADONIS"){
			# =============== PERMANOVA ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.ADONIS.new.r $otufile $trefile $groupfile "$select" $path{"ADONISPath"} &>> $tmpPath/IRbeta.log};
		}
		
		if ($f eq "Lefse"){
			# =============== Lefse ===============
			#system qq{perl $scriptdir/16S.BetaDiversity.Lefse.pl $otustatRAfile $groupfile $path{"LefsePath"}};
			############ lefse
			my $lefsepath = qq{$path{"LefsePath"}};

			statlefse($otustatRAfile, $groupfile, $lefsepath);

			system qq{$lefse format_input.py $lefsepath/taxon.stat.for.lefse.xls $lefsepath/taxon.stat.format.xls -c 1 -o 1000000 &>> $tmpPath/beta.log};
			system qq{$lefse run_lefse.py $lefsepath/taxon.stat.format.xls $lefsepath/taxon.stat.res.xls &>> $tmpPath/beta.log};

			statres($lefsepath);
			system qq{$Rscript $scriptdir/lefse4plot.r $lefsepath/taxon.stat.res1.xls 50 $lefsepath &>> $tmpPath/beta.log};
			system qq{$lefse plot_cladogram.py $lefsepath/50.taxon.stat.res4plot.xls $lefsepath/taxon.stat.cladogram.pdf --format pdf --dpi 1000 --title_font_size 8 --label_font_size 3 --class_legend_font_size 5 --right_space_prop 0.25 --left_space_prop 0.1 &>> $tmpPath/beta.log};
			
			statres2($lefsepath);
			system qq{$Rscript $scriptdir/lefse4plot.r $lefsepath/taxon.stat.res2.xls 10 $lefsepath &>> $tmpPath/beta.log};
			system qq{$lefse plot_res.py $lefsepath/10.taxon.stat.res4plot.xls $lefsepath/taxon.stat.res.pdf --width 10 --left_space 0.35 --format pdf --dpi 1000 &>> $tmpPath/beta.log};

			#system qq{rm -rf $lefsepath/taxon.stat.for.lefse.xls};
			system qq{rm -rf $lefsepath/10.taxon.stat.res4plot.xls};
			system qq{rm -rf $lefsepath/50.taxon.stat.res4plot.xls};
			system qq{rm -rf $lefsepath/taxon.stat.format.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res1.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res2.xls};
		}

		# OTU 水平
		if ($f eq "Venn"){
			# =============== Venn ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.Venn.r $otufile $groupfile "$select" $path{"VennPath"} $scriptdir &>> $tmpPath/beta.log};

		}

		if ($f eq "BetaNTI"){
			# =============== BetaNTI ===============
			system qq{$Rscript $scriptdir/16S.BetaDiversity.betaNTI.r $otufile $groupfile $trefile $path{"BetaNTIPath"} &>> $tmpPath/beta.log};

		}
		
		
		$pm->finish; # Terminates the child process
		
	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

	### check: tell you Where is the error!
	foreach my $p (keys %path){

		my $errorpath = $path{$p};
		$p =~ s/Path//;
		my $errorinfo = qq{ERROR:**** mistakes when generate pdf file into $errorpath, 
		Please check it using 16S.BetaDiversity.$p.r *** $!\n};

		if ($p eq "SamplesTree"){
			my @list = `ls -l $errorpath|grep "pdf"`;
			if (scalar @list == 0){
				print $errorinfo;
			}else{
				my $minsize = `ls -l $errorpath|grep "pdf"|awk '{print \$5}'| sort -n |head -n 1`;
				print $errorinfo if $minsize < 4000;   # exists & size>4000
			}
		}
		
		if ($p eq "PCA" or $p eq "PCoA" or $p eq "NMDS" or $p eq "Lefse"){
			my @list = `ls -lR $errorpath|grep "pdf"`;
			if (scalar @list == 0){
				print $errorinfo;
			}
		}
	}

}


sub precheck
{
	my $groupfile = shift;

	open GROUP,$groupfile or die "Can not open $groupfile: $!\n";      # 读取样本分组信息
	my %GroupInfo;
	while(<GROUP>){

		chomp;
		next if /^\#|^\s+$/;    # 跳过注释行和空行
		my ($sample, $Group) = split /\s+/;
		$GroupInfo{$sample} = $Group if not exists $GroupInfo{$sample};

	}
	close GROUP;
	my @uniq_Groups = uniq values %GroupInfo;       # 去重
	die "ERROR: *** Less than two groups!: $!\n" if scalar @uniq_Groups < 2 ;    # 组数小于2
}

sub statlefse
{
	my $otustatRAfile = shift;	
	my $groupfile     = shift;
	my $lefsepath    = shift;

	open GROUP, $groupfile or die "Can not open $groupfile, Please Check It!\n";
	open STAT, $otustatRAfile or die "Can not open $otustatRAfile, Please Check It!\n";
	open FILE, qq{>$lefsepath/taxon.stat.for.lefse.xls} or die "Can not write to file!\n";

	############ prepare data for lefse
	my %sampleinfo;
	while(<GROUP>){

	    s/\s+$//g;
		chomp;
		next if /^\s+$/;
		
		my @data = split /\s+/;
				
		die "Sample name duplication, please check sample.groups!\n" if exists $sampleinfo{$data[0]};
		
		$sampleinfo{$data[0]} = $data[1];

	}
	close GROUP;

	while(<STAT>){

		chomp;
		s/\s+$//;
		my @names = split /\t/;
		
		if (/Taxlevel/){
		
			for my $i (1..$#names){
			
				if (exists $sampleinfo{$names[$i]}){
				
					$names[$i] = $sampleinfo{$names[$i]};       # change sample name to groups name 
				
				}else{
				
					die "Sample names in $groupfile and $otustatRAfile do not match, Check It !!!\n";
				
				}
			
			}
			
		}

		my $out = join "\t",@names[0..$#names];    # do not export Size column
	 	print FILE qq{$out\n};

	}

	close STAT;
	close FILE;
}

sub statres
{
	my $lefsepath = shift;

	open(STAT,"$lefsepath/taxon.stat.res.xls") || die "cannot open taxon.stat.res.xls\n";
	open(OUT,">$lefsepath/taxon.stat.res1.xls") || die "cannot make taxon.stat.res1.xls\n";
	while(my $line = <STAT>){
		my $taxon = (split/\t/,$line)[0];
		my @arr = split/\./,$taxon;
		next if $arr[$#arr] =~/Unassigned/;
		print OUT "$line";
	}
	close STAT;
	close OUT;
}

sub statres2
{
	my $lefsepath = shift;

	open(STAT,"$lefsepath/taxon.stat.res.xls") || die "cannot open taxon.stat.res.xls\n";
	open(OUT,">$lefsepath/taxon.stat.res2.xls") || die "cannot make taxon.stat.res2.xls\n";
	while(my $line = <STAT>){
		my $taxon = (split/\t/,$line)[0];
		my @arr = split/\./,$taxon;
		next if $arr[$#arr] =~/Unassigned/;
		next if $arr[$#arr] =~/s__/;
		print OUT "$line";
	}
	close STAT;
	close OUT;
}

1;
