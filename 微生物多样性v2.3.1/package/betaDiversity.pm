package betaDiversity;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata      = shift;
	my $base          = shift;

	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);

	for my $group (@groups){
		
		my $groupfile = (split /:/,$group)[0];
		my $groupbase = (split /:/,$group)[1];

		die "ERROR: *** There is no $groupfile, please check it: $!***\n" if not -e $groupfile;
		
		my @groupname     = split /\//,$groupfile;
		my ($groupname)   = $groupname[-1] =~ /(.*).xls/;
		my $intm_result   = qq{$metadata->{intm_result}/$groupname};
		my $projectpath   = qq{$metadata->{result}/$groupname};
		
		if (-e qq{$intm_result/tmp/beta.finish}){
			print qq{$groupname Beta分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/beta.finish};
			system qq{rm $intm_result/tmp/beta.log};
		}

		# Forks and returns the pid for the child:
		$manager->start and next;
		
		### Beta分析
		betaDiversity($intm_result,$projectpath,$groupfile,$groupbase,$metadata,$base);
		
		system qq{touch $intm_result/tmp/beta.finish};
		print qq{$groupname Beta分析已经运行完成\n};

		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{Beta分析已经运行完成\n};
}

sub betaDiversity
{
	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $groupbase   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi     = qq{$metadata->{dissi}};
	my $lefse     = qq{$base->{lefse_bin}};
	my $Rscript   = qq{$base->{Rscript_bin}};
	my $Rscript_3_5_3  = qq{$base->{Rscript_3_5_3}};
	my $scriptdir = qq{$base->{scriptdir}};

	### prepare output path
	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $BetaDiversity = qq{$projectpath/BetaDiversity};
	system qq{mkdir -p $BetaDiversity} if not -d $BetaDiversity;
	
	my $otufile       = -e qq{$tmpPath/subsample_otu.tax.$dissi.xls} ? qq{$tmpPath/subsample_otu.tax.$dissi.xls} : qq{$tmpPath/otu.tax.$dissi.xls};
	my $trefile       = -e qq{$tmpPath/subsample_otu.repseq.fasta.tre} ? qq{$tmpPath/subsample_otu.repseq.fasta.tre} : qq{$tmpPath/otu.repseq.fasta.tre};
	my $otustatRAfile = -e qq{$tmpPath/subsample_otu.tax.$dissi.stat.RA.xls} ? qq{$tmpPath/subsample_otu.tax.$dissi.stat.RA.xls} : qq{$tmpPath/otu.tax.$dissi.stat.RA.xls};

	die "ERROR: *** There is no $otufile, please check it: $!***\n" if not -e $otufile;
	die "ERROR: *** There is no $trefile, please check it: $!***\n" if not -e $trefile;
	die "ERROR: *** There is no $otustatRAfile, please check it: $!***\n" if not -e $otustatRAfile;

	

	#### use functions
	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	###差异分析
	my @samplegroup = samplegroup($groupfile,$groupbase,$intm_result);
	
	for my $samplegroup (@samplegroup){
		
		my $path = "";
		
		($path) = (split /\//,$samplegroup)[-1] =~ /(.*).xls/ if ((split /\//,$samplegroup)[-1] =~ /_vs_/);
		
		my %path;
		my @function = ("Venn", "SamplesTree", "PCA", "PLS_DA", "PCoA", "NMDS", "ADONIS", "Lefse", "BetaNTI");
		foreach (@function){

			my $str = $_."Path"; 
			$path{$str} = $BetaDiversity."/".$_;
			#system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
		}
		my @function = ("SamplesTree", "Lefse", "BetaNTI");
		foreach (@function){

			my $str = $_."Path"; 
			$path{$str} = $BetaDiversity."/".$_."/".$path;
			#system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
		}
		
		
		my $num = `less $samplegroup | cut -f 2 | sort -u | wc -l`;
		my $num_sample =  `less $samplegroup | cut -f 1 | wc -l`;

		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;

		# =============== Venn ===============
		system qq{mkdir -p $path{"VennPath"}} if not -d $path{"VennPath"};
		system qq{$Rscript $scriptdir/16S.BetaDiversity.Venn.r $otufile $samplegroup $path{"VennPath"} $scriptdir &>> $tmpPath/beta.log};
	
		# # =============== SamplesTree ===============
		system qq{mkdir -p $path{"SamplesTreePath"}} if not -d $path{"SamplesTreePath"};
		system qq{$Rscript $scriptdir/16S.BetaDiversity.SamplesTree.r $otufile $samplegroup $path{"SamplesTreePath"} &>> $tmpPath/beta.log};

	
		# # =============== PCA ===============
		system qq{mkdir -p $path{"PCAPath"}} if not -d $path{"PCAPath"};
		system qq{$Rscript $scriptdir/16S.BetaDiversity.PCA.r $otufile $samplegroup $path{"PCAPath"} &>> $tmpPath/beta.log};

		# =============== PLS_DA ===============
		system qq{mkdir -p $path{"PLS_DAPath"}} if not -d $path{"PLS_DAPath"};
		system qq{$Rscript_3_5_3 $scriptdir/16S.BetaDiversity.PLS_DA.r $otufile $samplegroup $path{"PLS_DAPath"} &>> $tmpPath/beta.log};

	
		# # =============== PCoA ===============
		system qq{mkdir -p $path{"PCoAPath"}} if not -d $path{"PCoAPath"};
		system qq{$Rscript $scriptdir/16S.BetaDiversity.PCoA.r $otufile $trefile $samplegroup $path{"PCoAPath"} &>> $tmpPath/beta.log};

	
		# # =============== NMDS ===============
		system qq{mkdir -p $path{"NMDSPath"}} if not -d $path{"NMDSPath"};
		system qq{$Rscript $scriptdir/16S.BetaDiversity.NMDS.r $otufile $trefile $samplegroup $path{"NMDSPath"} &>> $tmpPath/beta.log};
	
		if ($num_sample > $num and $num >= 2){
			# =============== ADONIS ===============
			system qq{mkdir -p $path{"ADONISPath"}} if not -d $path{"ADONISPath"};
			# system qq{$Rscript $scriptdir/16S.BetaDiversity.PERMANOVA.r $otufile $samplegroup $path{"ADONISPath"} &>> $tmpPath/beta.log};
			system qq{$Rscript $scriptdir/16S.BetaDiversity.ADONIS.new.r $otufile $trefile $samplegroup $path{"ADONISPath"} &>> $tmpPath/beta.log};
		
			# =============== Lefse ===============
			system qq{mkdir -p $path{"LefsePath"}} if not -d $path{"LefsePath"};
			#system qq{perl $scriptdir/16S.BetaDiversity.Lefse.pl $otustatRAfile $samplegroup $path{"LefsePath"}};
			############ lefse
			my $lefsepath = qq{$path{"LefsePath"}};
			
			system qq{$Rscript $scriptdir/16S.BetaDiversity.select_stat.RA.r $otustatRAfile $samplegroup $lefsepath};
			
			statlefse("$lefsepath/stat.RA.xls", $samplegroup, $lefsepath);

			system qq{$lefse format_input.py $lefsepath/taxon.stat.for.lefse.xls $lefsepath/taxon.stat.format.xls -c 1 -o 1000000 &>> $tmpPath/beta.log};
			system qq{$lefse run_lefse.py $lefsepath/taxon.stat.format.xls $lefsepath/taxon.stat.res.xls &>> $tmpPath/beta.log};

			statres($lefsepath);
			system qq{$Rscript $scriptdir/lefse4plot.r $lefsepath/taxon.stat.res1.xls 50 $lefsepath &>> $tmpPath/beta.log};
			system qq{$lefse plot_cladogram.py $lefsepath/taxon.stat.res4plot.xls $lefsepath/taxon.stat.cladogram.pdf --format pdf --dpi 1000 --title_font_size 8 --label_font_size 3 --class_legend_font_size 5 --right_space_prop 0.25 --left_space_prop 0.1 &>> $tmpPath/beta.log};
			
			statres2($lefsepath);
			system qq{$Rscript $scriptdir/lefse4plot.r $lefsepath/taxon.stat.res2.xls 10 $lefsepath &>> $tmpPath/beta.log};
			system qq{$lefse plot_res.py $lefsepath/taxon.stat.res4plot.xls $lefsepath/taxon.stat.res.pdf --width 10 --left_space 0.35 --format pdf --dpi 1000 &>> $tmpPath/beta.log};

			#system qq{rm -rf $lefsepath/taxon.stat.for.lefse.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res4plot.xls};
			system qq{rm -rf $lefsepath/taxon.stat.format.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res.xls};
			system qq{rm -rf $lefsepath/taxon.stat.res1.xls};

			
		
			# # =============== BetaNTI ===============
			system qq{mkdir -p $path{"BetaNTIPath"}} if not -d $path{"BetaNTIPath"};
			system qq{$Rscript $scriptdir/16S.BetaDiversity.betaNTI.r $otufile $samplegroup $trefile $path{"BetaNTIPath"} &>> $tmpPath/beta.log};
			#system qq{perl $scriptdir/16S.BetaDiversity.betaNTI.pl $otufile $groupfile $trefile $path{"BetaNTIPath"} &>> $tmpPath/beta.log};
		}
		
		$pm->finish; # Terminates the child process
	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

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
