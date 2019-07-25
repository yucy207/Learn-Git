package community;

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

		if (-e qq{$intm_result/tmp/community.finish}){
			print qq{$groupname Community分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/community.finish};
			system qq{rm $intm_result/tmp/community.log};
		}

		# Forks and returns the pid for the child:
		$manager->start and next;

		### Community分析
		community($intm_result,$projectpath,$groupfile,$groupbase,$metadata,$base);

		system qq{touch $intm_result/tmp/community.finish};
		print qq{$groupname Community分析已经运行完成\n};

		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{Community分析已经运行完成\n};
}


sub community
{
	
	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $groupbase	= shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $kronatools  = qq{$base->{kronatools_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}};
	
	
	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $Community = qq{$projectpath/Community};
	system qq{mkdir -p $Community} if not -d $Community;

	#system qq{touch $tmpPath/community.finish};
	return if -e qq{$tmpPath/community.finish};
	
	#my $otufile     = qq{$tmpPath/subsample_otu.tax.$dissi.xls};
	#my $otustatfile = qq{$tmpPath/subsample_otu.tax.$dissi.stat.xls};
	my $otufile     = -e qq{$tmpPath/subsample_otu.tax.$dissi.xls} ? qq{$tmpPath/subsample_otu.tax.$dissi.xls} : qq{$tmpPath/otu.tax.$dissi.xls};
	my $otustatfile = -e qq{$tmpPath/subsample_otu.tax.$dissi.stat.xls} ? qq{$tmpPath/subsample_otu.tax.$dissi.stat.xls} : qq{$tmpPath/otu.tax.$dissi.stat.xls};

	die "ERROR: There is no $otufile, Please Check It! $!\n" if not -e $otufile;
	die "ERROR: There is no $groupfile, Please Check It! $!\n" if not -e $groupfile;

	# prepare output path
	my %path;
	foreach("TaxonTree", "Community_Structure", "Pieplot", "Barplot", "Heatmap", "Treebar", "Bubble", "KronaPlot"){

		my $str = $_."Path";
		$path{$str} = $Community."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}

	# =============== TaxonTree =============== 
	system qq{$Rscript $scriptdir/16S.Community.TaxonTree.r $otustatfile $path{"TaxonTreePath"} &>> $tmpPath/community.log};

	# =============== OTUSplit =============== 
	system qq{$Rscript $scriptdir/16S.CommunityStructure.OTUSplit.r $otufile $path{"Community_StructurePath"} &>> $tmpPath/community.log};

	opendir Dir, $path{"Community_StructurePath"} or die "Can not cd to $path{'Community_StructurePath'}: $!\n";
	my @fileslist = grep { /\.Abundance.xls/ } readdir(Dir);
	die qq{ERROR: Can not generate *.Abundance.xls into $path{Community_StructurePath}, Please check it using 16S.CommunityStructure.OTUSplit.r $!\n} if scalar @fileslist == 0;      # 注意格式 die...if...
	close Dir;

	my $num = `less $groupfile | cut -f 2 | sort -u | wc -l`;
	my $num_sample =  `less $groupfile | cut -f 1 | wc -l`;

	my $MAX_PROCESSES = 10;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	foreach my $file (@fileslist) {

		my $dir = (split /\./, $file)[0];      # use taxon level as dir
		my $taxondir = qq{$path{Community_StructurePath}/$dir};
		my $taxondir2 = qq{$path{PieplotPath}/$dir};
		
		system qq{mkdir -p $taxondir} if not -d qq{$taxondir};
		system qq{mkdir -p $taxondir2} if not -d qq{$taxondir2};
		
		system qq{mv $path{Community_StructurePath}/$dir.* $taxondir};
		
		my $taxonfile = qq{$taxondir/$file};
		
		# =============== CommunityStructure.Barplot =============== 
		system qq{$Rscript $scriptdir/16S.CommunityStructure.Barplot.r $taxonfile $taxondir &>> $tmpPath/community.log};
		
		# =============== Pieplot =============== 
		system qq{$Rscript $scriptdir/16S.Community.Pieplot.r $taxonfile $groupfile $taxondir2 &>> $tmpPath/community.log};

		# =============== Barplot =============== 
		system qq{$Rscript $scriptdir/16S.Community.Barplot.r $taxonfile $groupfile $path{"BarplotPath"} &>> $tmpPath/community.log};
		
		system qq{$Rscript $scriptdir/16S.Community.Barplot.group.r $taxonfile $groupfile $path{"BarplotPath"} &>> $tmpPath/community.log} if ($num_sample > $num and $num != 1);
	
		# =============== Heatmap =============== 
		system qq{$Rscript $scriptdir/16S.Community.Heatmap.r $taxonfile 100 $path{"HeatmapPath"} &>> $tmpPath/community.log};
		
		# =============== TreeBar =============== 
		system qq{$Rscript $scriptdir/16S.Community.TreeBar.r $taxonfile $groupfile $path{"TreebarPath"} &>> $tmpPath/community.log};

		# =============== Bubble =============== 
		system qq{$Rscript $scriptdir/16S.Community.Bubble.r $taxonfile $groupfile 30 $path{"BubblePath"} &>> $tmpPath/community.log};
		
		###差异分析
		my @samplegroup = samplegroup($groupfile,$groupbase,$intm_result);
		
		for my $samplegroup (@samplegroup){
			
			my $path = "";
			
			($path) = (split /\//,$samplegroup)[-1] =~ /(.*).xls/ if ((split /\//,$samplegroup)[-1] =~ /_vs_/);
		
			my $num = `less $samplegroup | cut -f 2 | sort -u | wc -l`;
			my $num_sample =  `less $samplegroup | cut -f 1 | wc -l`;
			
			# =============== Metastats ===============
			if ($num_sample > $num and $num >= 2){
				system qq{mkdir -p $Community/Metastats} if not -d qq{$Community/Metastats};
				system qq{$Rscript $scriptdir/16S.Community.Metastats.r $taxonfile $samplegroup "$Community/Metastats" &>> $tmpPath/community.log};
			}
		
			# =============== Wilcoxon_ANOVA.DEA ===============
			if ($num_sample > $num and $num >= 3){
				system qq{mkdir -p $Community/ANOVA/$path} if not -d qq{$Community/ANOVA/$path};
				system qq{$Rscript $scriptdir/16S.Community.Wilcoxon_ANOVA.DEA.r $taxonfile $samplegroup P 0.05 "$Community/ANOVA/$path" &>> $tmpPath/community.log};
			}
		}
		
	}

	# =============== KronaPlot =============== 
	my $speciesfile = qq{$path{Community_StructurePath}/species/species.taxon.Abundance.xls};
	#system qq{$Rscript $scriptdir/16S.Community.KronaPlot.r $speciesfile $path{"KronaPlotPath"}};

	Krona($speciesfile, $path{"KronaPlotPath"});
	my @inputfile = `find $path{"KronaPlotPath"} -name *\_data_for_Krona.xls`;
	foreach my $inputfile (@inputfile) {
		chomp $inputfile;
		my ($sample) = $inputfile =~ m/([^\/]*)\_data_for_Krona.xls/;
		my $outputfile = qq{$path{"KronaPlotPath"}/$sample\.krona.html};
		system qq{$kronatools ktImportText $inputfile -o $outputfile &>> $tmpPath/community.log};
	}


	# check: tell you Where is the error! 
	foreach my $key (keys %path){

		my $errorpath = qq{$path{$key}};
		$key =~ s/Path//;
		
		if ($key eq "Barplot" or $key eq "Heatmap" or $key eq "TaxonTree"){
		
			my $minsize = `ls -l $errorpath|grep "pdf"|awk '{print \$5}'| sort -n |head -n 1`;

			print qq{WARINING/ERROR:
			**** Please check if mistakes when generate pdf file into $errorpath, 
			Please check it using 16S.Community.$key.r *** $!\n} if $minsize < 4500;   # exists & size>4500
		
		}else{
		
			foreach my $tax ("phylum", "class", "order", "family", "genus", "species"){
		
				if ($key eq "Community_Structure" or $key eq "Pieplot"){
			
					my $p = qq{$errorpath/$tax};
					my $minsize = `ls -l $p|grep "pdf"|awk '{print \$5}'|sort -n|head -n 1`;
					print qq{WARINING/ERROR:**** Please check if mistakes when generate pdf into $p, 
					Please check it using 16S.Community.$key.r *** $!\n} if $minsize < 4000;
			
				}
			
			}
		
		}

	}

}

sub Krona
{
	my $in = shift;
	my $out = shift;

	my $index = 0;
	my @sample;
	open IN ,"$in" || die "can not open $in";
	my $title = <IN>;
	my @arr = split /\t/, $title;
	$index++ until $arr[$index] eq "Abundance";
	@sample = @arr[1..$index-1];
	close IN;

	for my $sample (@sample){

		my $i = 0;
		my $j = 0;
		open IN ,"$in" || die "can not open $in";
		open OUT ,">$out/$sample\_data_for_Krona.xls" || die "can not write to $out/$sample\_data_for_Krona.xls";
		while (<IN>){
			chomp;
			$i++;
			if ($i == 1){
				my @str = split /\t/, $_;
				$j++ until $str[$j] eq $sample;
			}
			else{
			
				next if /^\bAll\b/;
				next if /^\bNo_Rank\b/;
				my @str = split /\t/, $_;
				my $taxon = join "\t", @str[$index+1..$index+6];
				
				print OUT qq{$str[$j]\t$taxon\t$str[0]\n};
			}
		}
		close IN;
		close OUT;
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
