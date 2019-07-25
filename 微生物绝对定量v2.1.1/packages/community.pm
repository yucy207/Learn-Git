package community;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata     = shift;
	my $base         = shift;
	my $quantitation = shift;
	# my $quantitation = qq{$metadata->{quantitation}};
	my $level 		 = qq{$metadata->{level}};
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
			Absolute_Quantitation($level,$projectpath,$intm_result,$groupfile,$select,$metadata,$base);
		}

		if (uc($quantitation) eq "R" || uc($quantitation) eq "B"){
			Relative_Quantitation($projectpath,$intm_result,$groupfile,$select,$metadata,$base);
		}

		print qq{$groupname Community分析完成\n};
		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{Community分析已经运行完成\n};
}




sub Absolute_Quantitation{
	### prepare  Absolute_Quantitation ###
	my $level 		= shift;
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $kronatools  = qq{$base->{kronatools_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	my $tmpPath = qq{$intm_result/Absolute_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $Community = qq{$projectpath/Absolute_Quantitation/Community};
	system qq{mkdir -p $Community} if not -d $Community;


	my $otufile     = qq{$tmpPath/otu.tax.$dissi.xls};
	my $otufasta    = qq{$tmpPath/otu.repseq.fasta};
	my $otustatfile = qq{$tmpPath/otu.tax.$dissi.stat.xls};
	my $DNA_quantity= qq{$metadata->{DNA_quantity}};

	# prepare output path
	my %path;
	foreach("Community_Structure", "Heatmap", "Pieplot", "Barplot", "TaxonTree", "KronaPlot", "Treebar", "ANOVA", "Kruskal_Wallis", "Bubble", "Wilcoxon", "ANOVA"){

		my $str = $_."Path";
		$path{$str} = $Community."/".$_;
		system qq{mkdir -p $path{$str}} if not -d qq{$path{$str}};
	}

	# =============== TaxonTree =============== 
	system qq{$Rscript $scriptdir/16S.Community.TaxonTree.r $otustatfile $path{"TaxonTreePath"} &>> $tmpPath/IRcommunity.log};

	# =============== OTUSplit =============== 
	system qq{$Rscript $scriptdir/16S.CommunityStructure.OTUSplit.r $otufile $path{"Community_StructurePath"} &>> $tmpPath/IRcommunity.log};

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

		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;
		
		my $taxonfile = qq{$taxondir/$file};
		
		# =============== CommunityStructure.Barplot =============== 
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.CommunityStructure.Barplot_IR.r $taxonfile $taxondir &>> $tmpPath/IRcommunity.log};
		
		# # =============== Pieplot =============== 
		system qq{$Rscript $scriptdir/16S.Community.Pieplot.r $taxonfile $groupfile $taxondir2 &>> $tmpPath/IRcommunity.log};

		# # =============== Barplot =============== 
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Barplot_IR.r $taxonfile $groupfile $path{"BarplotPath"} &>> $tmpPath/IRcommunity.log};
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Barplot.group_IR.r $taxonfile $groupfile $path{"BarplotPath"} &>> $tmpPath/IRcommunity.log} if ($num_sample > $num and $num != 1);
		
		# # # =============== Heatmap =============== 
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Heatmap_IR.r $taxonfile 100 $path{"HeatmapPath"} &>> $tmpPath/IRcommunity.log};
		
		# # # =============== TreeBar =============== 
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.TreeBar_IR.r $taxonfile $groupfile $path{"TreebarPath"} &>> $tmpPath/IRcommunity.log};

		# # # # =============== Kruskal_Wallis =============== 
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Kruskal_Wallis_IR.r $taxonfile $groupfile $path{"Kruskal_WallisPath"} &>> $tmpPath/IRcommunity.log} if ($num >2);


		# # # =============== Wilcoxon =============== 
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Wilcoxon_IR.r $taxonfile $groupfile "$select" $path{"WilcoxonPath"} &>> $tmpPath/IRcommunity.log} if ($num >1);


		# # # =============== Wilcoxon_ANOVA.DEA ===============
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Wilcoxon_ANOVA.DEA_IR.r $taxonfile $groupfile  P 0.05 $path{"ANOVAPath"} &>> $tmpPath/IRcommunity.log} if ($num >2);
		
		
		# # =============== Bubble =============== 
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Bubble_IR.r $taxonfile $groupfile 30 $path{"BubblePath"} &>> $tmpPath/IRcommunity.log};

		$pm->finish;


	}

	$pm->wait_all_children;


	# =============== KronaPlot =============== 
	my $speciesfile = qq{$path{Community_StructurePath}/species/species.taxon.Abundance.xls};
	system qq{$Rscript $scriptdir/16S.Community.KronaPlot.r $speciesfile $path{"KronaPlotPath"} &>> $tmpPath/IRcommunity.log};

	# Krona($speciesfile, $path{"KronaPlotPath"});

	# my @inputfile = `find $path{"KronaPlotPath"} -name *\_data_for_Krona.xls`;

	# foreach my $inputfile (@inputfile) {
	# 	chomp $inputfile;
	# 	my ($sample) = $inputfile =~ m/([^\/]*)\_data_for_Krona.xls/;
	# 	my $outputfile = qq{$path{"KronaPlotPath"}/$sample\.krona.html};
	# 	system qq{$kronatools ktImportText $inputfile -o $outputfile &>> $tmpPath/IRcommunity.log};
	# }

	# =============== Gene_Copies_Rectification && Copies_Unit_Sample =============== 
	my $Gene_Copies_Correction = qq{$Community/Gene_Copies_Correction};
	system qq{mkdir -p $Gene_Copies_Correction} if not -d $Gene_Copies_Correction;
	system qq{$Rscript $scriptdir/Absolute_Abundance/otu_genecopy_select.r $tmpPath/otu_genecopy_count.xls $tmpPath/otu_copies_unit_DNA.xls $Gene_Copies_Correction/otu_genecopy_count.xls};
	# system qq{cp $tmpPath/otu_genecopy_count.xls $Gene_Copies_Rectification};


	if ($level eq "DNA") {
		my $otufile1     = qq{$tmpPath/otu_copies_unit_DNA.xls};
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.rrnDB.otu.rectification.r $otufile1 "$Gene_Copies_Correction/otu_genecopy_count.xls" $Gene_Copies_Correction &>> $tmpPath/IRcommunity.log};

		if (-d "$projectpath/Absolute_Quantitation/OTU/Copies_Unit_Sample"){
			my $otu_copies = qq{$projectpath/Absolute_Quantitation/OTU/Copies_Unit_Sample/otu_copies_unit_sample.xls};
			system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.rrnDB.otu.rectification.r $otu_copies "$Gene_Copies_Correction/otu_genecopy_count.xls" $Gene_Copies_Correction &>> $tmpPath/IRcommunity.log};
		}

	}else{
		my $otufile1     = qq{$tmpPath/otu_copies_unit_DNA.xls};
		my $otufile2     = qq{$tmpPath/otu_copies_unit_sample.xls};
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.rrnDB.otu.rectification.r $otufile1 "$Gene_Copies_Correction/otu_genecopy_count.xls" $Gene_Copies_Correction &>> $tmpPath/IRcommunity.log};
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.rrnDB.otu.rectification.r $otufile2 "$Gene_Copies_Correction/otu_genecopy_count.xls" $Gene_Copies_Correction &>> $tmpPath/IRcommunity.log};
	}



	###----------------------检查结果-----------------------
	my $Wilcoxonresult = `ls $Community/Wilcoxon|wc -l`;
	if ($Wilcoxonresult == 0) {
		system qq{rm -rf $Community/Wilcoxon};
	}

	if ($num < 3) {
		system qq{rm -rf $Community/Kruskal_Wallis};
		system qq{rm -rf $Community/ANOVA};
	}

	if ($num_sample == 1) {
		system qq{rm -rf $Community/Bubble};
		system qq{rm -rf $Community/Heatmap};
		system qq{rm -rf $Community/TaxonTree};
		system qq{rm -rf $Community/Treebar};
	}

	# check: tell you Where is the error! 
	foreach my $key (keys %path){

		my $errorpath = qq{$path{$key}};
		$key =~ s/Path//;
		
		if ($key eq "Barplot" or $key eq "Heatmap" or $key eq "TaxonTree"){
		
			my $minsize = `ls -l $errorpath|grep "pdf"|awk '{print \$5}'| sort -n |head -n 1`;

			print qq{WARINING/ERROR:
			**** Please check if mistakes when generate pdf file into $errorpath, 
			Please check it using 16S.Community.$key\_IR.r *** $!\n} if $minsize < 4500;   # exists & size>4500
		
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

sub Relative_Quantitation {
	### prepare  Relative_Quantitation ###
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $kronatools  = qq{$base->{kronatools_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};
	
	my $tmpPath = qq{$intm_result/Relative_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $Community = qq{$projectpath/Relative_Quantitation/Community};
	system qq{mkdir -p $Community} if not -d $Community;


	my $otufile     = qq{$tmpPath/subsample_otu.tax.$dissi.xls};
	my $otustatfile = qq{$tmpPath/subsample_otu.tax.$dissi.stat.xls};

	# prepare output path
	my %path;
	foreach("Community_Structure", "Heatmap", "Pieplot", "Barplot", "TaxonTree", "KronaPlot", "Treebar", "Metastats", "Bubble", "Kruskal_Wallis", "Wilcoxon", "ANOVA"){

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
		my $taxondir3 = qq{$path{EnterotypePath}/$dir};
		
		system qq{mkdir -p $taxondir} if not -d qq{$taxondir};
		system qq{mkdir -p $taxondir2} if not -d qq{$taxondir2};
		
		system qq{mv $path{Community_StructurePath}/$dir.* $taxondir};

		# Forks and returns the pid for the child:
		my $pid = $pm->start and next;
		
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

		# =============== Metastats =============== 
		system qq{$Rscript $scriptdir/16S.Community.Metastats.r $taxonfile $groupfile "$select" $path{"MetastatsPath"} &>> $tmpPath/community.log} if ($num >1);

		# # # =============== Wilcoxon =============== 
		system qq{$Rscript $scriptdir/16S.Community.Wilcoxon.r $taxonfile $groupfile "$select" $path{"WilcoxonPath"} &>> $tmpPath/community.log} if ($num >1);

		# # # # =============== Kruskal_Wallis =============== 
		system qq{$Rscript $scriptdir/16S.Community.Kruskal_Wallis.r $taxonfile $groupfile $path{"Kruskal_WallisPath"} &>> $tmpPath/community.log};

		# =============== Wilcoxon_ANOVA.DEA ===============
		system qq{$Rscript $scriptdir/16S.Community.Wilcoxon_ANOVA.DEA.r $taxonfile $groupfile  P 0.05 $path{"ANOVAPath"} &>> $tmpPath/community.log} if ($num >2);
		
		# =============== Bubble =============== 
		system qq{$Rscript $scriptdir/16S.Community.Bubble.r $taxonfile $groupfile 30 $path{"BubblePath"} &>> $tmpPath/community.log};

		$pm->finish;
	}

	$pm->wait_all_children;

	###------------------------检查结果-------------------
	my $Wilcoxonresult = `ls $Community/Wilcoxon|wc -l`;
	if ($Wilcoxonresult == 0) {
		system qq{rm -rf $Community/Wilcoxon};
	}

	my $Metastatsresult = `ls $Community/Metastats|wc -l`;
	if ($Metastatsresult == 0) {
		system qq{rm -rf $Community/Metastats};
	}


	if ($num < 3) {
		system qq{rm -rf $Community/Kruskal_Wallis};
		system qq{rm -rf $Community/ANOVA};
	}
	

	if ($num_sample == 1) {
		system qq{rm -rf $Community/Bubble};
		system qq{rm -rf $Community/Heatmap};
		system qq{rm -rf $Community/TaxonTree};
		system qq{rm -rf $Community/Treebar};
	}


	# =============== KronaPlot =============== 
	my $speciesfile = qq{$path{Community_StructurePath}/species/species.taxon.Abundance.xls};
	system qq{$Rscript $scriptdir/16S.Community.KronaPlot.r $speciesfile $path{"KronaPlotPath"} &>> $tmpPath/community.log};

	# Krona($speciesfile, $path{"KronaPlotPath"});

	# my @inputfile = `find $path{"KronaPlotPath"} -name *\_data_for_Krona.xls`;

	# foreach my $inputfile (@inputfile) {
	# 	chomp $inputfile;
	# 	my ($sample) = $inputfile =~ m/([^\/]*)\_data_for_Krona.xls/;
	# 	my $outputfile = qq{$path{"KronaPlotPath"}/$sample\.krona.html};
	# 	system qq{$kronatools ktImportText $inputfile -o $outputfile &>> $tmpPath/community.log};
	# }

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
	open IN ,"$in" || die "can open $in";
	my $title = <IN>;
	my @arr = split /\t/, $title;
	$index++ until $arr[$index] eq "Abundance";
	@sample = @arr[1..$index-1];
	close IN;

	for my $sample (@sample){

		my $i = 0;
		my $j = 0;
		open IN ,"$in" || die "can open $in";
		open OUT ,">$out/$sample\_data_for_Krona.xls" || die "can write to $out/$sample\_data_for_Krona.xls";
		while (<IN>){
			chomp;
			$i++;
			if ($i == 1){
				my @str = split /\t/, $_;
				$j++ until $str[$j] eq $sample;
			}
			else{
			
				next if /\bAll\b/;
				next if /\bNo_Rank\b/;
				my @str = split /\t/, $_;
				my $taxon = join "\t", @str[$index+1..$index+6];
				
				print OUT qq{$str[$j]\t$taxon\t$str[0]\n};
			}
		}
		close IN;
		close OUT;
	}
}

1;
