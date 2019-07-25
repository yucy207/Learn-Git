package picrust;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata      = shift;
	my $base          = shift;
	my $quantitation  = shift;
	# my $quantitation  = qq{$metadata->{quantitation}};
	my $level 		 = qq{$metadata->{level}};
	my $IR_copeis    = qq{$metadata->{IR_copeis}};
	my $DNA_quantity = qq{$metadata->{DNA_quantity}};

	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);

	for my $groupfile (@groups){
		
		my @groupname     = split /\//,$groupfile;
		# my ($groupname)   = $groupname[-1] =~ /(.*)\.[\w+]/;
		my ($groupname) = ($groupname[-1] =~ /::/)? $groupname[-1] =~ /(.*)\.\w\w\w::/ : $groupname[-1] =~ /(.*)\.[\w+]/;
		my ($groupfile, $select) = split /::/, $groupfile;
		my $intm_result   = qq{$metadata->{intm_result}/$groupname};
		my $intm_result_fasta   = qq{$metadata->{intm_result}};
		my $projectpath   = qq{$metadata->{result}/Report/$groupname};

		die "ERROR: *** There is no $groupfile, please check it: $!***\n" if not -e $groupfile;
		
		# Forks and returns the pid for the child:
		$manager->start and next;

		my $tmpPath = qq{$intm_result/tmp};
		system qq{mkdir -p $tmpPath} if not -d $tmpPath;

		if (uc($quantitation) eq "A" || uc($quantitation) eq "B"){
			Absolute_Quantitation($level, $projectpath, $intm_result, $intm_result_fasta, $groupfile, $select, $metadata,$base);
		}
		
		if (uc($quantitation) eq "R" || uc($quantitation) eq "B"){
			Relative_Quantitation($projectpath, $intm_result, $groupfile, $select, $metadata, $base);
		}


		# system qq{touch $tmpPath/picrust.finish};
		print qq{$groupname picrust分析完成\n};
		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{picrust分析已经运行完成\n};
}



sub Absolute_Quantitation{
	### prepare Absolute_Quantitation ###
	my $level 		= shift;
	my $projectpath = shift;
	my $intm_result = shift;
	my $intm_result_fasta = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $qiime		= qq{$base->{qiime_bin}};
	my $picrust     = qq{$base->{picrust_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $Rscript_xlsx = qq{$base->{Rscript_xlsx}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};
	my $otufasta    = qq{$base->{Greengenes}{fasta}};
	my $otutax      = qq{$base->{Greengenes}{Taxonomy}};
	my $IR_copeis    = qq{$metadata->{IR_copeis}};
	my $DNA_quantity = qq{$metadata->{DNA_quantity}};
	
	my $otufasta    = "/home/zhengy/bin/modules/database/SpikeIN_database/Greengenes_PI/97_otus.fasta";
	my $otutax      = "/home/zhengy/bin/modules/database/SpikeIN_database/Greengenes_PI/97_otu_taxonomy.txt";



	my $PIpath = qq{$projectpath/Absolute_Quantitation/PICRUSt};
	system qq{mkdir -p $PIpath} if not -d $PIpath;

	my $tmpPath = qq{$intm_result/Absolute_Quantitation/tmp/PI};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	
	###提取拷贝删除内参前的fasta
	system qq{cp $intm_result_fasta/tmp/Back_Before_IR_Detection/bakIR_reads.clean.fasta $intm_result/Absolute_Quantitation/tmp/reads.clean.fasta};
	my $fasta = qq{$intm_result/Absolute_Quantitation/tmp/reads.clean.fasta};  

	### Run 
	##1、fasta重命名
	rename_seq($fasta, $tmpPath);

	##2、clust otu
	my $id = `id`;
	my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;
	system qq{$qiime bash -c 'pick_closed_reference_otus.py -i $tmpPath/reads.clean.pi.fasta -r $otufasta -t $otutax -f -a -O 10 -o $tmpPath/pick_otus;chown -R $uid:$gid $tmpPath/pick_otus' >> $tmpPath/picrust.log};

	##3、计算OTU绝对丰度
	my $AA_Calculate = "$tmpPath/pick_otus/AA_Calculate";
	system qq{mkdir -p $AA_Calculate } if not -d $AA_Calculate;
	system qq{cp $tmpPath/pick_otus/otu_table.biom $AA_Calculate} if not -e "$AA_Calculate/otu_table.biom";

	system qq{$picrust biom convert -i $AA_Calculate/otu_table.biom -o $AA_Calculate/otu_table.txt --to-tsv --header-key taxonomy};
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.PI.Absolute_Abundance.r $AA_Calculate/otu_table.txt $IR_copeis $AA_Calculate &>> $tmpPath/picrust.log};

	###按group提取OTU表
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.PI.split_otu_newOTU.r $AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $groupfile $AA_Calculate/no.IR.absolute.otu.tax.0.03.xls &>> $tmpPath/picrust.log};

	if ($level eq "DNA") {
		system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.PI.Copies.unit.DNA.r $AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $DNA_quantity $AA_Calculate &>> $tmpPath/picrust.log};
		system qq{$picrust biom convert -i $AA_Calculate/otu_copies_unit_DNA.xls -o $tmpPath/pick_otus/otu_table.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy};
	}else{
		system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/16S.PI.Copies.unit.sample.r $AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $DNA_quantity $AA_Calculate &>> $tmpPath/picrust.log};
		system qq{$picrust biom convert -i $AA_Calculate/otu_copies_unit_sample.xls -o $tmpPath/pick_otus/otu_table.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy};
	}

	##4、normalize
	system qq{$picrust normalize_by_copy_number.py -i $tmpPath/pick_otus/otu_table.biom -o $tmpPath/normalized_otus.biom &>> $tmpPath/picrust.log};

	##5、Predict Functions
	system qq{$picrust predict_metagenomes.py -i $tmpPath/normalized_otus.biom -t ko -o $tmpPath/ko_predicted_metagenomes.biom &>> $tmpPath/picrust.log};
	system qq{$picrust predict_metagenomes.py -i $tmpPath/normalized_otus.biom -t cog -o $tmpPath/cog_predicted_metagenomes.biom &>> $tmpPath/picrust.log};

	##6、collapse
	system qq{$picrust categorize_by_function.py -f -i $tmpPath/cog_predicted_metagenomes.biom -c COG_Category -l 2 -o $PIpath/PICRUSt_cog.txt &>> $tmpPath/picrust.log};
	system qq{$picrust categorize_by_function.py -f -i $tmpPath/ko_predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o $PIpath/PICRUSt_kegg.txt &>> $tmpPath/picrust.log};

	system qq{sed -i 's/\#OTU ID/COG_Category/' $PIpath/PICRUSt_cog.txt};
	system qq{sed -i 's/\#OTU ID/KEGG_Pathways/' $PIpath/PICRUSt_kegg.txt};

	##7、prepare for STAMP
	pi2stamp("$PIpath/PICRUSt_cog.txt", $PIpath);
	pi2stamp("$PIpath/PICRUSt_kegg.txt", $PIpath);
	system qq{rm -rf $PIpath/PICRUSt_cog.txt};
	system qq{rm -rf $PIpath/PICRUSt_kegg.txt};

	##8、diff_analysis
	system qq{$Rscript $scriptdir/Absolute_Abundance/16S.PI_barplot_IR.r $PIpath/cog_for_stamp.txt $PIpath &>> $tmpPath/picrust.log};
	
	my $num = `less $groupfile | cut -f 2 | sort |  uniq | wc -l`;

	if ($num >=2 ){
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.PI_pair_diff_analysis_IR.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $groupfile "$select" $PIpath &>> $tmpPath/picrust.log};
	}

	if ($num >2 ){
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.PI_ANOVA_DEA_IR.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $groupfile P 0.05 $PIpath &>> $tmpPath/picrust.log};
	}


}

sub Relative_Quantitation{
	### prepare  Rbsolute_Quantitation ###
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $select      = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $qiime		= qq{$base->{qiime_bin}};
	my $picrust     = qq{$base->{picrust_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};
	my $otufasta    = qq{$base->{Greengenes}{fasta}};
	my $otutax      = qq{$base->{Greengenes}{Taxonomy}};
	

	my $PIpath = qq{$projectpath/Relative_Quantitation/PICRUSt};
	system qq{mkdir -p $PIpath} if not -d $PIpath;

	my $tmpPath = qq{$intm_result/Relative_Quantitation/tmp/PI};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $fasta = qq{$intm_result/Relative_Quantitation/tmp/reads.clean.fasta};

	rename_seq($fasta, $tmpPath);

	#2、clust otu
	my $id = `id`;
	my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;
	system qq{$qiime bash -c 'pick_closed_reference_otus.py -i $tmpPath/reads.clean.pi.fasta -r $otufasta -t $otutax -f -a -O 10 -o $tmpPath/pick_otus;chown -R $uid:$gid $tmpPath/pick_otus' &>> $tmpPath/picrust.log};

	#3、normalize
	system qq{$picrust normalize_by_copy_number.py -i $tmpPath/pick_otus/otu_table.biom -o $tmpPath/normalized_otus.biom &>> $tmpPath/picrust.log};

	#4、Predict Functions
	system qq{$picrust predict_metagenomes.py -i $tmpPath/normalized_otus.biom -t ko -o $tmpPath/ko_predicted_metagenomes.biom &>> $tmpPath/picrust.log};
	system qq{$picrust predict_metagenomes.py -i $tmpPath/normalized_otus.biom -t cog -o $tmpPath/cog_predicted_metagenomes.biom &>> $tmpPath/picrust.log};

	#5、collapse
	system qq{$picrust categorize_by_function.py -f -i $tmpPath/cog_predicted_metagenomes.biom -c COG_Category -l 2 -o $PIpath/PICRUSt_cog.txt &>> $tmpPath/picrust.log};
	system qq{$picrust categorize_by_function.py -f -i $tmpPath/ko_predicted_metagenomes.biom -c KEGG_Pathways -l 3 -o $PIpath/PICRUSt_kegg.txt &>> $tmpPath/picrust.log};

	system qq{sed -i 's/\#OTU ID/COG_Category/' $PIpath/PICRUSt_cog.txt};
	system qq{sed -i 's/\#OTU ID/KEGG_Pathways/' $PIpath/PICRUSt_kegg.txt};

	#6、prepare for STAMP
	pi2stamp("$PIpath/PICRUSt_cog.txt", $PIpath);
	pi2stamp("$PIpath/PICRUSt_kegg.txt", $PIpath);
	system qq{rm -rf $PIpath/PICRUSt_cog.txt};
	system qq{rm -rf $PIpath/PICRUSt_kegg.txt};

	#7、diff_analysis
	system qq{$Rscript $scriptdir/PI_barplot.r $PIpath/cog_for_stamp.txt $PIpath &>> $tmpPath/picrust.log};
	my $num = `less $groupfile | cut -f 2 | sort |  uniq | wc -l`;

	if ($num >=2 ){
		system qq{$Rscript $scriptdir/PI_pair_diff_analysis.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $groupfile "$select" $PIpath &>> $tmpPath/picrust.log};
	}
	if ($num >2 ){
		system qq{$Rscript $scriptdir/PI_ANOVA_DEA.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $groupfile P 0.05 $PIpath &>> $tmpPath/picrust.log};
	}

	
}


sub rename_seq
{
	my $fasta = shift @_;
	my $outputpath = shift @_;
	my ($name) = $fasta =~ /([^\/]*)\.fa/;
	open FA, $fasta or die "Can not open $fasta, please check it!\n";
	open OUT, qq{>$outputpath/$name\.pi.fasta} or die "Can not write into $outputpath, please check it!\n";
	$/ = "\n>"; 
	my %sample2ID;
	while(<FA>){
		
		chomp;
		next if /^\s*$/;
		
		my @data = split /\n/, $_, 2;
		my $id = (split /=/, $data[0])[1];
		$id =~ s/_/-/g;
		$id =~ s/;//;
		my $seq = $data[1];
		
		if (exists $sample2ID{$id}){
		
			$sample2ID{$id} = $sample2ID{$id} + 1;
		
		}else{
		
			$sample2ID{$id} = 1;
		
		}
		
		print OUT qq{>$id\_$sample2ID{$id}\n$seq\n};
		
	}

	$/ = "\n";
	close FA;
	close OUT;
}

sub pi2stamp 
{

	my $inputfile = shift @_;
	my $PIpath = shift @_;
	
	open IN, $inputfile or die "Can not open $inputfile, please check it!\n";
	
	if($inputfile =~ /cog/){
		open OUT,qq{>$PIpath/cog_for_stamp.txt} or die "Can not write into $PIpath, please check it!\n";
	}elsif($inputfile =~ /kegg/){
		open OUT,qq{>$PIpath/kegg_for_stamp.txt} or die "Can not write into $PIpath, please check it!\n";
	}
	
	while(<IN>){
		
		chomp;
		next if /biom/;
		my @data = split /\t/;
		my $str = join "\t", @data[0..($#data-1)];
		print OUT qq{$str\n};
		
	}	
}


sub splitfasta
{
	my $in   = shift;
	my $list = shift;
	my $out  = shift;

	open (LIST,"$list") || die "cannot open";
	my %hash = ();
	while(<LIST>){
		chomp;	
		my($id,$name) = split/\t/,$_;
		$hash{$id} = $name;
	}
	close LIST;

	open(IN,"$in") || die "cannot open";
	open(OUT,">$out");
	while(my $seq = <IN>){
		if ($seq =~ /^\>/){
			my ($id) = $seq =~ /barcodelabel=(.*);/;	
			if (exists $hash{$id}){	
				#$seq =~ s/barcodelabel=$id;/barcodelabel=$hash{$id};/;		
				print OUT "$seq";	
				my $line = <IN>;	
				print OUT "$line";

			}
		}	
	}
	close IN;
	close OUT;
}

1;
