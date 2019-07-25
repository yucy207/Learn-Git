package picrust;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata   = shift;
	my $base       = shift;

	my @groups     = @{$metadata->{'group'}};
	my $group_num  = @groups;
	my $manager    = new Parallel::ForkManager($group_num);
	my $renamelist = qq{$metadata->{rename_sample}};
	my $fasta      = qq{$metadata->{result}/reads.clean.fasta};
	my $result     = qq{$metadata->{intm_result}};
	die "ERROR: Can not find $fasta, Please Check It!\n" if not -e $fasta;

	my $tmpPath = qq{$result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	if (-e qq{$renamelist}){
		renamefasta($fasta, $renamelist, "$tmpPath/reads.clean.fasta") if not -e "$tmpPath/reads.clean.fasta";
	}else{
		system qq{cp $fasta "$tmpPath/reads.clean.fasta"} if not -e "$tmpPath/reads.clean.fasta";
	}

	for my $group (@groups){
		
		my $groupfile = (split /:/,$group)[0];
		my $groupbase = (split /:/,$group)[1];

		die "ERROR: Can not find $groupfile, Please Check It!\n" if not -e $groupfile;
		
		my @groupname   = split /\//,$groupfile;
		my ($groupname) = $groupname[-1] =~ /(.*).xls/;
		my $intm_result = qq{$metadata->{intm_result}/$groupname};
		my $projectpath = qq{$metadata->{result}/$groupname};
		
		
		if (-e qq{$intm_result/tmp/PI/picrust.finish}){
			print qq{$groupname PICRUSt分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/PI/picrust.finish};
			system qq{rm $intm_result/tmp/PI/picrust.log};
		}

		# Forks and returns the pid for the child:
		$manager->start and next;

		### PI分析
		picrust($intm_result,$projectpath,$groupfile,$groupbase,$metadata,$base);

		system qq{touch $intm_result/tmp/PI/picrust.finish};
		print qq{$groupname PI分析已经运行完成!\n};

		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{PI分析已经运行完成!\n};
}

sub picrust
{
	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $groupbase   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $qiime		= qq{$base->{qiime_bin}};
	my $picrust     = qq{$base->{picrust_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $Rscript_3_5_3     = qq{$base->{Rscript_3_5_3}};
	my $scriptdir   = qq{$base->{scriptdir}};
	my $otufasta    = qq{$base->{Greengenes}{fasta}};
	my $otutax      = qq{$base->{Greengenes}{Taxonomy}};
	

	my $PIpath = qq{$projectpath/PICRUSt};
	system qq{mkdir -p $PIpath} if not -d $PIpath;

	my $tmpPath = qq{$intm_result/tmp/PI};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	#system qq{touch $tmpPath/picrust.finish};
	return if -e qq{$tmpPath/picrust.finish};
	
	my $result = qq{$metadata->{intm_result}};
	#system qq{cp $result/tmp/reads.clean.fasta $tmpPath/reads.clean.fasta};
	splitfasta("$result/tmp/reads.clean.fasta", $groupfile, "$tmpPath/reads.clean.fasta");

	### Run 
	#1、rename fasta
	rename_seq("$tmpPath/reads.clean.fasta", $tmpPath);

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

	# #6、prepare for STAMP
	pi2stamp("$PIpath/PICRUSt_cog.txt", $groupfile, $PIpath);
	pi2stamp("$PIpath/PICRUSt_kegg.txt", $groupfile, $PIpath);
	system qq{rm -rf $PIpath/PICRUSt_*.txt};



	#7、 barplot
	system qq{$Rscript $scriptdir/PI_barplot.r $PIpath/cog_for_stamp.txt $PIpath &>> $tmpPath/picrust.log};

	###
	my $num_sample1 = `less $groupfile | cut -f 1 | wc -l`;
	if ($num_sample1 > 1 ){
		system qq{mkdir -p $PIpath/Heatmap};
		system qq{$Rscript_3_5_3 $scriptdir/16S.PI.heatmap.r $PIpath/cog_for_stamp.txt  $groupfile 30 "$PIpath/Heatmap" &>> $tmpPath/picrust.log};
		system qq{$Rscript_3_5_3 $scriptdir/16S.PI.heatmap.r $PIpath/kegg_for_stamp.txt  $groupfile 30 "$PIpath/Heatmap" &>> $tmpPath/picrust.log};

	}

	
	#8、diff_analysis
	###差异分析
	my @samplegroup = samplegroup($groupfile,$groupbase,$intm_result);
	
	for my $samplegroup (@samplegroup){
		
		my $path = "";
		
		($path) = (split /\//,$samplegroup)[-1] =~ /(.*).xls/ if ((split /\//,$samplegroup)[-1] =~ /_vs_/);
	
		my $num = `less $samplegroup | cut -f 2 | sort -u | wc -l`;
		my $num_sample =  `less $samplegroup | cut -f 1 | wc -l`;
		
		# system qq{$Rscript $scriptdir/PI_pair_diff_analysis.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $samplegroup $PIpath &>> $tmpPath/picrust.log} if ($num_sample > $num and $num >= 2);
		# system qq{$Rscript $scriptdir/PI_ANOVA_DEA.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $samplegroup P 0.05 $PIpath/$path &>> $tmpPath/picrust.log} if ($num_sample > $num and $num >= 3);
		system qq{$Rscript $scriptdir/PI_pair_diff_analysis.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $samplegroup $PIpath} if ($num_sample > $num and $num >= 2);
		system qq{$Rscript $scriptdir/PI_ANOVA_DEA.r $PIpath/cog_for_stamp.txt $PIpath/kegg_for_stamp.txt $samplegroup P 0.05 $PIpath/$path} if ($num_sample > $num and $num >= 3);

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
		s/[\r\n]//g;
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

sub renamefasta
{
	my $in = shift;
	my $list = shift;
	my $out = shift;
	
	open (LIST,"$list") || die "cannot open";

	my %hash = ();
	while(<LIST>){
		chomp;
		s/[\r\n]//g;
		my($id,$name) = split/\t/,$_;
		$hash{$id} = $name;
	}
	close LIST;


	open(IN,"$in") || die "cannot open";
	open(OUT,">$out");
	while(my $seq = <IN>){
		if ($seq =~ /^\>/){
			my ($id) = $seq =~/=(.*);/;
			if (exists $hash{$id}){
				$seq=~s/barcodelabel=$id;/barcodelabel=$hash{$id};/;
			}
			print OUT "$seq";
			my $line = <IN>;
			print OUT "$line";
		}	
	}
	close IN;
	close OUT;
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
	my $groupfile = shift @_;
	my $PIpath = shift @_;
	
	open IN, $inputfile or die "Can not open $inputfile, please check it!\n";
	
	if($inputfile =~ /cog/){
		open OUT,qq{>$PIpath/cog_for_stamp.txt} or die "Can not write into $PIpath, please check it!\n";

	}elsif($inputfile =~ /kegg/){
		open OUT,qq{>$PIpath/kegg_for_stamp.txt} or die "Can not write into $PIpath, please check it!\n";
	}
	
	my %hash;
	open GN, $groupfile;
	while (<GN>) {
		chomp;
		my ($samplename, $tmp) = split /\t/;
		my $samplename1 = $samplename;
		$samplename1 =~ s/_/-/g;
		$hash{$samplename1} = $samplename;
		
	}
	close GN;


	while(<IN>){
		
		chomp;
		next if /biom/;
		my @data = split /\t/;
		if ($data[0] eq "KEGG_Pathways" or $data[0] eq "COG_Category" ) {
			for (my $i = 1; $i < $#data; $i++) {
				$data[$i] = $hash{$data[$i]};
				# print $hash{$data[$i]}."\n";
			}
		}
		my $str = join "\t", @data[0..($#data-1)];
		print OUT qq{$str\n};
		
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
