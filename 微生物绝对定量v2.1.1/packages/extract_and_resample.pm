package extract_and_resample;

use Sort::Key::Natural qw(natsort);
use Parallel::ForkManager;

sub run
{
	my $metadata     = shift;
	my $base         = shift;
	my $quantitation = shift;
	# my $quantitation = qq{$metadata->{quantitation}};
	my $level 		 = qq{$metadata->{level}};
	my $IR_copeis    = qq{$metadata->{IR_copeis}};
	my $DNA_quantity = qq{$metadata->{DNA_quantity}};
	my $path         = qq{$metadata->{intm_result}};
	my $result       = qq{$metadata->{result}};

	####根据分组文件提取所需数据和抽平
	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);
	
	for my $groupfile (@groups){
		
		my @groupname = split /\//,$groupfile;
		my ($groupname) = ($groupname[-1] =~ /::/)? $groupname[-1] =~ /(.*)\.\w\w\w::/ : $groupname[-1] =~ /(.*)\.[\w+]/;
		my ($groupfile, $select) = split /::/, $groupfile;

		# print "$groupname\t$groupfile\n";
		# my ($groupname) = $groupname[-1] =~ /(.*)\.[\w+]/;
		my $intm_result = qq{$metadata->{intm_result}/$groupname};
		my $projectpath = qq{$metadata->{result}/Report/$groupname};
		#pre_check
		die "ERROR: There is no $groupfile, Please Check It! $!\n" if not -e $groupfile;
		
		# Forks and returns the pid for the child:
		$manager->start and next;
		my $tmpPath = qq{$intm_result/tmp};
		system qq{mkdir -p $tmpPath} if not -d $tmpPath;
		
		if (uc($quantitation) eq "A" or uc($quantitation) eq "B"){
			Absolute_Quantitation($level, $projectpath,$intm_result,$groupfile,$metadata,$base);
		}
		
		if (uc($quantitation) eq "R" or uc($quantitation) eq "B"){
			Relative_Quantitation($projectpath,$intm_result,$groupfile,$metadata,$base);
		}

		# system qq{touch $tmpPath/extract_and_resample.finish};
		print qq{$groupname 样本提取和抽平分析完成\n};

		$manager->finish;
		
	}

	$manager->wait_all_children;
	# print qq{样本提取和抽平分析完成\n};
}

###################################################子程序

sub Absolute_Quantitation{

	### prepare  Absolute_Quantitation ###
	my $level 		= shift;
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;
	
	my $dissi       = qq{$metadata->{dissi}};
	my $cutoff      = qq{$metadata->{cutoff}};
	my $trimOTUs    = qq{$metadata->{trimOTUs}};
	my $quantitation= qq{$metadata->{quantitation}};
	my $DNA_quantity= qq{$metadata->{DNA_quantity}};
	my $IR_resultpath = qq{$metadata->{intm_result}};
	
	my $muscle      = qq{$base->{muscle_bin}};
	my $FastTree    = qq{$base->{FastTree_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $Rscript_xlsx = qq{$base->{Rscript_xlsx}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	### prepare
	my $OTUpath = qq{$projectpath/Absolute_Quantitation/OTU};
	system qq{mkdir -p $OTUpath} if not -d $OTUpath;
	
	my $tmpPath = qq{$intm_result/Absolute_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;
	

	system qq{cp $groupfile $projectpath/sample_groups.xls} if not -e "$projectpath/sample_groups.xls";

	#####基于level水平选取OTU拷贝数表
	if ($level eq "DNA") {
		system qq{$Rscript $scriptdir/split_otu.r $IR_resultpath/tmp/AA_Calculate/otu_copies_unit_DNA.xls $groupfile $tmpPath &>> $tmpPath/IRnoresample.log};
	}else{

		system qq{$Rscript $scriptdir/split_otu.r $IR_resultpath/tmp/AA_Calculate/otu_copies_unit_sample.xls $groupfile $tmpPath &>> $tmpPath/IRnoresample.log};
	}
	repseq("$IR_resultpath/tmp/IR_Detection/no.IR.otu.repseq.fasta", "$tmpPath/map.ID", "$tmpPath/otu.repseq.fasta");

	
	###OTU stat
	###重新计算 otu.tax.$dissi.stat.xls
	system qq{perl $scriptdir/Absolute_Abundance/16S.OTU.stat_IR.pl $tmpPath/otu.tax.$dissi.xls $tmpPath & >> $tmpPath/IRnoresample.log};

	###level count
	system qq{$Rscript $scriptdir/16S.OTU.Count.r $tmpPath/otu.tax.$dissi.xls $tmpPath &>> $tmpPath/IRnoresample.log};

	# stat.RA for lefse
	# system qq{$Rscript $scriptdir/stat2relativeAbundance.r $tmpPath/otu.tax.$dissi.stat.xls $tmpPath &>> $tmpPath/IRnoresample.log};

	# PhyloTree
	phylotree("$tmpPath/otu.tax.$dissi.xls", "$tmpPath/otu.repseq.fasta", 100, $tmpPath);

	# TopAbundantOTUs.seq.fasta.tre
	# muscle
	system qq{$muscle  muscle -maxiters 2 -in $tmpPath/TopAbundantOTUs.seq.fasta -out $tmpPath/TopAbundantOTUs.seq.muscle.fasta &>> $tmpPath/IRnoresample.log};
	
	# FastTree generate tree
	system qq{$FastTree FastTree -nt $tmpPath/TopAbundantOTUs.seq.muscle.fasta > $tmpPath/TopAbundantOTUs.seq.fasta.tre 2>> $tmpPath/IRnoresample.log};
	
	# phylotree.pdf
	system qq{$Rscript $scriptdir/16S.OTU.PhyloTree.r $tmpPath/otu.tax.$dissi.xls $tmpPath/TopAbundantOTUs.seq.fasta.tre $tmpPath &>> $tmpPath/IRnoresample.log};
	
	# OTU count
	system qq{$Rscript $scriptdir/16S.OTU.Count.r $tmpPath/otu.tax.$dissi.xls $tmpPath &>> $tmpPath/IRnoresample.log};
	
	# subsample_otu.repseq.fasta.tre
	# muscle
	system qq{$muscle muscle -maxiters 2 -maxmb 10000 -in $tmpPath/otu.repseq.fasta -out $tmpPath/otu.repseq.muscle.fasta &>> $tmpPath/IRnoresample.log};
	
	# FastTree generate tree
	system qq{$FastTree FastTree -nosupport -nt $tmpPath/otu.repseq.muscle.fasta > $tmpPath/otu.repseq.fasta.tre 2>> $tmpPath/IRnoresample.log};

	
	# mv 
	# system qq{cp $result/tmp/AA_Calculate/Standard_curve_formula.xls $OTUpath};
	if ($level eq "DNA") {
		system qq{cp $tmpPath/otu.tax.$dissi.xls $OTUpath/otu_copies_unit_DNA.xls};  ####DNA水平结果
		system qq{cp $tmpPath/otu.tax.$dissi.xls $tmpPath/otu_copies_unit_DNA.xls}; ####DNA水平结果,ALpha多样性计算需要

		###若存在，只计算样本水平数值，不用于计算后续多样性分析
		my $ncol = `head -n 2 $DNA_quantity |awk 'NR==2 {print NF}'`;
		$ncol =~ s/[\r\n]//g;

		if ($ncol == 4) {
			my $Copies_Unit_Sample = qq{$OTUpath/Copies_Unit_Sample};
			system qq{mkdir -p $Copies_Unit_Sample} if not -d $Copies_Unit_Sample;
			system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Copies.1g.sample.r $tmpPath/otu_copies_unit_DNA.xls $DNA_quantity $groupfile $Copies_Unit_Sample &>> $tmpPath/IRcommunity.log};
		}

	}else{
		
		system qq{rm -rf $OTUpath/Copies_Unit_Sample} if -d "$OTUpath/Copies_Unit_Sample";
		system qq{cp $tmpPath/otu.tax.$dissi.xls $OTUpath/otu_copies_unit_sample.xls};
		system qq{cp $tmpPath/otu.tax.$dissi.xls $tmpPath/otu_copies_unit_sample.xls};####拷贝数矫正分析用
		system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $IR_resultpath/tmp/AA_Calculate/otu_copies_unit_DNA.xls $groupfile $OTUpath/otu_copies_unit_DNA.xls};  ####DNA水平结果
		system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $IR_resultpath/tmp/AA_Calculate/otu_copies_unit_DNA.xls $groupfile $tmpPath/otu_copies_unit_DNA.xls};  ####DNA水平结果,ALpha多样性计算需要	
	}

	###跟计算水平无关的结果
	system qq{cp $tmpPath/otu.repseq.fasta $OTUpath/otu.copies.repseq.fasta};	
	system qq{cp $tmpPath/otu.repseq.fasta.tre $OTUpath/otu.copies.repseq.fasta.tre};
	system qq{cp $tmpPath/otu.tax.$dissi\_count.xls $OTUpath/otu.copies.tax.$dissi\_count.xls};
	system qq{cp $tmpPath/otu.tax.$dissi\_phylotree.pdf $OTUpath/otu.copies.tax.$dissi\_phylotree.pdf};
	system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $IR_resultpath/tmp/IR_Detection/with.IR.otu.tax.0.03.xls $groupfile $OTUpath/spike-in_otu.tax.0.03.xls};      ####未计算绝对丰度，未删除内参序列的OTU表
	system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $IR_resultpath/tmp/IR_Detection/no.IR.otu.tax.0.03.xls $groupfile $OTUpath/otu.tax.$dissi.xls};               ####未计算绝对丰度，已删除内参序列的OTU表
	system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $IR_resultpath/tmp/AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $groupfile $OTUpath/otu.copies.tax.$dissi.xls};          ###绝对拷贝数
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/spikein_persent_select.r $IR_resultpath/tmp/IR_Detection/IR_reads_persent.xls $groupfile $OTUpath/spike-in_reads_persent.xls &>> $tmpPath/IRnoresample.log};               ###spikein占比统计表格拷贝
	system qq{$Rscript_xlsx $scriptdir/Absolute_Abundance/standard_curve_formula_select.r $IR_resultpath/tmp/AA_Calculate/Standard_curve_formula.xls $groupfile $OTUpath/Standard_curve_formula.xls &>> $tmpPath/IRnoresample.log};  ###标准曲线表格拷贝
 	system qq{cp $IR_resultpath/tmp/otu_genecopy_count.xls $tmpPath};    ##复制基因拷贝数表用于拷贝数矫正





	# system qq{touch $tmpPath/IRnoresample.finish};
	#print qq{Absolute_Quantitation分析完成\n};
}

sub Relative_Quantitation{
	### prepare  Relative_Quantitation ###
	my $projectpath = shift;
	my $intm_result = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $cutoff      = qq{$metadata->{cutoff}};
	my $trimOTUs    = qq{$metadata->{trimOTUs}};
	my $quantitation= qq{$metadata->{quantitation}};
	my $IR_resultpath = qq{$metadata->{intm_result}};
	
	my $muscle      = qq{$base->{muscle_bin}};
	my $FastTree    = qq{$base->{FastTree_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	my $OTUpath = qq{$projectpath/Relative_Quantitation/OTU};
	system qq{mkdir -p $OTUpath} if not -d $OTUpath;
	
	my $tmpPath = qq{$intm_result/Relative_Quantitation/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	#system qq{touch $tmpPath/resample.finish};
	# return if -e qq{$tmpPath/resample.finish};

	system qq{cp $groupfile $projectpath/sample_groups.xls};
	system qq{cp $IR_resultpath/tmp/IR_Detection/no.IR.otu.repseq.fasta $tmpPath/otu.repseq.fasta};
	splitfasta("$IR_resultpath/tmp/IR_Detection/no.IR.reads.clean.fasta", $groupfile, "$tmpPath/reads.clean.fasta");

	### Run
	# split_otu并生成mapID
	system qq{$Rscript $scriptdir/split_otu.r $IR_resultpath/tmp/IR_Detection/no.IR.otu.tax.0.03.xls $groupfile $tmpPath &>> $tmpPath/resample.log};

	my $sample_num = `wc -l $groupfile`;
	####一个样本不抽平
	if ($sample_num == 1) {
		system qq{cp $tmpPath/otu.tax.$dissi.xls $tmpPath/subsample_otu.tax.$dissi.xls};
		repseq("$IR_resultpath/tmp/IR_Detection/no.IR.otu.repseq.fasta", "$tmpPath/map.ID", "$tmpPath/otu.repseq.fasta");
		repseq("$IR_resultpath/tmp/IR_Detection/no.IR.otu.repseq.fasta", "$tmpPath/map.ID", "$tmpPath/subsample_otu.repseq.fasta");

	}else{
		# resample
		system qq{$Rscript $scriptdir/16S.Resample.r $tmpPath/otu.tax.$dissi.xls $groupfile $trimOTUs $tmpPath &>> $tmpPath/resample.log};

		# OTU Modify
		system qq{mv $tmpPath/subsample_otu.wang.taxonomy $tmpPath/subsample_otu.taxonomy};
		modify("$tmpPath/subsample_otu_table.txt", "$tmpPath/subsample_otu.taxonomy", $dissi, $cutoff, $tmpPath);

		# stat.RA for lefse
		system qq{$Rscript $scriptdir/stat2relativeAbundance.r $tmpPath/subsample_otu.tax.$dissi.stat.xls $tmpPath &>> $tmpPath/resample.log};

		# pick resample otu representative sequences
		repseq("$tmpPath/otu.repseq.fasta", "$tmpPath/map.ID", "$tmpPath/subsample_otu.repseq.fasta");
	}

	# PhyloTree
	phylotree("$tmpPath/subsample_otu.tax.$dissi.xls", "$tmpPath/subsample_otu.repseq.fasta", 100, $tmpPath);


	# TopAbundantOTUs.seq.fasta.tre
	# muscle
	system qq{$muscle  muscle -maxiters 2 -in $tmpPath/TopAbundantOTUs.seq.fasta -out $tmpPath/TopAbundantOTUs.seq.muscle.fasta &>> $tmpPath/resample.log};
	
	# FastTree generate tree
	system qq{$FastTree FastTree -nt $tmpPath/TopAbundantOTUs.seq.muscle.fasta > $tmpPath/TopAbundantOTUs.seq.fasta.tre 2>> $tmpPath/resample.log};
	
	# phylotree.pdf
	system qq{$Rscript $scriptdir/16S.OTU.PhyloTree.r $tmpPath/subsample_otu.tax.$dissi.xls $tmpPath/TopAbundantOTUs.seq.fasta.tre $tmpPath &>> $tmpPath/resample.log};
	
	# OTU count
	system qq{$Rscript $scriptdir/16S.OTU.Count.r $tmpPath/subsample_otu.tax.$dissi.xls $tmpPath &>> $tmpPath/resample.log};
	
	# subsample_otu.repseq.fasta.tre
	# muscle
	system qq{$muscle muscle -maxiters 2 -maxmb 10000 -in $tmpPath/subsample_otu.repseq.fasta -out $tmpPath/subsample_otu.repseq.muscle.fasta &>> $tmpPath/resample.log};
	
	# FastTree generate tree
	system qq{$FastTree FastTree -nosupport -nt $tmpPath/subsample_otu.repseq.muscle.fasta > $tmpPath/subsample_otu.repseq.fasta.tre 2>> $tmpPath/resample.log};

	# mv 
	system qq{cp $tmpPath/otu.repseq.fasta $OTUpath};
	system qq{cp $tmpPath/otu.tax.$dissi.xls $OTUpath};
	system qq{cp $tmpPath/subsample_otu.tax.$dissi.xls $OTUpath};
	system qq{cp $tmpPath/subsample_otu.repseq.fasta $OTUpath};
	system qq{cp $tmpPath/subsample_otu.repseq.fasta.tre $OTUpath};
	system qq{cp $tmpPath/subsample_otu.tax.$dissi\_count.xls $OTUpath};
	system qq{cp $tmpPath/subsample_otu.tax.$dissi\_phylotree.pdf $OTUpath};
	system qq{cp $tmpPath/map.ID $OTUpath};
	
	# system qq{touch $tmpPath/resample.finish};
	#print qq{Relative_Quantitation分析完成\n};
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


sub modify
{
	my $otutable     = shift;
	my $taxonomyfile = shift;
	my $dissimi      = shift;
	my $cutoff       = shift;
	my $outputpath   = shift;

	my ($name) = $otutable =~ /([^\/]*)_table/;

	open OTU, $otutable or die "Can not open $otutable, Please check it!\n";
	open TAX, $taxonomyfile or die "Can not open $taxonomyfile, Please check it!\n";
	open FILE, ">$outputpath/$name.tax.$dissimi.unmodify.xls" or die "Can not write to file, Please check it!\n";
	open MAP, ">$outputpath/map.ID" or die "Can not write to file, Please check it!\n";

	# Taxonomy
	my %taxs;
	while(<TAX>){

		chomp;
		next if /unknown|Unclassified/;
		
		$_ =~ s/"//g;
		$_ =~ s/\bkingdom\b/superkingdom/g;
		my @tax = split /\t/;
		#Unassigned
		my @arr = split /;/,$tax[1];
		my $flag = 1;

		foreach my $i(@arr){

			my ($value) = $i =~/.*\((\d+)\)/;
			if ($value < $cutoff){
				if ($i=~/\S+\{/){
					$i =~s/\S+\{/Unassigned\{/;
				}else{
					$i =~s/\{/Unassigned\{/;
				}
				
			}

			if ($i =~ /species/){
				if ($flag > 1){
					$i =~s/.*//;
				}
				$flag++;
			}
		}
		my $line = join";", grep { /\S/ } @arr;
		$taxs{$tax[0]} = $line;
	}
	close TAX;

	my $num = 0;
	# combine otu and Taxonomy
	while(<OTU>){

		chomp;
		
		if (/OTU.?I/){
		
			print FILE qq{$_\tOTUsize\tTaxonomy\n};
			
		}else{
		
			my @otu = split /\t/;
			
			if (exists $taxs{$otu[0]}){
			
				$num = $num + 1;
			
				my $OTUsize = 0;
				for my $i (1..$#otu){
					$OTUsize+=$otu[$i];
				}
				
				next if $OTUsize == 0;

				my $samplesotu = join "\t",@otu[1..$#otu];
				print FILE qq{OTU$num\t$samplesotu\t$OTUsize\t$taxs{$otu[0]}\n};
				print MAP qq{$otu[0]\tOTU$num\n};
			
			}
		
		}

	}
	close OTU;
	close FILE;

	# split Taxonomy
	open INPUT, qq{$outputpath/$name.tax.$dissimi.unmodify.xls} or die "Can not open this file: $!\n";
	open OUTPUT1, qq{>$outputpath/$name.tax.$dissimi.xls} or die "Can not write to file: $!\n";
	open OUTPUT2, qq{>$outputpath/$name.tax.$dissimi.for.stat.xls} or die "Can not write to file: $!\n";

	my @taxonclass = ("superkingdom", "phylum", "class", "order", "family", "genus", "species");
	my @tags = ("k", "p", "c", "o", "f", "g", "s");
	my %count;
	my %otutable;
	while(<INPUT>){

		chomp;
		if (/OTUsize/){
		
			print OUTPUT1 $_."\t".(join "\t",@taxonclass)."\n";     # 输出到 otu.tax.0.03.xls 表头
			
			my @title = split /\t/;
			my @samples = @title[1..($#title-2)];
			
			print OUTPUT2 "taxonomy"."\t".(join "\t",@samples)."\t"."Abundance"."\n";     # 输出到 otu.tax.0.03.for.stat.xls 表头
			
		}else{
			
			my @data = split /\t/; 
			$data[$#data] =~ s/"//g;
			
			my $key = join "\t",@data[0..$#data];
			my @taxon;
			
			# 将 otu.tax.0.03.unmodify.xls 最后一列按注释类别拆分, 输出到 otu.tax.0.03.xls
			for my $i (0..$#taxonclass){       # 数组下标
			
				my $x;
				if ( $data[$#data] =~ /__/){
					($x) = $data[$#data] =~ /$tags[$i]__([^__;]+)/;		# QIIME 
				}else{
					($x) = $data[$#data] =~ /([^;]+)\{$taxonclass[$i]\}/;		# {}
				}
				push @{$otutable{$key}}, $x;      # hash数组赋值
				$taxon[$i] = $x;
				
			}
			
			# otu.tax.0.03.stat.xls
			my ($Subscript) = grep{$taxon[$_] eq ""} 0..$#taxon;    # 数组中第一个空值元素的下标
		
			my @taxonomy;
			
			# 注释层中断,取中断位置前的注释内容
			for my $s (0..($Subscript-1)){
				$taxonomy[$s] = substr($tags[$s],0,1)."__".$taxon[$s];
			}
			
			# 注释层无中断,取全部注释内容
			if ($Subscript eq ""){
				for my $s (0..$#taxon){
					$taxonomy[$s] = substr($tags[$s],0,1)."__".$taxon[$s];
				}
			}

			my $str = join "|",@taxonomy;	
			if ( exists $count{$str} ) {       
			
				# 如果hash键值对已存在, 则对数据进行相加
				for my $i (0..($#data-2)) {
					$count{$str}->[$i] = $count{$str}->[$i] + $data[$i+1];
				}
				
			}else{
			
				@{$count{$str}} = @data[1..($#data-1)];      # 各样本分类数目及加和
			}
		}
	}

	foreach my $k (natsort (keys %otutable)){            # Sorting hash keys by Alphanumeric sort
		print OUTPUT1 $k."\t".(join "\t",@{$otutable{$k}})."\n";
	}

	foreach my $k (sort keys %count){
		print OUTPUT2 $k."\t".(join "\t",@{$count{$k}})."\n";
	}

	close INPUT;
	close OUTPUT1;
	close OUTPUT2;

	# prepare tree
	open FORSTAT, "$outputpath/$name.tax.$dissimi.for.stat.xls" or die "Can not open file: $!\n";
	open STAT, ">$outputpath/$name.tax.$dissimi.stat.xls" or die "Can not write to file: $!\n";

	my %stats;
	while(<FORSTAT>){

		chomp;
		next if /^\s+$/;
		
		if (/taxonomy/){
			print STAT qq{$_\n};

		}else{
		
			my @st = split /\t/;
			my @Tax = split /\|/, $st[0];
			
			if ($#Tax < 1){       # Only kingdom
				if (exists $stats{$Tax[0]}){
					for my $i (0..($#st-1)) {
						$stats{$Tax[0]}->[$i] = $stats{$Tax[0]}->[$i] + $st[$i+1];	
					}
				
				}else{
					@{$stats{$Tax[0]}} = @st[1..$#st];
				}
			
			}else{
			
				for my $t (0..$#Tax){
					my $test = join "|",@Tax[0..$t];      # 存在上层注释则数据合并
					if (exists $stats{$test}){    
						for my $i (0..($#st-1)) {
							$stats{$test}->[$i] = $stats{$test}->[$i] + $st[$i+1];	
						}
				
					}else{
						@{$stats{$test}} = @st[1..$#st];
					}	
				}	
			}	
		}
	}
	close FORSTAT;

	foreach my $s (sort keys %stats){

		print STAT $s."\t".(join "\t",@{$stats{$s}})."\n";
		
	}
	close STAT;


	system qq{rm -rf $outputpath/$name.tax.$dissimi.for.stat.xls};
	system qq{rm -rf $outputpath/$name.tax.$dissimi.unmodify.xls};
}


sub repseq
{
	my $fasta     = shift;
	my $mapid     = shift;
	my $rep_fasta = shift;

	open FASTA,$fasta or die "ERROR: can not open $fasta, please check it\n";
	open ID,$mapid or die "ERROR: can not open $mapid, please check it\n";

	my %IDS;
	while(<ID>){

		chomp;
		my @data = split /\t/;
		$IDS{$data[0]} = $data[1];

	}
	close ID;

	my %fastas;
	$/ = ">";       # 以 > 作为行的分隔
	while(<FASTA>){

		$_ =~ s/>//;
		my @data = split /[\r\n]/, $_, 2;
		my $otu = $data[0];
		
		next if /^\s+$/;

		if (exists $IDS{$otu}){
			$fastas{$IDS{$otu}} = $data[1];
		}
	}

	$/ = "\n";
	close FASTA;
	
	open OUT, ">$rep_fasta" or die "ERROR: can not write $rep_fasta, please check it\n";
	foreach my $k (natsort (keys %fastas)){       # Sorting hash keys by Alphanumeric sort
		print OUT qq{>$k\n$fastas{$k}};
	}
	close OUT;
}


sub phylotree
{
	my $otufile    = shift;
	my $fasta      = shift;
	my $N          = shift;
	my $outputpath = shift;
	
	my ($name) = $otufile =~ /([^\/]*)\.xls/;

	# OTUsize 所在的列
	my $colum = `head -1 $otufile | awk -F '\t' '{for(i=1;i<NF;i++){if(\$i ~ /OTUsize|Abundance/){print i}}}'`;
	chomp $colum;

	# Most N abundance OTUs
	my $OTUs = `less $otufile | sort -nrk $colum | head -n $N | awk -F '\\t' '{print \$1}'`;
	my @MAOTUs = split /\n/, $OTUs;

	open FASTA, $fasta or die "ERROR: can not open $fasta, please check it\n";

	my %fastas;
	$/ = ">";       # 以 > 作为行的分隔
	while(<FASTA>){

		$_ =~ s/>//;
		my @data = split /[\r\n]/, $_, 2;
		my $otu = $data[0];
		
		next if /^\s+$/;

		if (grep {$_ eq $otu} @MAOTUs){
			$fastas{$otu} = $data[1];
		}

	}

	$/ = "\n";
	close FASTA;

	open OTU, $otufile or die "ERROR: can not open $otufile, please check it\n";
	open TopOTU, qq{>$outputpath/TopAbundantOTUs.xls} or die "ERROR: can not write into $outputpath, please check it\n";
	open MAOTU, qq{>$outputpath/TopAbundantOTUs.seq.fasta};
	while(<OTU>){

		chomp;
		my @data = split /\t/;
		
		print TopOTU "$_\n" if /#OTU/;
		
		if (exists $fastas{$data[0]}) {
		
			print TopOTU "$_\n";
			print MAOTU ">".$data[0]."\n".$fastas{$data[0]};
			
		}

	}
	close OTU;
	close TopOTU;
	close MAOTU;

}


1

