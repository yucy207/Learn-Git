package otu_fungene_resample;

use Parallel::ForkManager;
use Sort::Key::Natural qw(natsort);

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
		
		my @groupname     = split /\//,$groupfile;
		my ($groupname)   = $groupname[-1] =~ /(.*).xls/;
		my $intm_result   = qq{$metadata->{intm_result}/$groupname};
		my $projectpath   = qq{$metadata->{result}/$groupname};
		
		if (-e qq{$intm_result/tmp/otu_fungene_resample.finish}){
			
			print qq{$groupname fungen抽平分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/fungene.finish} if -e qq{$intm_result/tmp/fungene.finish};
			system qq{rm $intm_result/tmp/fungene.log} if -e qq{$intm_result/tmp/fungene.log};
		}

		# Forks and returns the pid for the child:
		$manager->start and next;

		### fungene_otu_resample聚类注释
		fungene_otu_resample($intm_result,$projectpath,$groupfile,$metadata,$base);
		
		system qq{touch $intm_result/tmp/otu_fungene_resample.finish};
		print qq{$groupname otu_fungene_resample聚类注释完成\n};

		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{fungen聚类注释完成\n};
}

sub fungene_otu_resample
{
	
	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $trimOTUs    = qq{$metadata->{trimOTUs}};
	my $muscle      = qq{$base->{muscle_bin}};
	my $FastTree    = qq{$base->{FastTree_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}};

	my $OTUpath = qq{$projectpath/OTU};
	system qq{mkdir -p $OTUpath} if not -d $OTUpath;

	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	my $result  = qq{$metadata->{intm_result}};

	system qq{cp $groupfile $projectpath/sample_groups.xls};	
	#system qq{cp $result/tmp/otu.tax.$dissi.xls $tmpPath};
	#system qq{cp $result/tmp/otu.repseq.fasta $tmpPath};
	
	### Run
	# split_otu
	system qq{$Rscript $scriptdir/split_otu.r $result/tmp/otu.tax.$dissi.xls $groupfile $tmpPath &>> $tmpPath/resample.log};
	repseq("$result/tmp/otu.repseq.fasta", "$tmpPath/map.ID", "$tmpPath/otu.repseq.fasta");

	if ( $trimOTUs eq "T" ){
	
		### Run resample
		system qq{$Rscript $scriptdir/16S.Resample.r $tmpPath/otu.tax.$dissi.xls $groupfile $trimOTUs $tmpPath &>> $tmpPath/fungene.log};

		# OTU Modify
		#system qq{perl $scriptdir/16S.OTU.Modify.bak.pl $tmpPath/subsample_otu_table.txt $tmpPath/subsample_otu.wang.taxonomy $dissi $tmpPath};
		system qq{mv $tmpPath/subsample_otu.wang.taxonomy $tmpPath/subsample_otu.taxonomy};

		modify("$tmpPath/subsample_otu_table.txt", "$tmpPath/subsample_otu.taxonomy", $dissi, $tmpPath);
		
		# stat.RA for lefse
		system qq{$Rscript $scriptdir/stat2relativeAbundance.r $tmpPath/subsample_otu.tax.$dissi.stat.xls $tmpPath &>> $tmpPath/fungene.log &>> $tmpPath/fungene.log};
		
		# pick resample otu representative sequences
		#system qq{perl $scriptdir/16S.OTU.Repseq.pl $tmpPath/otu.repseq.fasta $tmpPath/map.ID > $tmpPath/subsample_otu.repseq.fasta};
		repseq("$tmpPath/otu.repseq.fasta", "$tmpPath/map.ID", "$tmpPath/subsample_otu.repseq.fasta");

		# tree
		# muscle
		system qq{$muscle muscle -maxiters 2 -maxmb 10000 -in $tmpPath/subsample_otu.repseq.fasta -out $tmpPath/subsample_otu.repseq.muscle.fasta &>> $tmpPath/fungene.log};
		system qq{$FastTree FastTree -nosupport $tmpPath/subsample_otu.repseq.muscle.fasta > $tmpPath/subsample_otu.repseq.fasta.tre 2>> $tmpPath/fungene.log};

		# PhyloTree
		# CirTreePlot
		#system qq{Rscript $scriptdir/16S.OTU.CirTreePlot.r $tmpPath/subsample_otu.tax.$dissi.xls $tmpPath &>> $tmpPath/fungene.log};
		# PhyloTree
		phylotree("$tmpPath/subsample_otu.tax.$dissi.xls", "$tmpPath/subsample_otu.repseq.fasta", 100, $tmpPath);
		# TopAbundantOTUs.seq.fasta.tre
		# muscle
		system qq{$muscle  muscle -maxiters 2 -in $tmpPath/TopAbundantOTUs.seq.fasta -out $tmpPath/TopAbundantOTUs.seq.muscle.fasta &>> $tmpPath/fungene.log};
		# FastTree generate tree
		system qq{$FastTree FastTree $tmpPath/TopAbundantOTUs.seq.muscle.fasta > $tmpPath/TopAbundantOTUs.seq.fasta.tre 2>> $tmpPath/fungene.log};
		
		# phylotree.pdf
		system qq{$Rscript $scriptdir/16S.OTU.PhyloTree.r "$tmpPath/subsample_otu.tax.$dissi.xls" $tmpPath/TopAbundantOTUs.seq.fasta.tre $tmpPath &>> $tmpPath/fungene.log};

		# OTU count
		system qq{$Rscript $scriptdir/16S.OTU.Count.r $tmpPath/subsample_otu.tax.$dissi.xls $tmpPath &>> $tmpPath/fungene.log};

		# mv
		#system qq{cp reads.clean.fasta $OTUpath};
		system qq{cp $tmpPath/otu.repseq.fasta $OTUpath};
		system qq{cp $tmpPath/otu.tax.$dissi.xls $OTUpath};
		system qq{cp $tmpPath/subsample_otu.repseq.fasta $OTUpath};
		system qq{cp $tmpPath/subsample_otu.tax.$dissi.xls $OTUpath};
		system qq{cp $tmpPath/subsample_otu.repseq.fasta.tre $OTUpath};
		system qq{cp $tmpPath/subsample_otu.tax.$dissi\_count.xls $OTUpath};
		system qq{cp $tmpPath/subsample_otu.tax.$dissi\_phylotree.pdf $OTUpath};
		system qq{cp $tmpPath/map.ID $OTUpath};
	}
	
	
	if ( $trimOTUs eq "F" ){
	
		# stat
		system qq{perl $scriptdir/otu_stat.pl $tmpPath/otu.tax.$dissi.xls $tmpPath};
	
		# stat.RA for lefse
		system qq{$Rscript $scriptdir/stat2relativeAbundance.r $tmpPath/otu.tax.$dissi.stat.xls $tmpPath};

		# PhyloTree
		phylotree("$tmpPath/otu.tax.$dissi.xls", "$tmpPath/otu.repseq.fasta", 100, $tmpPath);

		# TopAbundantOTUs.seq.fasta.tre
		# muscle
		system qq{$muscle muscle -maxiters 2 -in $tmpPath/TopAbundantOTUs.seq.fasta -out $tmpPath/TopAbundantOTUs.seq.muscle.fasta &>> $tmpPath/otu.log};
		
		# FastTree generate tree
		system qq{$FastTree FastTree $tmpPath/TopAbundantOTUs.seq.muscle.fasta > $tmpPath/TopAbundantOTUs.seq.fasta.tre 2>> $tmpPath/otu.log};
		
		# PhyloTree.pdf
		system qq{$Rscript $scriptdir/16S.OTU.PhyloTree.r $tmpPath/otu.tax.$dissi.xls $tmpPath/TopAbundantOTUs.seq.fasta.tre $tmpPath &>> $tmpPath/otu.log};
		
		# OTU count
		system qq{$Rscript $scriptdir/16S.OTU.Count.r $tmpPath/otu.tax.$dissi.xls $tmpPath &>> $tmpPath/otu.log};

		# otu.repseq.fasta.tre
		# muscle
		system qq{$muscle muscle -maxiters 2 -maxmb 10000 -in $tmpPath/otu.repseq.fasta -out $tmpPath/otu.repseq.muscle.fasta &>> $tmpPath/otu.log};
		# FastTree generate tree
		system qq{$FastTree FastTree -nosupport $tmpPath/otu.repseq.muscle.fasta > $tmpPath/otu.repseq.fasta.tre 2>> $tmpPath/otu.log};
			
		# mv 
		#system qq{cp $tmpPath/reads.clean.fasta $projectpath};
		system qq{cp $tmpPath/otu.repseq.fasta $OTUpath};
		system qq{cp $tmpPath/otu.tax.$dissi.xls $OTUpath};
		system qq{cp $tmpPath/otu.repseq.fasta.tre $OTUpath};
		system qq{cp $tmpPath/otu.tax.$dissi\_count.xls $OTUpath};
		system qq{cp $tmpPath/otu.tax.$dissi\_phylotree.pdf $OTUpath};
	}
	
}


sub modify
{
	#my ($otutable, $taxonomyfile, $dissimi, $outputpath) = @ARGV;
	my $otutable     = shift;
	my $taxonomyfile = shift;
	my $dissimi      = shift;
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
		my @tax = split /\t/;
		$taxs{$tax[0]} = $tax[1];

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
			
				#my ($x) = $data[$#data] =~ /([^;]+){$taxonclass[$i]}/;     # 注意返回值解列表
				my $x;
				if ( $data[$#data] =~ /__/){
					($x) = $data[$#data] =~ /$tags[$i]__([^__;]+)/;		# QIIME 
				}else{
					($x) = $data[$#data] =~ /([^;]+){$taxonclass[$i]}/;		# {}
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

	# stat.RA for lefse
	#system qq{Rscript /home/panrf/LearningPipeline/Metagenomics/stat2relativeAbundance.r $outputpath/$name.tax.$dissimi.stat.xls $outputpath};

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

1;