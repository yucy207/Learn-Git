package otu_noresample;

use Sort::Key::Natural qw(natsort);

sub run
{
	my $metadata    = shift;
	my $base        = shift;

	my $fasta       = qq{$metadata->{result}/reads.clean.fasta};
	my $dissi       = qq{$metadata->{dissi}};
	my $type        = qq{$metadata->{type}};
	my $cutoff      = qq{$metadata->{cutoff}};
	my $trimOTUs    = qq{$metadata->{trimOTUs}};	
	my $database    = qq{$metadata->{database}};
	my $soft        = qq{$metadata->{soft}};
	my $soft_tax    = qq{$metadata->{soft_tax}};
	my $renamelist  = qq{$metadata->{rename_sample}};
	my $intm_result = qq{$metadata->{intm_result}};
	my $projectpath = qq{$metadata->{result}};

	my $usearch     = qq{$base->{usearch_bin}};
	my $mothur      = qq{$base->{mothur_bin}};
	my $qiime       = qq{$base->{qiime_bin}};
	my $fasta_number= qq{$base->{fasta_number_bin}};
	my $muscle      = qq{$base->{muscle_bin}};
	my $FastTree    = qq{$base->{FastTree_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}};
	
	#pre_check
	die "ERROR: There is no $fasta, Please Check It! $!\n" if not -e $fasta;
	
	if (-e qq{$intm_result/tmp/otu.finish}){
		print qq{OTU聚类已经运行完成!\n 是否重新运行? [y/n]};
		my $option = <STDIN>;
		$option =~ s/[\r\n]//g;
		return if $option ne "y";
		system qq{rm $intm_result/tmp/otu.finish} if -e qq{$intm_result/tmp/otu.finish};
		system qq{rm $intm_result/tmp/otu.log} if -e qq{$intm_result/tmp/otu.log};
	}

	### prepare 
	my $OTUpath = qq{$projectpath/OTU};
	system qq{mkdir -p $OTUpath} if not -d $OTUpath;
	
	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	if (-e qq{$renamelist}){
		renamefasta($fasta, $renamelist, "$tmpPath/reads.clean.fasta");
	}else{
		system qq{cp $fasta "$tmpPath/reads.clean.fasta"};
	}

	### Run 
	system qq{$usearch /usearch -fastx_uniques $tmpPath/reads.clean.fasta -fastaout $tmpPath/derep.fasta -sizeout &>> $tmpPath/otu.log};
	system qq{$usearch /usearch -sortbysize $tmpPath/derep.fasta -fastaout $tmpPath/sorted.uniq.fasta -minsize 2 &>> $tmpPath/otu.log};

	my $in = 100*$dissi;
	my $similarity = 1-$dissi;

	# 聚OTU && 去嵌合体
	if (uc($soft) eq "U"){
	
		if ($dissi == "0.03"){
		
			system qq{$usearch /usearch -cluster_otus $tmpPath/sorted.uniq.fasta -otus $tmpPath/otu.fasta -relabel OTU -uparseout $tmpPath/results.txt &>> $tmpPath/otu.log};
		
		}else{
		
			system qq{/home/panrf/Softwares/usearch-8.1/usearch8.1.1861_i86linux64 -cluster_otus $tmpPath/sorted.uniq.fasta -otu_radius_pct $in -otus $tmpPath/otu.fasta -relabel OTU -uparseout $tmpPath/results.txt &>> $tmpPath/otu.log};
		
		}
		
	}elsif(uc($soft) eq "N"){
	
		system qq{$usearch /usearch -unoise3 $tmpPath/sorted.uniq.fasta -zotus $tmpPath/zotu.fasta &>> $tmpPath/otu.log};
	    system qq{less $tmpPath/zotu.fasta | sed -e 's/Zotu/OTU/g' > $tmpPath/otu.fasta};
		
	}elsif (uc($soft) eq "Q"){
	
		my ($uid, $gid) = `id` =~ /uid=(\d+).+?gid=(\d+)/;
		system qq{$qiime bash -c 'pick_otus.py -i $tmpPath/sorted.uniq.fasta -o $tmpPath/picked_otus -s $similarity;chown -R $uid:$gid $tmpPath/picked_otus'};
		system qq{$qiime bash -c 'pick_rep_set.py -i $tmpPath/picked_otus/sorted.uniq_otus.txt -f $tmpPath/sorted.uniq.fasta -o $tmpPath/otu.tmp.fasta;chown -R $uid:$gid $tmpPath/otu.tmp.fasta' &>> $tmpPath/otu.log};
		system qq{$fasta_number $tmpPath/otu.tmp.fasta OTU > $tmpPath/otu.fasta 2>> $tmpPath/otu.log};
		#system qq{rm -rf $tmpPath/picked_otus $tmpPath/otu.tmp.fasta};
	}

	# OTU丰度表
	system qq{$usearch /usearch -otutab $tmpPath/reads.clean.fasta -otus $tmpPath/otu.fasta -otutabout $tmpPath/otu_table.txt -biomout $tmpPath/otu_table.json -mapout $tmpPath/map.txt -notmatched $tmpPath/unmapped.fa &>> $tmpPath/otu.log};
	
	#database
	my ($ChimeraFA, $TaxonomyFA, $Taxonomy, $Taxonomy_udb) = database($database, $type, $base);
	# 物种注释
	if (uc($soft_tax) eq "M"){
	
		# FASTA路径不能含"-",参数传递会中断
		system qq{$mothur bash -c 'cd $tmpPath; mothur "#classify.seqs(fasta=otu.fasta, template=$TaxonomyFA, taxonomy=$Taxonomy, processors=10, cutoff=0)"' &>> $tmpPath/otu.log};
		
		my $taxonomyfile = `ls $tmpPath | grep 'wang.taxonomy'`;
		chomp $taxonomyfile;    # !!!
		system qq{mv $tmpPath/$taxonomyfile $tmpPath/otu.taxonomy} if -e qq{$tmpPath/$taxonomyfile};
		
	}elsif (uc($soft_tax) eq "U"){
	
		system qq{$usearch /usearch -sintax $tmpPath/otu.fasta -db $Taxonomy_udb -tabbedout $tmpPath/otu.tax -strand both &>> $tmpPath/otu.log};
		
		format_tax("$tmpPath/otu.tax", "$tmpPath/otu.taxonomy");
	}

	# OTU Modify
	modify("$tmpPath/otu_table.txt", "$tmpPath/otu.taxonomy", $dissi, $cutoff, $tmpPath);
	
	# out otu representative sequences
	repseq("$tmpPath/otu.fasta", "$tmpPath/map.ID", "$tmpPath/otu.repseq.fasta");	
	
	system qq{touch $intm_result/tmp/otu.finish};
		
	print qq{OTU聚类注释完成\n};
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


sub database
{
	my $database = shift;
	my $type     = shift;
	my $base     = shift;

	#database
	$database = uc $database;
	my $ChimeraFA;
	my $TaxonomyFA;
	my $Taxonomy;
	my $Taxonomy_udb;
	if ($database eq "RDP"){

		if (uc $type eq "16S"){

			$ChimeraFA    = qq{$base->{RDP_16S}{ChimeraFA}};
			$TaxonomyFA   = qq{$base->{RDP_16S}{TaxonomyFA}};
			$Taxonomy     = qq{$base->{RDP_16S}{Taxonomy}};
			$Taxonomy_udb = qq{$base->{RDP_16S}{Taxonomy_udb}};

		}else{

			die "ERROR and Sorry: *** For RDP, only support 16S now, can not support your type option, please check it! ***\n";
		}

	}elsif ($database eq "SILVA"){

		if (uc $type eq "16S"){

			$ChimeraFA    = qq{$base->{SILVA_16S}{ChimeraFA}};
			$TaxonomyFA   = qq{$base->{SILVA_16S}{TaxonomyFA}};
			$Taxonomy     = qq{$base->{SILVA_16S}{Taxonomy}};
			$Taxonomy_udb = qq{$base->{SILVA_16S}{Taxonomy_udb}};

		}elsif (uc $type eq "18S"){

			$ChimeraFA    = qq{$base->{SILVA_18S}{ChimeraFA}};
			$TaxonomyFA   = qq{$base->{SILVA_18S}{TaxonomyFA}};
			$Taxonomy     = qq{$base->{SILVA_18S}{Taxonomy}};
			$Taxonomy_udb = qq{$base->{SILVA_18S}{Taxonomy_udb}};

		}else{
		
			die "ERROR and Sorry: *** For SILVA , only support 16S and 18S now, can not support your type option, please check it! ***\n";
		}

	}elsif ($database eq "UNITE"){

		$TaxonomyFA   = qq{$base->{UNITE_ITS}{TaxonomyFA}};
		$Taxonomy     = qq{$base->{UNITE_ITS}{Taxonomy}};
		$Taxonomy_udb = qq{$base->{UNITE_ITS}{Taxonomy_udb}};
		
		if (uc $type eq "ITS1"){

			$ChimeraFA  = qq{$base->{UNITE_ITS1}{ChimeraFA}};
		}elsif (uc $type eq "ITS2"){

			$ChimeraFA  = qq{$base->{UNITE_ITS2}{ChimeraFA}};
		}

		else{

			die "ERROR and Sorry: *** For UNITE , only support ITS1 and ITS2 now, can not support your type option, please check it! ***\n";
		}

	}elsif ($database eq "GREENGENE"){

		if (uc $type eq "16S"){

			$ChimeraFA    = qq{$base->{GREENGENE_16S}{ChimeraFA}};
			$TaxonomyFA   = qq{$base->{GREENGENE_16S}{TaxonomyFA}};
			$Taxonomy     = qq{$base->{GREENGENE_16S}{Taxonomy}};
			$Taxonomy_udb = qq{$base->{GREENGENE_16S}{Taxonomy_udb}};

			
		}else{
			die "ERROR and Sorry: *** For Greengene, only support 16S now, can not support your type option, please check it! ***\n";
		}

	}elsif ($database eq "PR2"){
		
		if (uc $type eq "18S"){

			$ChimeraFA    = qq{$base->{PR2_18S}{ChimeraFA}};
			$TaxonomyFA   = qq{$base->{PR2_18S}{TaxonomyFA}};
			$Taxonomy     = qq{$base->{PR2_18S}{Taxonomy}};
			$Taxonomy_udb = qq{$base->{PR2_18S}{Taxonomy_udb}};

		}else{
			die "ERROR and Sorry: *** For PR2 , only support 18S now, can not support your type option, please check it! ***\n";
		}

	}else{

			die qq{ERROR and Sorry：*** We can not support the database you want, please add it to 16S.OTU.pl if necessary ***\n};

	}

	#die qq{Can not find database file $ChimeraFA, please chack it!\n} if not -e $ChimeraFA;
	die qq{Can not find database file $TaxonomyFA, please chack it!\n} if not -e $TaxonomyFA;
	die qq{Can not find database file $Taxonomy, please chack it!\n} if not -e $Taxonomy;
	#die qq{Can not find database file $Taxonomy_udb, please chack it!\n} if not -e $Taxonomy_udb;

	return ($ChimeraFA, $TaxonomyFA, $Taxonomy, $Taxonomy_udb);
}

sub format_tax
{
	$tax = shift;
	$taxonomy  = shift;

	open(TAX, $tax) or die "Can't open $tax\n";
	open(OUT, ">$taxonomy") or die "Can't write $taxonomy\n";

	while (<TAX>) {
	    chomp;
		my ($id, $tax) = split /\t/, $_;
		
		#$tax =~ s/,/;/g;
		
		$tax =~ s/k:/{superkingdom\}/g;
		$tax =~ s/d:/{kingdom\}/g;
		$tax =~ s/p:/{phylum\}/g;
		$tax =~ s/c:/{class\}/g;
		$tax =~ s/o:/{order\}/g;
		$tax =~ s/f:/{family\}/g;
		$tax =~ s/g:/{genus\}/g;
		$tax =~ s/s:/{species\}/g;
		my @tax = split /,/, $tax;
		for my $t (@tax){
			my ($a) = $t =~ /(\{.*\})/;
			my ($b) = $t =~ /\}(.*)\(/;
			my ($c) = $t =~ /\((.*)\)/;
			$c *= 100;
			$c = sprintf "%.f", $c;
			$t = qq{$b$a($c)};
		}
		my $newtax = join ";",@tax;
		
		# print qq{$id\t$newtax;\n};
		print OUT qq{$id\t$newtax;\n};

	}
	close TAX;
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
			
				#my ($x) = $data[$#data] =~ /([^;]+){$taxonclass[$i]}/;     # 注意返回值解列表
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

1
