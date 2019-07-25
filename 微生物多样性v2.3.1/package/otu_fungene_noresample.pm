package otu_fungene_noresample;

use Parallel::ForkManager;
use Sort::Key::Natural qw(natsort);

sub run
{
	
	my $metadata    = shift;
	my $base        = shift;

	my $intm_result = qq{$metadata->{intm_result}};
	my $projectpath = qq{$metadata->{result}};

	my $fasta       = qq{$metadata->{result}/reads.clean.fasta};
	my $dissi       = qq{$metadata->{dissi}};
	my $type        = qq{$metadata->{type}};
	my $soft        = qq{$metadata->{soft}};
	my $database    = qq{$metadata->{database}};
	my $renamelist  = qq{$metadata->{rename_sample}};
	my $trimOTUs    = qq{$metadata->{trimOTUs}};
	
	my $seq_crumbs  = qq{$base->{seq_crumbs_bin}};
	my $usearch     = qq{$base->{usearch_bin}};
	my $mothur      = qq{$base->{mothur_bin}};
	my $muscle      = qq{$base->{muscle_bin}};
	my $FastTree    = qq{$base->{FastTree_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}};
	my $rdptools    = qq{$base->{rdptools}};
	my $blast       = qq{$base->{blast}};

	#pre_check
	die "ERROR: There is no $fasta, Please Check It! $!\n" if not -e $fasta;

	if (-e qq{$intm_result/tmp/fungene.finish}){
	
			print qq{$groupname fungen聚类分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/fungene.finish} if -e qq{$intm_result/tmp/fungene.finish};
			system qq{rm $intm_result/tmp/fungene.log} if -e qq{$intm_result/tmp/fungene.log};
	}

	### prepare
	my $FrameBotpath = qq{$projectpath/FrameBot};
	system qq{mkdir -p $FrameBotpath} if not -d $FrameBotpath;

	#my $OTUpath = qq{$projectpath/OTU};
	#system qq{mkdir -p $OTUpath} if not -d $OTUpath;

	my $tmpPath = qq{$intm_result/tmp};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	if (-e qq{$renamelist}){
		renamefasta($fasta, $renamelist, "$tmpPath/reads.clean.fasta");
	}else{
		system qq{cp $fasta "$tmpPath/reads.clean.fasta"};
	}

	#my $length
	system qq{$seq_crumbs calculate_stats $tmpPath/reads.clean.fasta > $tmpPath/reads.clean.stat};
	my $average = seq_stat("$tmpPath/reads.clean.stat");
	my $length = int(0.8*$average/3);
	# fungene database
	my ($database, $blastdb) = database($type, $base);

	### Run FrameBot
	system qq{$mothur mothur "#unique.seqs(fasta=$tmpPath/reads.clean.fasta)" &>> $tmpPath/fungene.log};

	derep_single("$tmpPath/reads.clean.unique.fasta", "$tmpPath/reads.clean.names", 1, $intm_result);

	my @files = `find $tmpPath -name '*_uniq.fasta'`;

	open Sta, qq{>$FrameBotpath/prot_seqStat.xls} or die "Can not write into $FrameBotpath, please check it!\n";
	my $MAX_PROCESSES = 30;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	foreach my $file (@files){

		chomp $file;
		my ($name) = $file =~ /\/([^\.\/]*)\_uniq.fasta/;
		system qq{mkdir -p $tmpPath/$name} if not -d "$tmpPath/$name";

		#Forks and returns the pid for the child:
		my $pid = $pm -> start and next;
		
		system qq{touch $tmpPath/$name/$name.finish} if -e "$tmpPath/$name/$name\_corr_prot.clean.upper.fasta";

		## FrameBot
		if (not -e qq{$tmpPath/$name/$name.finish}){
			## FrameBot
			system qq{java -jar /home/panrf/Softwares/RDPTools/FrameBot.jar framebot -o $tmpPath/$name/$name -l $length -i 0.8 -N $database $file &>> $tmpPath/fungene.log};
			#system qq{$rdptools FrameBot framebot -o $tmpPath/$name/$name -l $length -i 0.8 -N $database $file &>> $tmpPath/fungene.log} if not -e qq{$tmpPath/$name/$name\_corr_prot.fasta};
			
			### change threshold
			#system qq{mv $tmpPath/$name/$name\_corr_prot.fasta  $tmpPath/$name/$name\_corr_prot.fasta.bak} if -e qq{$tmpPath/$name/$name\_corr_prot.fasta};
			###system qq{perl $scriptdir/merge_framebot.pl -successed $tmpPath/$name/$name\_framebot.txt -failed $tmpPath/$name/$name\_failed_framebot.txt -len_threshold $length -identity_threshold 30 -outdir $tmpPath/$name};
			#merge_framebot("$tmpPath/$name/$name\_framebot.txt", "$tmpPath/$name/$name\_failed_framebot.txt", $length, 30, "$tmpPath/$name");
			#system qq{mv $tmpPath/$name/framebot.prot.fasta  $tmpPath/$name/$name\_corr_prot.fasta} if -e qq{$tmpPath/$name/framebot.prot.fasta};

			## split_size
			###system qq{perl $scriptdir/split_size.pl $tmpPath/$name/$name\_corr_prot.fasta $tmpPath/$name};
			split_size("$tmpPath/$name/$name\_corr_prot.fasta", "$tmpPath/$name/$name\_corr_prot.clean.fasta");

			## change case
			system qq{$seq_crumbs change_case $tmpPath/$name/$name\_corr_prot.clean.fasta -o $tmpPath/$name/$name\_corr_prot.clean.upper.fasta -a upper &>> $tmpPath/fungene.log};
			
			system qq{touch $tmpPath/$name/$name.finish};
		}

		## raw reads 
		my $rawseqs = `less -S $tmpPath/reads.clean.fasta|grep "barcodelabel=$name;"|wc -l`;
		chomp $rawseqs;
		## prot reads
		my $protseqs = `less -S $tmpPath/$name/$name\_corr_prot.clean.upper.fasta|grep '>'|wc -l`;
		chomp $protseqs;

		my $perc_prot = sprintf "%0.2f",100*$protseqs/$rawseqs;
		my $str = $name."\t".$rawseqs."\t".$protseqs."\t".$perc_prot;
		
		print Sta qq{$str\n};
		
		$pm->finish;   # Terminates the child process	
	}

	$pm->wait_all_children;  # blocks until all forked processes have exited

	#close Sta;
	system qq{sed -i '1iSamples\tRaw_reads\tprot_reads\tprot_reads(%)' $FrameBotpath/prot_seqStat.xls};
	system qq{cat $tmpPath/*/*.clean.upper.fasta > $tmpPath/translated_prot.fasta};
	system qq{cp $tmpPath/translated_prot.fasta $FrameBotpath/translated_prot.fasta};

	### Run OTU
	system qq{$usearch /usearch -fastx_uniques $tmpPath/translated_prot.fasta -fastaout $tmpPath/derep.fasta -sizeout &>> $tmpPath/fungene.log};
	system qq{$usearch /usearch -sortbysize $tmpPath/derep.fasta -fastaout $tmpPath/sorted.uniq.fasta -minsize 2 &>> $tmpPath/fungene.log};

	my $in = 100*$dissi;
	my $similarity = 1-$dissi;

	if (uc($soft) eq "U"){
	    system qq{/home/panrf/Softwares/usearch-8.1/usearch8.1.1861_i86linux64 -cluster_otus $tmpPath/sorted.uniq.fasta -otu_radius_pct $in -otus $tmpPath/otu.fasta -relabel OTU -uparseout $tmpPath/results.txt &>> $tmpPath/fungene.log};
		#system qq{$usearch /usearch -cluster_otus $tmpPath/sorted.uniq.fasta -otus $tmpPath/otu.fasta -relabel OTU -uparseout $tmpPath/results.txt &>> $tmpPath/fungene.log};
	}elsif(uc($soft) eq "N"){
		system qq{$usearch /usearch -unoise3 $tmpPath/sorted.uniq.fasta -zotus $tmpPath/zotu.fasta &>> $tmpPath/fungene.log};
		system qq{less $tmpPath/zotu.fasta | sed -e 's/Zotu/OTU/g' > $tmpPath/otu.fasta};
	}

	# otu_table.txt
	system qq{$usearch /usearch -otutab $tmpPath/translated_prot.fasta -otus $tmpPath/otu.fasta -otutabout $tmpPath/otu_table.txt -biomout $tmpPath/otu_table.json -mapout $tmpPath/map.txt -notmatched $tmpPath/unmapped.fa &>> $tmpPath/fungene.log};

	# blast
	#system qq{/home/panrf/Softwares/ncbi-blast-2.5.0+/bin/blastp -query $tmpPath/otu.fasta -db $blastdb -out $tmpPath/blastout.res -outfmt 7 -evalue 1e-5 -num_threads 20 -max_target_seqs 1 &>> $tmpPath/fungene.log};
	system qq{$blast blastp -query $tmpPath/otu.fasta -db $blastdb -out $tmpPath/blastout.res -outfmt 7 -evalue 1e-5 -num_threads 20 -max_target_seqs 1 &>> $tmpPath/fungene.log};
	
	#system qq{perl $scriptdir/blastp2tax.pl $tmpPath/blastout.res $database $tmpPath};
	blastp2tax("$tmpPath/blastout.res", $database, $tmpPath);

	# OTU Modify
	#system qq{perl $scriptdir/16S.OTU.Modify.bak.pl $tmpPath/otu_table.txt $tmpPath/otu.taxonomy $dissi $tmpPath/};
	modify("$tmpPath/otu_table.txt", "$tmpPath/otu.taxonomy", $dissi, $tmpPath);

	# out otu representative sequences
	#system qq{perl $scriptdir/16S.OTU.Repseq.pl $tmpPath/otu.fasta $tmpPath/map.ID > $tmpPath/otu.repseq.fasta};
	repseq("$tmpPath/otu.fasta", "$tmpPath/map.ID", "$tmpPath/otu.repseq.fasta");
	

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

	my $type = shift;
	my $base = shift;

	my $database;
	my $blastdb;

	# fungene database
	if ($type eq "nifh"){
		$database = qq{$base->{nifh}{database}};
		$blastdb  = qq{$base->{nifh}{blastdb}};
		
	}elsif ($type eq "nirS"){
		$database = qq{$base->{nirS}{database}};
		$blastdb  = qq{$base->{nirS}{blastdb}};
		
	}elsif ($type eq "nosZ"){
		$database = qq{$base->{nosZ}{database}};
		$blastdb  = qq{$base->{nosZ}{blastdb}};
		
	}elsif ($type eq "nirK"){
		$database = qq{$base->{nirK}{database}};
		$blastdb  = qq{$base->{nirK}{blastdb}};
		
	}elsif ($type eq "narG"){
		$database = qq{$base->{narG}{database}};
		$blastdb  = qq{$base->{narG}{blastdb}};
		
	}elsif ($type eq "norB"){
		$database = qq{$base->{norB}{database}};
		$blastdb  = qq{$base->{norB}{blastdb}};
		
	}elsif ($type eq "phoD"){
		$database = qq{$base->{phoD}{database}};
		$blastdb  = qq{$base->{phoD}{blastdb}};
		
	}else{
		die "============= ERROR:sorry, we can not support $type now! =============\n";
	}

	die qq{Can not find database file $database, please chack it!\n} if not -e $database;
	
	return ($database, $blastdb);
}

sub seq_stat
{
	my $seq = shift;

	open FILE, $seq or die "Can not open $seq, please check it!\n";

		my $totbases;
		my $totreads;
		my $average;
		my $Q20;
		my $Q30;
		my $a = 1;
		while(<FILE>){
			
			chomp;
			next if /^\s*$/;
			
			if (/tot. residues/){
				($totbases) = $_ =~ /^tot. residues: (.*)/;
			}
			if (/num. seqs/){
				($totreads) = $_ =~ /^num. seqs.: (.*)/;
			}
			if (/average/ and $a == 1){
				($average) = $_ =~ /^average: (.*)/;
				$a = $a + 1;
			}
			if (/Q20/){
				($Q20) = $_ =~ /^Q20: (.*)/;
			}
			if (/Q30/){
				($Q30) = $_ =~ /^Q30: (.*)/;
			}	
		}
		
		my $stat = qq{$totreads\t$Q20\t$Q30};
	close FILE;

	return ($average);
}

sub derep_single
{
	#my ($fasta, $names, $numbe, $outputpath) = @ARGV;
	my $fasta = shift;
	my $names =	shift;
	my $numbe = shift;
	my $outputpath = shift;

	#$outputpath = abs_path($outputpath);

	system qq{rm -rf $outputpath/tmp/*/*uniq.fasta};

	open FASTA, $fasta or die "ERROR: can not open $fasta, please check it!\n";
	open NAMES, $names or die "ERROR: can not open $names, please check it!\n";

	# readin reads.clean.unique.fasta
	my %id2seq;
	$/ = ">";
	while(<FASTA>){

		chomp;
		next if /^\s*$/;

		my @data = split /\n/;
		$id2seq{$data[0]} = $data[1];
		
	}
	close FASTA;
	$/ = "\n";

	# readin reads.clean.names
	while(<NAMES>){

		my %sample2size;
		my %sample2seqid;

		chomp;
		my @data = split /\s+/;
		my $id = $data[0];
		my @samples = split /,/, $data[1];
		next if (scalar @samples) < $numbe;		# size = 1

		foreach my $s (@samples){

			my ($sample) = $s =~ /barcodelabel=([^;]*);/;
			
			if (exists $sample2size{$sample}){

				$sample2size{$sample} = $sample2size{$sample} + 1;

			}else{

				$sample2size{$sample} = 1;
				$sample2seqid{$sample} = $s;
				system qq{mkdir -p $outputpath/tmp/$sample} if not -d qq{$outputpath/tmp/$sample};
				
			}

		}

		foreach my $key (keys %sample2size){

			open OUT, qq{>>$outputpath/tmp/$key/$key\_uniq.fasta} or die "Can not write into $outputpath/tmp/$key, please check it!\n";
			print OUT qq{>$sample2seqid{$key}size=$sample2size{$key};\n$id2seq{$id}\n};
			close OUT;

		}
		
	}
	close NAMES;
}

sub merge_framebot
{
	#-successed $tmpPath/$name/$name\_framebot.txt -failed $tmpPath/$name/$name\_failed_framebot.txt -len_threshold $length -identity_threshold 30 -outdir $tmpPath/$name
	my $successed = shift;
	my $failed = shift;
	my $len_threshold = shift;
	my $identity_threshold = shift;
	my $outdir = shift;

	system qq{mkdir -p $outdir} if not -d $outdir;
	parse_frame_out($successed);
	parse_frame_out($failed);

	sub parse_frame_out
	{
		my $file = shift;
		local $/ = ">";
		open FILE, $file or die "Can't open $file!\n";
		while (<FILE>) {
			chomp;
			next if /^\s*$/;
			my ($stats) = $_ =~ /STATS\s+(.+)/;
			my @querys  = $_ =~ /Query\s*\d+\s*(.+?)\s*\d+/g;
			my  @arr    = split /\s+/, $stats;
			my $seq_name   = $arr[1];
			my $align_len  = $arr[3];
			my $identity   = $arr[4];
			$prot_stats{$seq_name} = \@arr;

			next if $align_len < $len_threshold;
			next if $identity  < $identity_threshold;
			foreach my $x (@querys) {
				$x =~ s/\s+//g;
				$x =~ s/-//g;
				$prot_seqs{$seq_name} .= lc($x);
			}
		}
		close FILE;
	}

	open STATS, qq{>$outdir/framebot.stats.txt} or die "Can't open $outdir/framebot.stats.txt!\n";
	foreach my $x (keys %prot_stats) {
		my $line = join "\t", @{$prot_stats{$x}};
		print STATS qq{$line\n};
	}
	close STATS;

	open FASTA, qq{>$outdir/framebot.prot.fasta} or die "Can't open $outdir/framebot.prot.fasta!\n";
	foreach my $x (keys %prot_seqs) {
		my $line = join "\t", $prot_seqs{$x};
		print FASTA qq{>$x\n$line\n};
	}
	close FASTA;
}

sub split_size
{
	#split_size("$tmpPath/$name/$name\_corr_prot.fasta", "$tmpPath/$name/$name\_corr_prot.clean.fasta");
	
	my $fasta       = shift;
	my $clean_fasta = shift;
	#my $outputpath = shift;
	#my ($name) = $fasta =~ /([^\/]*)\.fa/;

	open FA, $fasta or die "Can not open $fasta, please check it!\n";
	open OUT, qq{>$clean_fasta} or die "Can not write into $clean_fasta, please check it!\n";
	#open OUT, qq{>$outputpath/$name\.clean.fasta} or die "Can not write into $outputpath, please check it!\n";
	$/ = ">"; 
	my %sample2ID;
	while(<FA>){
		
		chomp;
		next if /^\s*$/;
		
		my @data = split /\n/, $_, 2;
		
		my $seqid = (split /;/, $data[0])[0];
		$seqid =~ s/-/_/g;
		my $sample = (split /;/, $data[0])[1];
		my ($size) = $data[0] =~ /size=([^;]+);/;
		my $seq = $data[1];
		
		for(my $i = 1; $i <= $size; $i++){
			
			print OUT qq{>$seqid\_$i;$sample;\n$seq};
		}	
	}

	$/ = "\n"; 
	close FA;
	close OUT;
}

sub blastp2tax
{
	#my ($blastresout, $fasta, $outputpath) = @ARGV;
	my $blastresout = shift;
	my $fasta       = shift;
	my $outputpath  = shift;

	open FA, $fasta or die "Can not open $fasta, please check it!\n";
	my %id2tax;
	while(<FA>){

		chomp;
		next if $_ !~ /^>/;
		
		$_ =~ s/>//g;
		my @data = split /\s+/;
		my $id = $data[0];
		my $tax = $data[1];
		
		$id2tax{$id} = $tax;

	}
	close FA;

	open BLA, $blastresout or die "Can not open $blastresout, please check it!\n";
	open OUT, ">$outputpath/otu.taxonomy" or die "Can not write into $outputpath, please check it!\n";
	while(<BLA>){

		chomp;
		next if /^#/;

		my @data = split /\t/;
		if (exists $id2tax{$data[1]}){
			
			print OUT qq{$data[0]\t$id2tax{$data[1]}\n};

		}

	}
	close BLA;
	close OUT;
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