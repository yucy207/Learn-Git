package package::circRNA::circRNA_profile;

use Parallel::ForkManager;

sub run
{
	my $metadata = shift;
	my $base     = shift;
	my @samples  = split /,/, $metadata->{'samples'};
	my $organ1   = qq{$metadata->{organ}};
	my $organ    = qq{$base->{$organ1}{common_name}};

	my $circrna  = qq{$metadata->{project}/circRNA/circRNA_predict/result};
	my $result   = qq{$metadata->{project}/circRNA/circrna_profile};
	my $report   = qq{$metadata->{report}/05_circRNA_Analysis/02_circRNA_Expression_Profile};

	my $util       = qq{$base->{util}};
	my $rscript    = qq{$base->{rscript_bin}};
	my $pca_plot   = qq{$base->{pca_bin}};
	my $iresfinder = qq{$base->{iresfinder_bin}};
	my $cpat       = qq{$base->{cpat_bin}};

	system qq{mkdir -p $result/log}    if not -d qq{$result/log};
	system qq{mkdir -p $result/run}    if not -d qq{$result/run};
	system qq{mkdir -p $result/result} if not -d qq{$result/result};
	system qq{mkdir -p $report}        if not -d $report;
	system qq{mkdir -p $report/circRNA_Expression_Figure}  if not -d qq{$report/circRNA_Expression_Figure};

	if (-e qq{$report/circRNA.expression.xlsx}) {

		print qq{cicRNA 定量分析已经运行完成!\n};
		#return;
	}

	my %cnt = ();
	my %ids = ();
	foreach my $x (@samples) {
		open EXO, qq{$circrna/$x.predict_circRNA.xls} or die "Can't open $circrna/$x.predict_circRNA.xls\n";
		while (<EXO>) {
			chomp;
			my @arr= split /\t/;
			next if /^circRNA_ID/;
			my $unique = qq{$arr[1]:$arr[2]-$arr[3]:$arr[10]};
			my $gene   = qq{$arr[9]};
			$gene      =~ s/,$//;
			my $reads  = $arr[4];
			$cnt{$x}{$unique} = $reads; 
			$ids{$unique} = $gene;
		}
		close EXO;
	}

	open SAVE, qq{>$result/result/circRNA.expression.xls} or die "Can't open $result/result/circRNA.expression.xls";
	my $title = join "\t", @samples;
	print SAVE qq{circRNA_ID\tgene_id\t$title\n};
	foreach my $x (sort keys %ids) {
		my @counts = ();
		foreach my $sample (@samples) {
			my $exp = exists $cnt{$sample}{$x} ? $cnt{$sample}{$x}: 0;
			push @counts, $exp;
		}
		my $res = join "\t", @counts;
		print SAVE qq{$x\t$ids{$x}\t$res\n};
	}
	close SAVE;
	
	open SAVE, qq{>$result/result/circRNA.count.xls} or die "Can't open $result/result/circRNA.count.xls";
	my $title = join "\t", @samples;
	print SAVE qq{circRNA_ID\t$title\n};
	foreach my $x (sort keys %ids) {
		my @counts = ();
		foreach my $sample (@samples) {
			my $exp = exists $cnt{$sample}{$x} ? $cnt{$sample}{$x}: 0;
			push @counts, $exp;
		}
		my $res = join "\t", @counts;
		print SAVE qq{$x\t$res\n};
	}
	close SAVE;

	#plot common figures
   	system qq{$rscript Rscript $util/Deseq2_PCA.R $result/result/circRNA.count.xls $result/result/ &> $result/log/deseq2.log};
	system qq{$rscript Rscript $util/circRNA_density.R $result/result/circRNA.expression.xls $report/circRNA_Expression_Figure/density.pdf &> $result/log/density.log};
	
	#plot pca 
	if (exists $metadata->{'group'}) {

		my @groups = @{$metadata->{'group'}};
		my %unique = ();
		foreach my $x (@groups) {

			my $control = $x->[0];
			my $case    = $x->[1];
			my $control_samples = (split /;/, $x->[2])[0];
			my $case_samples    = (split /;/, $x->[2])[1];
			$unique{$control}   = $control_samples;
			$unique{$case}      = $case_samples;
		}

		open OUT, qq{>$result/result/sample.groups} or die "Can't open $result/result/sample.groups!\n";
		foreach my $x (keys %unique) {
			my @arr = split /\,/,$unique{$x};
			foreach my $y(@arr){
				print OUT qq{$y\t$x\n};
			}	
		}
		close OUT;

		system qq{$pca_plot Rscript $util/pca_plot.R $result/result/count.rlog.xls $result/result/sample.groups $report/circRNA_Expression_Figure/ &> $result/log/pca.log};
	}

	#添加注释
		#IRESfinder
	system qq{$iresfinder python /home/chengsy/softwares/IRESfinder-master/IRESfinder.py -f $circrna/circrnas_splice_transcript.fasta -o $result/result/IRESfinder.out &> $result/log/IRESfinder.log};
	my %ires_out = parse_ires(qq{$result/result/IRESfinder.out});

		#CPAT
	my @Organ = ("human","mouse","fly","zebrafish");
	my %cpat_out;
	if(grep {$_ eq $organ} @Organ){
		system qq{$cpat python /home/chengsy/softwares/CPAT-1.2.4/bin/cpat.py -g $circrna/circrnas_splice_transcript.fasta -d /home/chengsy/softwares/CPAT-1.2.4/dat/$organ\_logitModel.RData -x /home/chengsy/softwares/CPAT-1.2.4/dat/$organ\_Hexamer.tsv -o $result/result/cpat.out &> $result/log/cpat.log};
		%cpat_out = parse_cpat(qq{$result/result/cpat.out});
	}

	open EXP, qq{$result/result/circRNA.expression.xls} or die "Can't open $result/result/circRNA.expression.xls";
	my $title = <EXP>;
	chomp $title;
	my $head;
	if(grep {$_ eq $organ} @Organ){
		$head = "$title\tIRESfinder_Index\tIRESfinder_Score\tORF_size\tFickett_score\tHexamer_score\tcoding_prob";
	}else{
		$head = "$title\tIRESfinder_Index\tIRESfinder_Score";
	}
	open SAVE, qq{>$result/result/circRNA.expression.add.annotation.xls} or die "Can't make $result/result/circRNA.expression.add.annotation.xls";
	print SAVE "$head\n";
	while(my $line = <EXP>){
		$line =~s/[\r\n]//g;
		my @arr = split/\t/,$line;
		my $res;
		if(exists $ires_out{$arr[0]}){
			$res = "$line\t$ires_out{$arr[0]}";
		}
		if(exists $cpat_out{$arr[0]}){
			$res = "$res\t$cpat_out{$arr[0]}";
		}
		print SAVE "$res\n";
	}
	close EXP;
	close SAVE;

	system qq{perl $util/convert_txt2excel.pl $result/result/circRNA.expression.add.annotation.xls $report/circRNA.expression.xlsx};
	
	print qq{cicRNA 定量分析已经运行完成!\n};

}

sub parse_ires
{
	my $in = shift;
	my %ires = ();
	open IRES, $in or die "Can't open $in!";
	while(<IRES>){
		next if $_ =~/^ID/;
		$_ =~s/[\r\n]//g;
		my @line = split/\t/,$_;
		my ($id) = $line[0] =~/(.*)\|[\+|\-]\|/;
		$id =~s/\|/\-/;
		$ires{$id} = "$line[1]\t$line[2]";
	}
	close IRES;
	return %ires;
}

sub parse_cpat
{
	my $in = shift;
	my %cpat = ();
	open CPAT, $in or die "Can't open $in!";
	while(<CPAT>){
		next if $_ =~/^mRNA_size/;
		$_ =~s/[\r\n]//g;
		my @line = split/\t/,$_;
		my ($id) = $line[0] =~/(.*)\|[\+|\-]\|/;
		$id =~s/\|/\-/;
		$id =~s/CHR/chr/;
		my $value = join("\t",@line[2..$#line]);
		$cpat{$id} = $value;
	}
	close CPAT;
	return %cpat;
}

1;
