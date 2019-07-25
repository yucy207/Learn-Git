package positive_control;

use strict;
# use Sort::Key::Natural qw(natsort);
use Parallel::ForkManager;

sub run
{
	my $metadata    = shift;
	my $base        = shift;
	my $intm_result  = qq{$metadata->{intm_result}};
	my $result       = qq{$metadata->{result}};
	my $pc_result     = "$result/Positive_Control";

	if (-d $pc_result) {

		##删除reads.clean.fasta中的PC序列
		Remove_PCseq("$intm_result/tmp/IR_Detection/no.IR.reads.clean.fasta", "$pc_result/sample_groups.xls", $result);	

		##分析
		Absolute_Quantitation($intm_result, $pc_result, "$pc_result/sample_groups.xls", $metadata, $base);
		Relative_Quantitation($intm_result, $pc_result, "$pc_result/sample_groups.xls", $metadata, $base);

	}

	print qq{PC样本分析完成\n};

}


sub Absolute_Quantitation
{
	my $tmpPath     = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	# return if -e qq{$tmpPath/Positive_Control.finish};
	### prepare
	my $OTUpath = qq{$projectpath/Absolute_Quantitation/OTU};
	system qq{mkdir -p $OTUpath} if not -d $OTUpath;

	my $Barplot = qq{$projectpath/Absolute_Quantitation/Barplot};
	system qq{mkdir -p $Barplot} if not -d $Barplot;

	system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $tmpPath/tmp/AA_Calculate/otu_copies_unit_DNA.xls $projectpath/sample_groups.xls $OTUpath/otu_copies_unit_DNA.xls &>> $tmpPath/IR_PC.log};
	
	system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $tmpPath/tmp/AA_Calculate/no.IR.absolute.otu.tax.0.03.xls $projectpath/sample_groups.xls $OTUpath/copies.otu.tax.$dissi.xls &>> $tmpPath/IR_PC.log};
	system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $tmpPath/tmp/AA_Calculate/with.IR.absolute.otu.tax.0.03.xls $projectpath/sample_groups.xls $OTUpath/spike-in_copies.otu.tax.$dissi.xls &>> $tmpPath/IR_PC.log};

	###计算各分类水平
	system qq{$Rscript $scriptdir/Absolute_Abundance/16S.CommunityStructure.OTUSplit_PC.r $OTUpath/otu_copies_unit_DNA.xls $Barplot &>> $tmpPath/IR_PC.log};

	###绘图
	opendir Dir, $Barplot;
	my @fileslist = grep { /\.Abundance.xls/ } readdir(Dir);
	close Dir;

	foreach my $file (@fileslist) {

		my $dir = (split /\./, $file)[0];      # use taxon level as dir
		my $taxondir = "$Barplot/$dir";
		system qq{mkdir -p $taxondir} if not -d qq{$taxondir};
		system qq{mv $Barplot/$dir.* $taxondir};
		my $taxonfile = qq{$taxondir/$file};
		system qq{$Rscript $scriptdir/Absolute_Abundance/16S.Community.Barplot_PC.r $taxonfile $projectpath/sample_groups.xls $taxondir &>> $tmpPath/IR_PC.log};

	}

}


sub Relative_Quantitation
{
	my $tmpPath     = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	# return if -e qq{$tmpPath/Positive_Control.finish};

	### prepare
	my $OTUpath = qq{$projectpath/Relative_Quantitation/OTU};
	system qq{mkdir -p $OTUpath} if not -d $OTUpath;

	my $Barplot = qq{$projectpath/Relative_Quantitation/Barplot};
	system qq{mkdir -p $Barplot} if not -d $Barplot;

	system qq{$Rscript $scriptdir/Absolute_Abundance/split_otu_newOTU.r $tmpPath/tmp/IR_Detection/no.IR.otu.tax.0.03.xls $projectpath/sample_groups.xls $OTUpath/otu.tax.$dissi.xls &>> $tmpPath/PC.log};

	###计算各分类水平
	system qq{$Rscript $scriptdir/Absolute_Abundance/16S.CommunityStructure.OTUSplit_PC.r $OTUpath/otu.tax.$dissi.xls $Barplot &>> $tmpPath/PC.log};

	###绘图
	opendir Dir, $Barplot;
	my @fileslist = grep { /\.Abundance.xls/ } readdir(Dir);
	close Dir;

	foreach my $file (@fileslist) {

		my $dir = (split /\./, $file)[0];      # use taxon level as dir
		my $taxondir = "$Barplot/$dir";
		system qq{mkdir -p $taxondir} if not -d qq{$taxondir};
		system qq{mv $Barplot/$dir.* $taxondir};
		my $taxonfile = qq{$taxondir/$file};
		system qq{$Rscript $scriptdir/16S.Community.Barplot_PC.r $taxonfile $projectpath/sample_groups.xls $taxondir &>> $tmpPath/PC.log};

	}


}


sub Remove_PCseq 
{

	my($fasta, $name, $outputpath) = @_;

	my %hash = ();
	open NAME, $name;  ##两列
	while (<NAME>) {
		$_ =~ s/[\r\n]//g;
		my ($name1,$name2) = split /\t/, $_;
		$hash{$name1} = $name2;
	}
	close NAME;

	####删除reads.clean.fasta中的PC序列
	my %tmp = ();
	open FASTA, $fasta;
	open NEW, ">$outputpath/reads.clean.fasta";
	while (my $sample = <FASTA>) {
		my $seq = <FASTA>;
		my $sample1 = $sample;
		my @arry = split /=/, $sample;
		my ($samplename) = $sample =~ /barcodelabel\=(.+)\;/;
		next if exists $hash{$samplename};   ###PC样本跳过
		print NEW qq{$sample1$seq};

	}
	close FASTA;
	close NEW;

}




1

