#!/usr/bin/env perl

my ($organ, $result) = @ARGV;
my @ORGAN  = ("mm10","tair10");

my $A3SS_JCEC = "$result/A3SS.MATS.JCEC.txt";
my $A3SS_JC = "$result/A3SS.MATS.JC.txt";

my $A5SS_JCEC = "$result/A5SS.MATS.JCEC.txt";
my $A5SS_JC = "$result/A5SS.MATS.JC.txt";;

my $MXE_JCEC = "$result/MXE.MATS.JCEC.txt";
my $MXE_JC = "$result/MXE.MATS.JC.txt";

my $RI_JCEC = "$result/RI.MATS.JCEC.txt";
my $RI_JC = "$result/RI.MATS.JC.txt";

my $SE_JCEC = "$result/SE.MATS.JCEC.txt";
my $SE_JC = "$result/SE.MATS.JC.txt";

my @A3SS_title = qw/ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @A5SS_title = qw/ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @MXE_title  = qw/ID	GeneID	geneSymbol	chr	strand	1stExonStart_0base	1stExonEnd	2ndExonStart_0base	2ndExonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @RI_title   = qw/ID	GeneID	geneSymbol	chr	strand	riExonStart_0base	riExonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;
my @SE_title = qw/ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference/;

my %A3SS_set = readfile($A3SS_JCEC,$A3SS_JC);
my %A5SS_set = readfile($A5SS_JCEC,$A5SS_JC);
my %MXE_set;
my %RI_set = readfile($RI_JCEC,$RI_JC);
my %SE_set = readfile($SE_JCEC,$SE_JC);

#筛选JC和JCEC两种算法交集结果中FDR最显著top10

#A3SS
open OUT,">$result/A3SS.overlap.top10.txt";
my $title = join "\t",  @A3SS_title;
print OUT qq{$title\n};
my $count = 0;
foreach my $FDR(sort {$a <=> $b} keys %A3SS_set){
	    foreach my $id (keys %{$A3SS_set{$FDR}}){
	    	if ($count < 10) {
	        print OUT $A3SS_set{$FDR}{$id}."\n";
	    	$count ++;
		}
	}
}
close OUT;

#A5SS
open OUT,">$result/A5SS.overlap.top10.txt";
my $title = join "\t",  @A5SS_title;
print OUT qq{$title\n};
my $count = 0;
foreach my $FDR(sort {$a <=> $b} keys %A5SS_set){
	    foreach my $id (keys %{$A5SS_set{$FDR}}){
	    	if ($count < 10) {
	        print OUT $A5SS_set{$FDR}{$id}."\n";
	    	$count ++;
		}
	}
}
close OUT;

#MXE
open JCEC,"$MXE_JCEC" or die "Can't open $MXE_JCEC!\n";

my @id =();
while (<JCEC>) {
	$_ =~ s/[\r\n]//g;
	next if (/FDR/);
	my @jcec_arr = split /\t/, $_;
	my $id = $jcec_arr[0];
	my $FDR = $jcec_arr[21];
	if ($FDR < 0.05) {
		push @id, $id;
	}
}
open JC, "$MXE_JC" or die "Can't open $MXE_JC!\n";
while (<JC>) {
	$_ =~ s/[\r\n]//g;
	next if (/FDR/);
	my @jc_arr = split /\t/, $_;
	if (grep{$_ eq $organ}@ORGAN){		
		if(/chr(\d+|[XY])\t[+-]/){
			$jc_arr[3] = $1;
			$_ = join "\t",@jc_arr;
		}	
	}
	my $id = $jc_arr[0];
	my $FDR = $jc_arr[21];
	if($FDR < 0.05){
		if (grep {$_ eq $id} @id){
			$MXE_set{$FDR}{$id} = $_;
		}
	}
}
close JCEC;
close JC;

open OUT,">$result/MXE.overlap.top10.txt";
my $title = join "\t",  @MXE_title;
print OUT qq{$title\n};
my $count = 0;
foreach my $FDR(sort {$a <=> $b} keys %MXE_set){
	    foreach my $id (keys %{$MXE_set{$FDR}}){
	    	if ($count < 10) {
	        print OUT $MXE_set{$FDR}{$id}."\n";
	    	$count ++;
		}
	}
}
close OUT;

#RI
open OUT,">$result/RI.overlap.top10.txt";
my $title = join "\t",  @RI_title;
print OUT qq{$title\n};
my $count = 0;
foreach my $FDR(sort {$a <=> $b} keys %RI_set){
	    foreach my $id (keys %{$RI_set{$FDR}}){
	    	if ($count < 10) {
	        print OUT $RI_set{$FDR}{$id}."\n";
	    	$count ++;
		}
	}
}
close OUT;

#SE
open OUT,">$result/SE.overlap.top10.txt";
my $title= join "\t",  @SE_title;
print OUT qq{$title\n};
my $count = 0;
foreach my $FDR(sort {$a <=> $b} keys %SE_set){
	    foreach my $id (keys %{$SE_set{$FDR}}){
	    	if ($count < 10) {
	        print OUT $SE_set{$FDR}{$id}."\n";
	    	$count ++;
		}
	}
}
close OUT;


#取JC和JCEC两种算法结果文件中FDR<0.05的ID交集，保存JC算法的结果
sub readfile{
	my $JCEC = shift;
	my $JC = shift;

	my %set;
	
	open JCEC,"$JCEC" or die "Can't open $JCEC!\n";

	my @id =();
	while (<JCEC>) {
		$_ =~ s/[\r\n]//g;
		next if (/FDR/);
		my @jcec_arr = split /\t/, $_;
		my $id = $jcec_arr[0];
		my $FDR = $jcec_arr[19];
		if ($FDR < 0.05) {
			push @id, $id;
		}
	}
    open JC, "$JC" or die "Can't open $JC!\n";
	while (<JC>) {
		$_ =~ s/[\r\n]//g;
		next if (/FDR/);	
		my @jc_arr = split /\t/, $_;
		if (grep{$_ eq $organ}@ORGAN){		
			if(/chr(\d+|[XY])\t[+-]/){
				$jc_arr[3] = $1;
				$_ = join "\t",@jc_arr;
			}	
		}
		my $id = $jc_arr[0];
		my $FDR = $jc_arr[19];
		if($FDR < 0.05){
			if (grep {$_ eq $id} @id){
				$set{$FDR}{$id} = $_;
			}
		}
	}
	close JCEC;
	close JC;
	return %set;
}