package package::snv;
use strict;
use warnings;

sub run{
	my $SNV=shift @_;
	my $format=shift @_;
	my $library=shift @_;
	my $annotation=shift @_;
	my $analysis=shift @_;
	my $filter=shift @_;
	my $run=shift @_;
	$SNV->set_row(0, 60);
	my @case=();@case=split /,/,$run->{"Case"}{1} if(exists $run->{"Case"}{1});
	my @control=();@control=split /,/,$run->{"Control"}{1} if(exists $run->{"Control"}{1});
	my @sample=(@case,@control);
	my @title=("SNV NO.","Gene","First Priority","SNP Calling Quality","Genotyping Quality","VCF Filter","Case_het","Case_hom");
	push @title,("Control_het","Control_hom") if(@control>0);
	push @title,(
        "Chrs",
		"Position",					
		"Ref Allele",
		"Alt Allele",
		"Gene Region",
		"Function",
		"Gene Strand Orientation",
		"Predicted Protein Variants",
		"Predicted Protein Variants(snpEFF)",
		"SNP ID",
		"Freq_Alt (1000g)",
		"1000g_chbs",
		"ExAC03",
		"ExAC03_EAS",
		"esp6500",
		"gnomAD_exome",
		"gnomAD_exome_EAS",
		"gnomAD_genome",
		"gnomAD_genome_EAS",
		"Hrcr1",
		"Kaviar",
		"NCI60",
		"GeneskyExonDB_Freq",
		"GeneskyGenomeDB_Freq",		
		"tfbsConsSites Score",
		"SIFT Score",
		"SIFT Score Pred",
		"POLYPhen V2 Score",
		"POLYPhen V2 Score Pred",
		"MutationTaster",
		"MutationTaster Pred",
		"Cadd_Raw",
		"Cadd_Phred",
		"Dann",		
		"Eigen",
		"dbscSNV_ADA_SCORE",
		"dbscSNV_RF_SCORE",				
		"OMIM",
		"HPO",
		"HGMD_site",
		"HGMD_site_class",
		"HGMD_gene",
		"InterVar",
		"InterVar_evidence",
		"Interpro_domain",
		"MGI",
		"ClinVar",
		"MalaCards",
		"COSMIC",
		"ICGC",
		"GO_BP",
		"GO_MF",
		"GO_CC",
		"KEGG_Pathway",
		"5' FLANKING SEQUENCE",
		"3' FLANKING SEQUENCE",
		"HOMOLOGY HITS",				
		"Case(0|1|2)",
	);
	push @title,"Control(0|1|2)" if(@control>0);
	push @title,"Alt Allele Freq";
	push @title,"Control Percentage" if(@control>0);
	if(not exists $run->{"Report"}{1} or $run->{"Report"}{1} ne "Simple"){
	  if(not exists $run->{"Report"}{1} or $run->{"Report"}{1} ne "HOMR"){ 
	  push @title,map{$_."_geno",$_."_reads"}@sample;	
	  }
	}
	my $row=0;
	for my $i(0..@title-1){
		$SNV->write($row,$i,$title[$i],$format->{"title"});
	}
	$row++;
	foreach my $id(sort{$a <=> $b}keys %{$filter}){
		foreach my $title(keys %{$filter->{$id}}){
			foreach my $alt(keys %{$filter->{$id}{$title}}){
				for my $i(0..@title-1){
					my $v=" ";
					$v=$library->{$title}{$alt}{$title[$i]} if(exists $library->{$title}{$alt}{$title[$i]});
					$v=$annotation->{$title}{$alt}{$title[$i]} if(exists $annotation->{$title}{$alt}{$title[$i]});
					$v=$analysis->{$title}{$alt}{$title[$i]} if(exists $analysis->{$title}{$alt}{$title[$i]});
					$v=$filter->{$id}{$title}{$alt}{$title[$i]} if(exists $filter->{$id}{$title}{$alt}{$title[$i]});
					$SNV->write($row,$i,$v,$format->{"normal"});
				}
				$row++;
			}
		}
	}
}


1