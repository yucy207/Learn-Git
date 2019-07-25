package package::genotype;
use strict;
use warnings;

sub run{
	my $Genotype=shift @_;
	my $Call=shift @_;
	my $format=shift @_;
	my $library=shift @_;
	my $annotation=shift @_;
	my $analysis=shift @_;
	my $filter=shift @_;
	my $run=shift @_;
	$Genotype->set_row(0, 60);
	$Call->set_row(0, 60);
	my @case=();@case=split /,/,$run->{"Case"}{1} if(exists $run->{"Case"}{1});
	my @control=();@control=split /,/,$run->{"Control"}{1} if(exists $run->{"Control"}{1});
	my @sample=(@case,@control);
	my %tmp=();
	foreach my $id(keys %{$filter}){
		foreach my $title(keys %{$filter->{$id}}){
			$tmp{$title}=1;
		}
	}
	foreach(@sample){$tmp{$_}=1;}
	
	my %genotype_info=();
	foreach my $title(keys %{$analysis}) {
	    foreach my $alt(keys %{$analysis->{$title}}) {
		    foreach my $i(0..2) {
	            my $info_list  = $analysis->{$title}{$alt}{"Mutation $i"};
	        	my @infos      = split /,/, $info_list;
	        	foreach my $info(@infos) {
	        	    my ($sample, $Zygotes, $GENOTYPE, $Ref_dp, $Alt_dp, $Quality, $PL_HOMR, $PL_HET, $PL_HOMA) = split /:/, $info;
	        		$genotype_info{$title}{$alt}{$sample}{'Zygotes'}  = $Zygotes;
					$genotype_info{$title}{$alt}{$sample}{'GENOTYPE'} = $GENOTYPE;
					$genotype_info{$title}{$alt}{$sample}{'Ref_dp'}   = $Ref_dp;
					$genotype_info{$title}{$alt}{$sample}{'Alt_dp'}   = $Alt_dp;
					$genotype_info{$title}{$alt}{$sample}{'GQ'}       = $Quality;
					$genotype_info{$title}{$alt}{$sample}{'PL_HOMR'}  = $PL_HOMR;
					$genotype_info{$title}{$alt}{$sample}{'PL_HET'}   = $PL_HET;
					$genotype_info{$title}{$alt}{$sample}{'PL_HOMA'}  = $PL_HOMA;
	        	}
	        }
		}
	}
				
	my $row_Genotype=0;
	my $row_Call=0;
	my @title_Genotype=("SNV NO.","Gene","Ref Allele","Alt Allele",@sample);
	for my $i(0..@title_Genotype-1){
		if($i>=4){
			$Genotype->write_string($row_Genotype,$i,$title_Genotype[$i],$format->{"title"});
		}else{
			$Genotype->write($row_Genotype,$i,$title_Genotype[$i],$format->{"title"});
		}
	}
	$row_Genotype++;
	my @title_Call=("SNV NO.","Gene","SAMPLE","Zygotes","GENOTYPE","Ref Allele","Alt Allele","Ref_dp","Alt_dp","Alt Allele Freq","GQ","PL_HOMR","PL_HET","PL_HOMA");
	for my $i(0..@title_Call-1){
		$Call->write($row_Call,$i,$title_Call[$i],$format->{"title"});
	}
	$row_Call++;
	foreach my $id(sort{$a <=> $b}keys %{$filter}){
		foreach my $title(keys %{$filter->{$id}}){
			foreach my $alt(keys %{$filter->{$id}{$title}}){
				my $ID=$filter->{$id}{$title}{$alt}{"SNV NO."};
				my $Gene=$library->{$title}{$alt}{"Gene"};
				# my $bj301GeneInfo=" ";$bj301GeneInfo=$annotation->{$title}{$alt}{"bj301 GeneInfo"} if(exists $annotation->{$title}{$alt}{"bj301 GeneInfo"});
				my $GeneLength=" ";$GeneLength=$annotation->{$title}{$alt}{"Gene Length"} if(exists $annotation->{$title}{$alt}{"Gene Length"});
				my $ref=$library->{$title}{$alt}{"Ref Allele"};
				my %geno=("$ref/$ref"=>1,"$ref/$alt"=>1,"$alt/$ref"=>1,"$alt/$alt"=>1);
				$Genotype->write($row_Genotype,0,$ID,$format->{"normal"});
				$Genotype->write($row_Genotype,1,$Gene,$format->{"normal"});
				$Genotype->write($row_Genotype,2,$ref,$format->{"normal"});
				$Genotype->write($row_Genotype,3,$alt,$format->{"normal"});
				for my $i(0..@sample-1){
					next if(not exists $genotype_info{$title}{$alt}{$sample[$i]});
					my $SAMPLE   = $sample[$i];
					my $Zygotes  = $genotype_info{$title}{$alt}{$SAMPLE}{'Zygotes'};
					my $GENOTYPE = $genotype_info{$title}{$alt}{$SAMPLE}{'GENOTYPE'};
					my $Ref_dp   = $genotype_info{$title}{$alt}{$SAMPLE}{'Ref_dp'};
					my $Alt_dp   = $genotype_info{$title}{$alt}{$SAMPLE}{'Alt_dp'};
					my $Quality  = $genotype_info{$title}{$alt}{$SAMPLE}{'GQ'};
					my $PL_HOMR  = $genotype_info{$title}{$alt}{$SAMPLE}{'PL_HOMR'};
					my $PL_HET   = $genotype_info{$title}{$alt}{$SAMPLE}{'PL_HET'};
					my $PL_HOMA  = $genotype_info{$title}{$alt}{$SAMPLE}{'PL_HOMA'};
					my $depth    = $Ref_dp + $Alt_dp;
					next if(not exists $geno{$GENOTYPE});
					my $freq=($depth > 0) ? $Alt_dp/$depth : 'NA';
					if($Zygotes ne "HOMR" or not exists $run->{"Report"}{1} or $run->{"Report"}{1} ne "Simple"){
						my @array=($ID,$Gene,$SAMPLE,$Zygotes,$GENOTYPE,$ref,$alt,$Ref_dp,$Alt_dp,$freq,$Quality,$PL_HOMR,$PL_HET,$PL_HOMA);
						for my $j(0..@array-1){
							if($j==2){
								$Call->write_string($row_Call,$j,$array[$j],$format->{"normal"});
							}else{
								$Call->write($row_Call,$j,$array[$j],$format->{"normal"});
							}
						}
						$row_Call++;
					}
					$Genotype->write($row_Genotype,$i+4,$GENOTYPE,$format->{"normal"});
				}
				$row_Genotype++;
			}
		}
	}
}



1