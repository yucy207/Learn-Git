package package::filter;
use strict;
use warnings;
$|=1; 

sub run{
	my $library=shift @_;
	my $annotation=shift @_;
	my $analysis=shift @_;
	my $run=shift @_;
	my @case=();@case=split /,/,$run->{"Case"}{1} if(exists $run->{"Case"}{1});
	my @control=();@control=split /,/,$run->{"Control"}{1} if(exists $run->{"Control"}{1});
	my @sample=(@case,@control);print "Lost Sample\n" if(@sample==0); 
	my %regionfilter=RunRegion($run);###区域过滤
	my %genefilter=RunFilterTarget($run,'Gene');# 基因过滤
	my %chrfilter=RunFilterTarget($run,'Chrs');# 染色体过滤
	my %gene=();
	my %geneMut=();
	my %geneAutoComHet=();
	my @titles=keys %$library;
	my $count=0;
	my $allsnvnum=@titles;
	foreach my $title(keys %{$library}){
	    $count++;
		my $p=sprintf "%0.4f",$count/$allsnvnum;
		print "\rFiltering $p";
		foreach my $alt(keys %{$library->{$title}}){
			next if(not exists $analysis->{$title}{$alt});
			my %mut=();
			foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 0"}){$mut{(split /:/,$_)[0]}="0:$_" if($_=~ /\w/);}
			foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 1"}){$mut{(split /:/,$_)[0]}="1:$_" if($_=~ /\w/);}
			foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 2"}){$mut{(split /:/,$_)[0]}="2:$_" if($_=~ /\w/);}
			my ($hc0,$hc1,$hc2)=(0,0,0);
			my ($sc1,$sc2)=(" "," ");
			foreach(@case){
				next if(not exists $mut{$_});
				my $mutN=(split /:/,$mut{$_})[0];
				$hc0++ if($mutN==0);
				$hc1++ if($mutN==1);
				$hc2++ if($mutN==2);
				$sc1.="$_," if($mutN==1);
				$sc2.="$_," if($mutN==2);
			}
			my ($hn0,$hn1,$hn2)=(0,0,0);
			my ($sn1,$sn2)=(" "," ");
			foreach(@control){
				next if(not exists $mut{$_});
				my $mutN=(split /:/,$mut{$_})[0];
				$hn0++ if($mutN==0);
				$hn1++ if($mutN==1);
				$hn2++ if($mutN==2);
				$sn1.="$_," if($mutN==1);
				$sn2.="$_," if($mutN==2);
			}
			my ($h0,$h1,$h2)=($hc0+$hn0,$hc1+$hn1,$hc2+$hn2);
			next if($h0+$h1+$h2==0);
			my $freq=0;$freq=($h1/2+$h2)/($h0+$h1+$h2) if($h0+$h1+$h2>0);
			my $contperc=0;$contperc=($hn1/2+$hn2)/($hn0+$hn1+$hn2) if($hn0+$hn1+$hn2>0);
			my $Gene=$library->{$title}{$alt}{"Gene"};
			my $blast=$library->{$title}{$alt}{"HOMOLOGY HITS"};
			my $snpcq=$library->{$title}{$alt}{"SNP Calling Quality"};
			my $Chrs=$library->{$title}{$alt}{"Chrs"};
			my $position2=$library->{$title}{$alt}{"Position2"};			
			my $regionfiltermark=0;###区间过滤标志，默认不过滤
			$regionfiltermark=CheckRegion($Chrs,$position2,\%regionfilter) if(exists $run->{"Region"}{1});			
			$snpcq-- if($blast>1);
			my $snvq="H";$snvq="M" if($snpcq==2);$snvq="L" if($snpcq<2);
			my $genoq=$analysis->{$title}{$alt}{'Genotyping Quality'};
			my $Function=$library->{$title}{$alt}{"Function"};
			my $g1000=0;$g1000=$library->{$title}{$alt}{"Freq_Alt (1000g)"} if($library->{$title}{$alt}{"Freq_Alt (1000g)"}=~ /\d/);
			my $ExAC03=0;$ExAC03=$annotation->{$title}{$alt}{"ExAC03"} if(exists $annotation->{$title}{$alt}{"ExAC03"});
			my $ExAC03_EAS=0;$ExAC03_EAS=$annotation->{$title}{$alt}{"ExAC03_EAS"} if(exists $annotation->{$title}{$alt}{"ExAC03_EAS"});
			my $gnomAD_EAS=0;$gnomAD_EAS=$annotation->{$title}{$alt}{"gnomAD_exome_EAS"} if(exists $annotation->{$title}{$alt}{"gnomAD_exome_EAS"});
			my $hgmd=0;$hgmd=1 if(exists $annotation->{$title}{$alt}{"HGMD_site_class"} and ($annotation->{$title}{$alt}{"HGMD_site_class"}=~ /^\s*DM\s*$/ or ($g1000<0.05 and $annotation->{$title}{$alt}{"HGMD_site_class"}=~ /^\s*DM\?\s*$/)));
			my $geneskyHit=0;$geneskyHit=$annotation->{$title}{$alt}{"GeneskyExonDB_SampeCount"} if(exists $annotation->{$title}{$alt}{"GeneskyExonDB_SampeCount"});
			my $geneskyHitFreq=0;$geneskyHitFreq=$annotation->{$title}{$alt}{"GeneskyExonDB_Freq"} if(exists $annotation->{$title}{$alt}{"GeneskyExonDB_Freq"});
			my $esp6500=0;$esp6500=$annotation->{$title}{$alt}{'esp6500'} if(exists $annotation->{$title}{$alt}{'esp6500'} and $annotation->{$title}{$alt}{'esp6500'}=~ /\d/);
			my $bj301Control=0;$bj301Control=$annotation->{$title}{$alt}{"bj301 ControlDB"} if(exists $annotation->{$title}{$alt}{"bj301 ControlDB"});
			my $sift=0;
			$sift=1 if(exists $annotation->{$title}{$alt}{'SIFT Score Pred'} and $annotation->{$title}{$alt}{'SIFT Score Pred'}=~/D/);
			$sift=1 if(exists $annotation->{$title}{$alt}{'POLYPhen V2 Score Pred'} and $annotation->{$title}{$alt}{'POLYPhen V2 Score Pred'}=~ /D/);
			$sift=1 if(exists $annotation->{$title}{$alt}{'MutationTaster Pred'} and $annotation->{$title}{$alt}{'MutationTaster Pred'}=~ /D/);			
			$sift=1 if(exists $annotation->{$title}{$alt}{'Cadd_Raw'} and $annotation->{$title}{$alt}{'Cadd_Raw'}>=4);
			$sift=1 if(exists $annotation->{$title}{$alt}{'Dann'} and $annotation->{$title}{$alt}{'Dann'}>=0.93);
			$sift=1 if( (exists $annotation->{$title}{$alt}{'dbscSNV_ADA_SCORE'} and $annotation->{$title}{$alt}{'dbscSNV_ADA_SCORE'}>=0.6) or (exists $annotation->{$title}{$alt}{'dbscSNV_RF_SCORE'} and $annotation->{$title}{$alt}{'dbscSNV_RF_SCORE'}>=0.6)); # 剪切位点
			$sift=1 if(($library->{$title}{$alt}{'Ref Allele'}=~ /\-/ or $alt=~ /\-/) and $library->{$title}{$alt}{"Gene Region"} =~/exonic/);
			$sift=1 if(($library->{$title}{$alt}{'Gene Region'} eq 'splicing' or $library->{$title}{$alt}{'Gene Region'} eq 'exonic;splicing') and $library->{$title}{$alt}{'Predicted Protein Variants'}=~ /[\+\-][1-2][ATCG-]/);
			my $Priority=package::priority::Site($hgmd,$sift,$g1000,$contperc,$esp6500,$snvq,$genoq,$blast,$geneskyHitFreq);
			# filter
			next if(exists $run->{"Rm_Function"}{1} and $run->{"Rm_Function"}{1}=~/^$Function$/i and $hgmd==0);# 移出指定功能位点，通常为UNKNOWN
			next if(exists $run->{"Region"}{1} and $regionfiltermark==1);
			next if(exists $run->{"SNP Calling Quality"}{1} and $run->{"SNP Calling Quality"}{1}=~/$snvq/);
            next if(exists $run->{"Genotyping Quality"}{1} and $genoq<$run->{"Genotyping Quality"}{1});	
			next if(exists $run->{"Freq_Alt (1000g)"}{1} and $g1000>$run->{"Freq_Alt (1000g)"}{1} and $hgmd==0);
			next if(exists $run->{"ExAC03"}{1} and $ExAC03>$run->{"ExAC03"}{1} and $hgmd==0);
			next if(exists $run->{"ExAC03_EAS"}{1} and $ExAC03_EAS=~/\d/ and $ExAC03_EAS>$run->{"ExAC03_EAS"}{1} and $hgmd==0); # ExAC03_EAS有为"."的情况
			next if(exists $run->{"gnomAD_exome_EAS"}{1} and $gnomAD_EAS=~/\d/ and $gnomAD_EAS>$run->{"gnomAD_exome_EAS"}{1} and $hgmd==0);
			next if(exists $run->{"esp6500"}{1} and $esp6500>$run->{"esp6500"}{1} and $hgmd==0);
			next if(exists $run->{"Control Percentage"}{1} and $contperc>$run->{"Control Percentage"}{1} and $hgmd==0);
			next if(exists $run->{"Alt Allele Freq"}{1} and $freq>$run->{"Alt Allele Freq"}{1} and $hgmd==0);
			next if(exists $run->{"GeneskyExonDB_Freq"}{1} and $geneskyHit>20 and $geneskyHitFreq>$run->{"GeneskyExonDB_Freq"}{1});
			next if(exists $run->{"Chrs"}{1} and !exists($chrfilter{$Chrs}) );###选取特定染色体
			next if(exists $run->{"Gene"}{1} and !exists($genefilter{$Gene}));###选取特定基因
			next if(exists $run->{"bj301 ControlDB"}{1} and $bj301Control>$run->{"bj301 ControlDB"}{1} and $hgmd==0);
			next if(exists $run->{"bj301 GeneInfo"}{1} and not exists $annotation->{$title}{$alt}{"bj301 GeneInfo"});
			next if(exists $run->{"OMIM"}{1} and not exists $annotation->{$title}{$alt}{"OMIM"});
			next if(exists $run->{"Priority"}{1} and $run->{"Priority"}{1}!~/$Priority/); # 保留指定优先级位点
			next if(exists $run->{"Danger_Predict"}{1} and $run->{"Danger_Predict"}{1}!=$sift); # 保留突变危险性预测。突变危险性预测包含危害性标注，包括剪切位点也要软件预测致病，或者cds区域的插入缺失
			next if(exists $run->{"VCF Filter"}{1} and exists $annotation->{$title}{$alt}{"VCF Filter"} and $run->{"VCF Filter"}{1}!~/$annotation->{$title}{$alt}{"VCF Filter"}/); # 保留VCF Filter信息满足条件的
			# filter 2
			my ($good,$all)=(0,0);
			if(exists $run->{"AND_Site"}){
				foreach my $weiyi(keys %{$run->{"AND_Site"}}){
					my $isok=1;
					my @info=split /,/,$run->{"AND_Site"}{$weiyi};
					for(my $i=0;$i<@info-1;$i+=2){
						$isok=0 if(not exists $mut{$info[$i]} or (split /:/,$mut{$info[$i]})[0]!=$info[$i+1]);
					}
					$good++ if($isok==1);
					$all++;
				}
			}
			if(exists $run->{"OR_Site"}){
				foreach my $weiyi(keys %{$run->{"OR_Site"}}){
					my $isok=0;
					my @info=split /,/,$run->{"OR_Site"}{$weiyi};
					for(my $i=0;$i<@info-1;$i+=2){
						$isok=1 if(exists $mut{$info[$i]} and (split /:/,$mut{$info[$i]})[0]==$info[$i+1]);
					}
					$good++ if($isok==1);
					$all++;
				}
			}
			if(exists $run->{"NOT_Site"}){
				foreach my $weiyi(keys %{$run->{"NOT_Site"}}){
					my $isok=0;
					my @info=split /,/,$run->{"NOT_Site"}{$weiyi};
					for(my $i=0;$i<@info-1;$i+=2){
						$isok=1 if(exists $mut{$info[$i]} and (split /:/,$mut{$info[$i]})[0]!=$info[$i+1]);
					}
					$good++ if($isok==1);
					$all++;
				}
			}
			next if(exists $run->{"MOD_SET"}{1} and $run->{"MOD_SET"}{1} eq "ALL" and $all>0 and $good!=$all);
			next if(exists $run->{"MOD_SET"}{1} and $run->{"MOD_SET"}{1} eq "ONE" and $all>0 and $good==0);
			
			# # 频率限制模式（对有深度信息的样本进行过滤）
			my $is_freq_ok = 1;
			if(exists $run->{"Freq_Limit_Model"}{1} and $run->{"Freq_Limit_Model"}{1} == 1) { # # 频率限制模式下对AND_Site的考察
				die "[ERR] Freq_Limit_Model without AND_Site\n" if(not exists $run->{"AND_Site"}); # # Freq_Limit_Model配置从AND_Site读取
				foreach my $weiyi(keys %{$run->{"AND_Site"}}) {
					my @info = split /,/, $run->{"AND_Site"}{$weiyi};
					for(my $i = 0; $i < @info - 1; $i += 2) {
						my $sample          = $info[$i];
						my $require_snv_num = $info[$i + 1];
						
						if(exists $mut{$sample}) {
						    my ($real_snv_num, $tmp1, $tmp2, $tmp3, $ref_depth, $alt_depth) = split /:/, $mut{$sample};
                            $ref_depth = 0 if($ref_depth !~ /\d/);
                            $alt_depth = 0 if($alt_depth !~ /\d/);
						    my $depth = $ref_depth + $alt_depth;
						    if($require_snv_num > 0) { # # 要求的突变数大于0，则视为case
						    	if($depth < 10 or $alt_depth < 5) { # # case最低测序深度10X,突变等位基因reads最低5X
						    	    $is_freq_ok = 0;
						    		next;
						    	}
						    	$is_freq_ok = 0 if(($alt_depth > 10) and (($alt_depth / $depth) < 0.2));   # # 如果突变深度大于10X，则突变reads比例>=0.2
						    	$is_freq_ok = 0 if(($alt_depth <= 10) and (($alt_depth / $depth) < 0.28)); # # 如果突变深度小于10X，则>=0.28
						    }
						    if($require_snv_num == 0) { # # 要求的突变数为0，则视为control
						    	$is_freq_ok = 0 if(($depth <= 10) or (($alt_depth / $depth) >= 0.035));    # # control的最低深10X,突变等位基因比小于3.5%
						    }
						}
					}
				}
			}
			next if($is_freq_ok == 0);
			
			# 记录
			$gene{$Gene}{$title}{$alt}{"SNP Calling Quality"}=$snvq;
			$gene{$Gene}{$title}{$alt}{"First Priority"}=$Priority;
			$gene{$Gene}{$title}{$alt}{"Control Percentage"}=$contperc;
			$gene{$Gene}{$title}{$alt}{"Case(0|1|2)"}="$hc0|$hc1|$hc2";
			$gene{$Gene}{$title}{$alt}{"Control(0|1|2)"}="$hn0|$hn1|$hn2";
			$gene{$Gene}{$title}{$alt}{"Alt Allele Freq"}=($h1/2+$h2)/($h0+$h1+$h2);
			$gene{$Gene}{$title}{$alt}{"Case_het"}=$sc1;
			$gene{$Gene}{$title}{$alt}{"Case_hom"}=$sc2;
			$gene{$Gene}{$title}{$alt}{"Control_het"}=$sn1;
			$gene{$Gene}{$title}{$alt}{"Control_hom"}=$sn2;
			$gene{$Gene}{$title}{$alt}{"gene"}=$Gene;###只为Union.pl存在
			foreach my $sample(keys %mut){
				my @splits=split /:/,$mut{$sample};
				$geneMut{$Gene}{$sample}+=$splits[0];
				$geneMut{$Gene}{$sample}=2 if($geneMut{$Gene}{$sample}>2);
				$gene{$Gene}{$title}{$alt}{$sample."_geno"}=$splits[3];
				$gene{$Gene}{$title}{$alt}{$sample."_reads"}="$splits[4]:$splits[5]";
				$geneAutoComHet{$Gene}{$sample}++ if($splits[0]==1);
			}
		}
	}
	print "\n";
	my $id=1;
	my %filter=();
	foreach my $Gene(keys %gene){
		my ($good,$all)=(0,0);
		if(exists $run->{"AND_Gene"}){
			foreach my $weiyi(keys %{$run->{"AND_Gene"}}){
				my $isok=1;
				my @info=split /,/,$run->{"AND_Gene"}{$weiyi};
				for(my $i=0;$i<@info-1;$i+=2){
					$isok=0 if(not exists $geneMut{$Gene}{$info[$i]} or $geneMut{$Gene}{$info[$i]}!=$info[$i+1]);
				}
				$good++ if($isok==1);
				$all++;
			}
		}
		if(exists $run->{"OR_Gene"}){
			foreach my $weiyi(keys %{$run->{"OR_Gene"}}){
				my $isok=0;
				my @info=split /,/,$run->{"OR_Gene"}{$weiyi};
				for(my $i=0;$i<@info-1;$i+=2){
					$isok=1 if(exists $geneMut{$Gene}{$info[$i]} and $geneMut{$Gene}{$info[$i]}==$info[$i+1]);
				}
				$good++ if($isok==1);
				$all++;
			}
		}
		if(exists $run->{"NOT_Gene"}){
			foreach my $weiyi(keys %{$run->{"NOT_Gene"}}){
				my $isok=0;
				my @info=split /,/,$run->{"NOT_Gene"}{$weiyi};
				for(my $i=0;$i<@info-1;$i+=2){
					$isok=1 if(exists $geneMut{$Gene}{$info[$i]} and $geneMut{$Gene}{$info[$i]}!=$info[$i+1]);
				}
				$good++ if($isok==1);
				$all++;
			}
		}
		next if(exists $run->{"MOD_SET"}{1} and $run->{"MOD_SET"}{1} eq "ALL" and $all>0 and $good!=$all);
		next if(exists $run->{"MOD_SET"}{1} and $run->{"MOD_SET"}{1} eq "ONE" and $all>0 and $good==0);
		if(exists $run->{"Model"}{1} and $run->{"Model"}{1} eq "AutoComHet"){
			my $isAutoComHet=0;
			foreach my $sample(keys %{$geneAutoComHet{$Gene}}){
				$isAutoComHet=1 if($geneAutoComHet{$Gene}{$sample}>=2);
			}
			next if($isAutoComHet==0);
		}
		foreach my $title(sort{$a cmp $b}keys %{$gene{$Gene}}){
			my $weiyi=1;
			foreach my $alt(keys %{$gene{$Gene}{$title}}){
				my $subID="";$subID=".".$weiyi if(keys %{$gene{$Gene}{$title}}>1);
				$filter{$id}{$title}{$alt}{"SNV NO."}="SNV"."0"x(5-length($id)).$id;
				foreach my $userInfo(keys %{$gene{$Gene}{$title}{$alt}}){
					$filter{$id}{$title}{$alt}{$userInfo}=$gene{$Gene}{$title}{$alt}{$userInfo};
				}
				$weiyi++;
				$id++;
			}
		}
	}
	return %filter;
}
sub RunRegion{
    my $run=shift @_;
	my %regionfilter=();
	if(exists $run->{"Region"}{1}){
	   my @regions=split /,/,$run->{"Region"}{1};
	   foreach(@regions){
	       my ($chr,$start,$end)=split /-/,$_;
		   $regionfilter{$chr}{'start'}=$start;
		   $regionfilter{$chr}{'end'}=$end;
	   }
	}
    return %regionfilter;
}

sub CheckRegion{
    my $chr=shift @_;
	my $position=shift @_;
	my $regionfilter=shift @_;
	my $mark=0;
	if(exists($regionfilter->{$chr}) and $position>=$regionfilter->{$chr}{'start'} and $position<=$regionfilter->{$chr}{'end'}){
	}else{
	   $mark=1;
	}
	return $mark;
}

sub RunFilterTarget{
    my $run=shift @_;
	my $type=shift @_;
	my %thisfilter=();
	if(exists $run->{$type}{1}){
	   my @values=split /,/,$run->{$type}{1};
	   map{$thisfilter{$_}++}@values;
	}
    return %thisfilter;
}

1
