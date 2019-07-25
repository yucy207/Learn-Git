$|=1;
use strict;
use warnings;
use Excel::Writer::XLSX;
use package::main;
use package::time;
use package::format;
use package::priority;
use package::readme;
use package::Math::p1p2p3p4;
use package::Math::logistic;
use package::Math::hwe;
use package::Math::leastSquare;
use package::Math::determinant;
use package::Math::basic;
use package::Math::normDist;
my $tab="##########";

# 读取配置
die "Usage perl $0 ConfigFile [outputName]\n" if(@ARGV!=1 and  @ARGV!=2);
my $configFile=shift @ARGV;
my $outputName = (exists($ARGV[0])) ? $ARGV[0] : "allSNV"; 
my %config=package::main::readConfig($configFile);
my @case=package::main::Sample(\%config,"case");
my @control=package::main::Sample(\%config,"control");
my @sample=(@case,@control);
# 开始流程
system "clear";
print "$tab Start All SNV Pipline ".package::time::ymd." ".package::time::hms." $tab\n";
die "[ERR] Loss Report\n" if(not exists $config{"Report"} or not -e $config{"Report"});
die "[ERR] Loss analysis\n" if(not -e $config{"Report"}."/analysis_all.txt");
die "[ERR] Loss annotation\n" if(not -e $config{"Report"}."/annotation.txt");
die "[ERR] Loss snv\n" if(not -e $config{"Report"}."/database.snv");
die "[ERR] Loss read me\n" if(not -e "readme_snvall.txt");
# read
my %library=package::main::readLibrary($config{"Report"}."/database.snv");
my %annotation=package::main::readAnnotation($config{"Report"}."/annotation.txt");
my %analysis=package::main::readAnalysis($config{"Report"}."/analysis_all.txt");
# Excel
mkdir $config{"Report"}."/document" if(not -e $config{"Report"}."/document");
my $workbook=Excel::Writer::XLSX->new($config{"Report"}."/document/$outputName.xlsx");
my %format=package::format::run($workbook);
my $SNV=$workbook->add_worksheet("SNV");
my $GENOTYPES=$workbook->add_worksheet("GENOTYPES");
my $High=$workbook->add_worksheet("High Freq. SNV");
my $Readme=$workbook->add_worksheet("READ ME");
$SNV->set_row(0, 60);
$GENOTYPES->set_row(0, 60);
$Readme->set_row(0, 60);
my @title=(
	"SNV NO.",
	"Gene",
	"First Priority",
	"SNP Calling Quality",
	"Genotyping Quality",
	"VCF Filter",
	"Case_het",
	"Case_hom",
	"Control_het",
	"Control_hom",	
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
	"Control(0|1|2)",
	"Alt Allele Freq",
	"Control Percentage",
);
my @title2=(
	"SNV NO.",
	"Gene",
	"First Priority",
	"SNP ID",
	"Freq_Alt (1000g)",
	"Chrs",
	"Position",
	"Ref Allele",
	"Alt Allele",	
	"Gene Region",
	"Function",
	"Case(0|1|2)",
	"Control(0|1|2)",
	"Alt Allele Freq",
	"HWE Case",
	"HWE Control",
	"HWE",
	"p1",
	"p2",
	"p3",
	"p4",
);
my $row=0;
for my $i(0..@title-1){
	$SNV->write($row,$i,$title[$i],$format{"title"});
	$GENOTYPES->write($row,$i,$title[$i],$format{"title"});
}
my $row_high=0;
for my $i(0..@title2-1){
	$High->merge_range($row_high,$i,$row_high+1,$i,$title2[$i],$format{"title"});
}
my @mode=("Additive","Dominant","Recessive","Allele");
for my $i(0..@mode-1){
	$High->merge_range($row_high,@title2+$i*2,$row_high,@title2+$i*2+1,$mode[$i],$format{"title"});
	$High->write($row_high+1,@title2+$i*2,"OR(95%CI)",$format{"title"});
	$High->write($row_high+1,@title2+$i*2+1,"P-value",$format{"title"});
}
for my $i(0..@sample-1){
	$SNV->write_string($row,@title+$i,$sample[$i],$format{"title"});
	$GENOTYPES->write_string($row,@title+$i,$sample[$i],$format{"title"});
	$High->merge_range($row_high,@title2+@mode*2+$i,$row_high+1,@title2+@mode*2+$i,$sample[$i],$format{"title"});
}
$row++;
$row_high+=2;
my $id=1;
my $count=0;
my @snvs=sort keys %library;
my $snvnum=@snvs;
foreach my $title(@snvs){
	$count++;
	my $p=sprintf "%0.4f",$count/$snvnum;
	print "\r".$config{"Report"}."/document/allSNV.xlsx ...$p";
	foreach my $alt(sort keys %{$library{$title}}){	
		my %mut=();
		foreach(split /,/,$analysis{$title}{$alt}{"Mutation 0"}){$mut{(split /:/,$_)[0]}="0:$_" if($_=~ /\w/);}
		foreach(split /,/,$analysis{$title}{$alt}{"Mutation 1"}){$mut{(split /:/,$_)[0]}="1:$_" if($_=~ /\w/);}
		foreach(split /,/,$analysis{$title}{$alt}{"Mutation 2"}){$mut{(split /:/,$_)[0]}="2:$_" if($_=~ /\w/);}
		my ($hc0,$hc1,$hc2)=(0,0,0);
		my ($sc1,$sc2)=("","");
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
		my ($sn1,$sn2)=("","");
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
		next if(($h1+$h2)==0);##跳过没有突变的位点
		next if(($h0+$h1+$h2)==0);##跳过没有信息的位点		
		my $SNVID="SNV"."0"x(5-length($id)).$id;
		$id++;
		my $altFreq=0;$altFreq=($h1/2+$h2)/($h0+$h1+$h2) if($h0+$h1+$h2>0);
		my $contperc=0;$contperc=($hn1/2+$hn2)/($hn0+$hn1+$hn2) if($hn0+$hn1+$hn2>0);
		my $Gene=$library{$title}{$alt}{"Gene"};
		my $blast=$library{$title}{$alt}{"HOMOLOGY HITS"};
		my $snpcq=$library{$title}{$alt}{"SNP Calling Quality"};
		$snpcq-- if($blast>1);
		my $snvq="H";$snvq="M" if($snpcq==2);$snvq="L" if($snpcq<2);
		my $genoq=$analysis{$title}{$alt}{'Genotyping Quality'};
		my $g1000=0;$g1000=$library{$title}{$alt}{"Freq_Alt (1000g)"} if($library{$title}{$alt}{"Freq_Alt (1000g)"}=~ /\d/);
		my $hgmd=0;$hgmd=1 if(exists $annotation{$title}{$alt}{"HGMD_site_class"} and ($annotation{$title}{$alt}{"HGMD_site_class"}=~ /^\s*DM\s*$/ or ($g1000<0.05 and $annotation{$title}{$alt}{"HGMD_site_class"}=~ /^\s*DM\?\s*$/)));
		my $geneskyHit=0;$geneskyHit=$annotation{$title}{$alt}{"GeneskyExonDB_SampeCount"} if(exists $annotation{$title}{$alt}{"GeneskyExonDB_SampeCount"});
		my $esp6500=0;$esp6500=$annotation{$title}{$alt}{'esp6500'} if(exists $annotation{$title}{$alt}{'esp6500'} and $annotation{$title}{$alt}{'esp6500'}=~ /\d/);
		my $sift=0;
		$sift=1 if(exists $annotation{$title}{$alt}{'SIFT Score Pred'} and $annotation{$title}{$alt}{'SIFT Score Pred'}=~/D/);
		$sift=1 if(exists $annotation{$title}{$alt}{'POLYPhen V2 Score Pred'} and $annotation{$title}{$alt}{'POLYPhen V2 Score Pred'}=~ /D/);
		$sift=1 if(exists $annotation{$title}{$alt}{'MutationTaster Pred'} and $annotation{$title}{$alt}{'MutationTaster Pred'}=~ /D/);
		$sift=1 if(exists $annotation{$title}{$alt}{'Cadd'} and $annotation{$title}{$alt}{'Cadd'}>=4);
		$sift=1 if(exists $annotation{$title}{$alt}{'Dann'} and $annotation{$title}{$alt}{'Dann'}>=0.93);
		$sift=1 if( (exists $annotation{$title}{$alt}{'dbscSNV_ADA_SCORE'} and $annotation{$title}{$alt}{'dbscSNV_ADA_SCORE'}>=0.6) or (exists $annotation{$title}{$alt}{'dbscSNV_RF_SCORE'} and $annotation{$title}{$alt}{'dbscSNV_RF_SCORE'}>=0.6)); # 剪切位点
		$sift=1 if(($library{$title}{$alt}{'Ref Allele'}=~ /\-/ or $alt=~ /\-/) and $library{$title}{$alt}{"Gene Region"} =~/exonic/);
		$sift=1 if(($library{$title}{$alt}{'Gene Region'} eq 'splicing' or $library{$title}{$alt}{'Gene Region'} eq 'exonic;splicing') and $library{$title}{$alt}{'Predicted Protein Variants'}=~ /[\+\-][1-2][ATCG-]/);
		my $Priority=package::priority::Site($hgmd,$sift,$g1000,$contperc,$esp6500,$snvq,$genoq,$blast,$geneskyHit);
		for my $i(0..@title-1){
			my $v=" ";
			$v=$library{$title}{$alt}{$title[$i]} if(exists $library{$title}{$alt}{$title[$i]});
			$v=$snvq if($title[$i] eq 'SNP Calling Quality');
			$v=$annotation{$title}{$alt}{$title[$i]} if(exists $annotation{$title}{$alt}{$title[$i]});
			$v=$analysis{$title}{$alt}{$title[$i]} if(exists $analysis{$title}{$alt}{$title[$i]});
			$v=$SNVID if($title[$i] eq "SNV NO.");
			$v=$Priority if($title[$i] eq "First Priority");
			$v="$hc0|$hc1|$hc2" if($title[$i] eq "Case(0|1|2)");
			$v="$hn0|$hn1|$hn2" if($title[$i] eq "Control(0|1|2)");
			$v=$altFreq if($title[$i] eq "Alt Allele Freq");
			$v=$sc1 if($title[$i] eq "Case_het");
			$v=$sc2 if($title[$i] eq "Case_hom");
			$v=$sn1 if($title[$i] eq "Control_het");
			$v=$sn2 if($title[$i] eq "Control_hom");
			$v=$contperc if($title[$i] eq "Control Percentage");
			$SNV->write($row,$i,$v,$format{"normal"});
			$GENOTYPES->write($row,$i,$v,$format{"normal"});
		}
		for my $i(0..@sample-1){
			my ($geno,$reads)=(" "," ");
			if(exists $mut{$sample[$i]}){
				my @splits=split /:/,$mut{$sample[$i]};
				$geno="$splits[3]";
				$reads="$splits[2],".$library{$title}{$alt}{'Ref Allele'}.":$splits[4],$alt:$splits[5]";
			}
			$SNV->write($row,@title+$i,$reads,$format{"left"});
			$GENOTYPES->write($row,@title+$i,$geno,$format{"normal"});
		}
		$row++;
		if($g1000>=0.05){
			my ($p1,$p2,$p3,$p4)=package::Math::p1p2p3p4::run($hc0,$hc1,$hc2,$hn0,$hn1,$hn2);
			for my $i(0..@title2-1){
				my $v=" ";
				$v=$library{$title}{$alt}{$title2[$i]} if(exists $library{$title}{$alt}{$title2[$i]});
				$v=$snvq if($title2[$i] eq 'SNP Calling Quality');
				$v=$annotation{$title}{$alt}{$title2[$i]} if(exists $annotation{$title}{$alt}{$title2[$i]});
				$v=$SNVID if($title2[$i] eq "SNV NO.");
				$v=$Priority if($title2[$i] eq "First Priority");
				$v="$hc0|$hc1|$hc2" if($title2[$i] eq "Case(0|1|2)");
				$v="$hn0|$hn1|$hn2" if($title2[$i] eq "Control(0|1|2)");
				$v=$altFreq if($title2[$i] eq "Alt Allele Freq");
				$v=package::Math::hwe::run($hc1,$hc0,$hc2) if($title2[$i] eq "HWE Case");
				$v=package::Math::hwe::run($hn1,$hn0,$hn2) if($title2[$i] eq "HWE Control");
				$v=package::Math::hwe::run($h1,$h0,$h2) if($title2[$i] eq "HWE");
				$v=$p1 if($title2[$i] eq "p1");
				$v=$p2 if($title2[$i] eq "p2");
				$v=$p3 if($title2[$i] eq "p3");
				$v=$p4 if($title2[$i] eq "p4");
				$High->write($row_high,$i,$v,$format{"normal"});
			}
			my %data=();
			# Additive
			%{$data{"Additive"}{1}}=("X1"=>0,"DV"=>1,"REP"=>$hc0);
			%{$data{"Additive"}{2}}=("X1"=>1,"DV"=>1,"REP"=>$hc1);
			%{$data{"Additive"}{3}}=("X1"=>2,"DV"=>1,"REP"=>$hc2);
			%{$data{"Additive"}{4}}=("X1"=>0,"DV"=>0,"REP"=>$hn0);
			%{$data{"Additive"}{5}}=("X1"=>1,"DV"=>0,"REP"=>$hn1);
			%{$data{"Additive"}{6}}=("X1"=>2,"DV"=>0,"REP"=>$hn2);
			# Dominant
			%{$data{"Dominant"}{1}}=("X1"=>0,"DV"=>1,"REP"=>$hc0);
			%{$data{"Dominant"}{2}}=("X1"=>1,"DV"=>1,"REP"=>$hc1+$hc2);
			%{$data{"Dominant"}{3}}=("X1"=>0,"DV"=>0,"REP"=>$hn0);
			%{$data{"Dominant"}{4}}=("X1"=>1,"DV"=>0,"REP"=>$hn1+$hn2);
			# Recessive
			%{$data{"Recessive"}{1}}=("X1"=>0,"DV"=>1,"REP"=>$hc0+$hc1);
			%{$data{"Recessive"}{2}}=("X1"=>1,"DV"=>1,"REP"=>$hc2);
			%{$data{"Recessive"}{3}}=("X1"=>0,"DV"=>0,"REP"=>$hn0+$hn1);
			%{$data{"Recessive"}{4}}=("X1"=>1,"DV"=>0,"REP"=>$hn2);
			# Allele
			%{$data{"Allele"}{1}}=("X1"=>0,"DV"=>1,"REP"=>$hc0*2);
			%{$data{"Allele"}{2}}=("X1"=>0,"DV"=>1,"REP"=>$hc1);
			%{$data{"Allele"}{3}}=("X1"=>1,"DV"=>1,"REP"=>$hc1);
			%{$data{"Allele"}{4}}=("X1"=>1,"DV"=>1,"REP"=>$hc2*2);
			%{$data{"Allele"}{5}}=("X1"=>0,"DV"=>0,"REP"=>$hn0*2);
			%{$data{"Allele"}{6}}=("X1"=>0,"DV"=>0,"REP"=>$hn1);
			%{$data{"Allele"}{7}}=("X1"=>1,"DV"=>0,"REP"=>$hn1);
			%{$data{"Allele"}{8}}=("X1"=>1,"DV"=>0,"REP"=>$hn2*2);
			for my $i(0..@mode-1){
				my %OR=();%OR=package::Math::logistic::run($data{$mode[$i]}) if(exists $data{$mode[$i]});
				if(exists $OR{"X1"}{"OR"} and $OR{"X1"}{"OR"}=~ /\d/){
					my ($or,$or_low,$or_up,$or_p)=($OR{"X1"}{"OR"},$OR{"X1"}{"OR95CI_Low"},$OR{"X1"}{"OR95CI_Up"},$OR{"X1"}{"Pvalue"});
					$or=sprintf("%.3f",$or);$or_low=sprintf("%.3f",$or_low);$or_up=sprintf("%.3f",$or_up);
					$High->write($row_high,@title2+$i*2,"$or($or_low-$or_up)",$format{"left"});
					$High->write($row_high,@title2+$i*2+1,$or_p,$format{"left"});
				}else{
					$High->write($row_high,@title2+$i*2,"unKnown",$format{"left"});
					$High->write($row_high,@title2+$i*2+1,"unKnown",$format{"left"});
				}
			}
			for my $i(0..@sample-1){
				my $geno=" ";
				if(exists $mut{$sample[$i]}){
					my @splits=split /:/,$mut{$sample[$i]};
					$geno="$splits[3]";
				}
				$High->write($row_high,@title2+@mode*2+$i,$geno,$format{"normal"});
			}
			$row_high++;
		}
	}
}
print "   OK\n";
package::readme::run($Readme,\%format,"readme_snvall.txt");

