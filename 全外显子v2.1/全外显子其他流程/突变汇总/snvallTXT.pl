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

# ��ȡ����
die "Usage perl $0 ConfigFile\n" if(@ARGV!=1);
my $configFile=shift @ARGV;
my %config=package::main::readConfig($configFile);
my @case=package::main::Sample(\%config,"case");
my @control=package::main::Sample(\%config,"control");
my @sample=(@case,@control);
# ��ʼ����
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
print $config{"Report"}."/document/allSNV.txt\n";
open OUT,">$config{'Report'}/document/allSNV.txt";
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
for my $i(0..@title-1){
	print OUT "$title[$i]\t";
}
for my $i(0..@sample-1){
    if($i==$#sample){
	   print OUT "$sample[$i]\n";
	}else{
	   print OUT "$sample[$i]\t";
	}    
}
my $id=1;
foreach my $title(sort keys %library){
	my $weiyi=1;
	foreach my $alt(sort keys %{$library{$title}}){
		my $subID="";$subID=".".$weiyi if(keys %{$library{$title}}>1);
		$weiyi++;
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
		next if(($h1+$h2)==0);##����û��ͻ���λ��
		next if(($h0+$h1+$h2)==0);##����û����Ϣ��λ��
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
		$sift=1 if( (exists $annotation{$title}{$alt}{'dbscSNV_ADA_SCORE'} and $annotation{$title}{$alt}{'dbscSNV_ADA_SCORE'}>=0.6) or (exists $annotation{$title}{$alt}{'dbscSNV_RF_SCORE'} and $annotation{$title}{$alt}{'dbscSNV_RF_SCORE'}>=0.6)); # ����λ��
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
			print OUT "$v\t";
		}
		for my $i(0..@sample-1){
			my ($geno,$reads)=(" "," ");
			if(exists $mut{$sample[$i]}){
				my @splits=split /:/,$mut{$sample[$i]};
				$geno="$splits[3]";
				$reads="$splits[2],".$library{$title}{$alt}{'Ref Allele'}.":$splits[4],$alt:$splits[5]";
			}
			if($i==$#sample){
			   print OUT "$geno\n";
			}else{
			   print OUT "$geno\t";
			}			
		}
	}
}

