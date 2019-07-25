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
die "Usage perl $0 ConfigFile\n" if(@ARGV!=1);
my $configFile=shift @ARGV;
my %config=package::main::readConfig($configFile);
my @case=package::main::Sample(\%config,"case");
my @control=package::main::Sample(\%config,"control");
my @sample=(@case,@control);
my $min_var_freq  = (exists $config{'MinVarFreq'})? $config{'MinVarFreq'} : 0.1;
# 开始流程
system "clear";
print "$tab Start All SNV Pipline ".package::time::ymd." ".package::time::hms." $tab\n";
die "[ERR] Loss Report\n" if(not exists $config{"Report"} or not -e $config{"Report"});
die "[ERR] Loss analysis\n" if(not -e $config{"Report"}."/analysis_all.txt");
die "[ERR] Loss annotation\n" if(not -e $config{"Report"}."/annotation.txt");
die "[ERR] Loss snv\n" if(not -e $config{"Report"}."/database.snv");
die "[ERR] Loss read me\n" if(not -e "readme_snvallMT.txt");
# read
my %library=package::main::readLibrary($config{"Report"}."/database.snv");
my %annotation=package::main::readAnnotation($config{"Report"}."/annotation.txt");
my %analysis=package::main::readAnalysis($config{"Report"}."/analysis_all.txt");
# Excel
mkdir $config{"Report"}."/document" if(not -e $config{"Report"}."/document");
print $config{"Report"}."/document/MT_allSNV.xlsx\n";
my $workbook=Excel::Writer::XLSX->new($config{"Report"}."/document/MT_allSNV.xlsx");
my %format=package::format::run($workbook);
my $SNV=$workbook->add_worksheet("SNV");
my $Readme=$workbook->add_worksheet("READ ME");
$SNV->set_row(0, 60);
$Readme->set_row(0, 60);
my @title=(
	"SNV NO.",
	"Gene",
	"Chrs",
	"Position",
	"Ref Allele",
	"Alt Allele",	
	"Gene Strand Orientation",
	"Gene Region",
    "Function",
	"Predicted Protein Variants",
	"Codon Mutation",
	"Mutation Pattern",
	"tRNA_anticodon",
	"phastCons100way",
	"phyloP100way",
	"HmtVar_Asia_healthy_freq",
	"HmtVar_Asia_pathologic_freq",
	"MT_Control_Freq",
	"mitomap_disease",
	"PolyPhen2",
    "SIFT",
	"FatHmm",
	"FatHmmW",
	"PROVEAN",
	"MutationAssessor",
	"EFIN_SP",
	"EFIN_HD",
	"CADD",
	"PANTHER",
	"PhD-SNP",
	"SNAP",
	"Meta-SNP",
	"CAROL",
	"Condel",
	"COVEC_WMV",
	"MtoolBox",
	"PolyPhen2_transf",
	"SIFT_transf",
	"MutationAssessor_transf",
	"MutationTaster",
	"PCT_AltFreq>=0.1",
	"PCT_AltFreq>=freq_cutoff",
	"AltFreq(0.1~0.9)",
	"AltFreq>=0.9",
);
my @title2=(
	"5' FLANKING SEQUENCE",
	"3' FLANKING SEQUENCE",
	"HOMOLOGY HITS",
);
my $row=0;
for my $i(0..@title-1){
	$SNV->write($row,$i,$title[$i],$format{"title"});
}
for my $i(0..@sample-1){
	$SNV->write_string($row,@title+$i,$sample[$i]."_AltFreq",$format{"title"});
}
for my $i(0..@sample-1){
	$SNV->write_string($row,@title+@sample+$i,$sample[$i]."_Reads",$format{"title"});
}
for my $i(0..@title2-1){
	$SNV->write($row,@title+@sample+@sample+$i,$title2[$i],$format{"title"});
}
$row++;
my $id=1;
foreach my $title(sort keys %library){
	my $weiyi=1;
	foreach my $alt(sort keys %{$library{$title}}){
		my $subID="";$subID=".".$weiyi if(keys %{$library{$title}}>1);
		my $SNVID="SNV"."0"x(5-length($id)).$id;
		$weiyi++;
		$id++;
		my %mut=();
		foreach(split /,/,$analysis{$title}{$alt}{"Mutation 0"}){$mut{(split /:/,$_)[0]}="0:$_" if($_=~ /\w/);}
		foreach(split /,/,$analysis{$title}{$alt}{"Mutation 1"}){$mut{(split /:/,$_)[0]}="1:$_" if($_=~ /\w/);}
		foreach(split /,/,$analysis{$title}{$alt}{"Mutation 2"}){$mut{(split /:/,$_)[0]}="2:$_" if($_=~ /\w/);}
		my $min_varfreq_num = 0; #满足频率过滤条件的样本数
		my $Esample=0;my $PCT10=0;my $PCT_spec=0;my @sampleMid=();my @sampleTop=();####统计突变频率大于特定阈值的样本比例，突变频率大于0.1的样本比例，突变比例在0.1~0.9的样本，突变比例大于0.9的样本
		foreach(@sample){
		    next if(not exists $mut{$_});
			$Esample++;
			my @splits=split /:/,$mut{$_};
			my $refdep=$splits[4];
			my $altdep=$splits[5];
			my $perc=$altdep/($refdep+$altdep);
			$PCT10++ if($perc>=0.1);
			$PCT_spec++ if($perc>=$min_var_freq);
			$min_varfreq_num++ if($perc>=$min_var_freq);
			push @sampleMid,$_ if($perc>0.1 and $perc<0.9);
			push @sampleTop,$_ if($perc>=0.9);		
		}
		next if($min_varfreq_num==0);
		$PCT10=$PCT10/$Esample if($Esample>0);
		$PCT_spec=$PCT_spec/$Esample if($Esample>0);
		my $Gene=$library{$title}{$alt}{"Gene"};
		my $blast=$library{$title}{$alt}{"HOMOLOGY HITS"};
		my $snpcq=$library{$title}{$alt}{"SNP Calling Quality"};
		$snpcq-- if($blast>1);
		my $snvq="H";$snvq="M" if($snpcq==2);$snvq="L" if($snpcq<2);
		my $genoq=$analysis{$title}{$alt}{'Genotyping Quality'};
		my $g1000=0;$g1000=$library{$title}{$alt}{"Freq_Alt (1000g)"} if($library{$title}{$alt}{"Freq_Alt (1000g)"}=~ /\d/);
		my $hgmd=0;$hgmd=1 if(exists $annotation{$title}{$alt}{"HGMD_site_class"} and ($annotation{$title}{$alt}{"HGMD_site_class"}=~ /^\s*DM\s*$/ or ($g1000<0.05 and $annotation{$title}{$alt}{"HGMD_site_class"}=~ /^\s*DM\?\s*$/)));
		my $geneskyHit=0;$geneskyHit=$annotation{$title}{$alt}{"GENESKYDBHITS"} if(exists $annotation{$title}{$alt}{"GENESKYDBHITS"});
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
		$annotation{$title}{$alt}{'phastCons100way'}=(split /\=/,$annotation{$title}{$alt}{'phastCons100way'})[1] if(exists $annotation{$title}{$alt}{'phastCons100way'});
		$annotation{$title}{$alt}{'phyloP100way'}=(split /\=/,$annotation{$title}{$alt}{'phyloP100way'})[1] if(exists $annotation{$title}{$alt}{'phyloP100way'});
		for my $i(0..@title-1){
			my $v=" ";
			$v=$library{$title}{$alt}{$title[$i]} if(exists $library{$title}{$alt}{$title[$i]});
			$v=$annotation{$title}{$alt}{$title[$i]} if(exists $annotation{$title}{$alt}{$title[$i]});
			$v=$analysis{$title}{$alt}{$title[$i]} if(exists $analysis{$title}{$alt}{$title[$i]});
			$v=$SNVID if($title[$i] eq "SNV NO.");
			$v=$PCT10 if($title[$i] eq "PCT_AltFreq>=0.1");
			$v=$PCT_spec if($title[$i] eq "PCT_AltFreq>=freq_cutoff");
			$v=join ",",@sampleMid if($title[$i] eq "AltFreq(0.1~0.9)" and exists($sampleMid[0]));
			$v=join ",",@sampleTop if($title[$i] eq "AltFreq>=0.9" and exists($sampleTop[0]));
			$SNV->write($row,$i,$v,$format{"normal"});
		}
		for my $i(0..@sample-1){
			my ($geno,$ALTP)=(" "," ");
			if(exists $mut{$sample[$i]}){
				my @splits=split /:/,$mut{$sample[$i]};
				my $refdep=$splits[4];
				my $altdep=$splits[5];
				$ALTP=sprintf "%0.5f",$altdep/($refdep+$altdep);
			}
			$SNV->write($row,@title+$i,$ALTP,$format{"left"});
		}
		for my $i(0..@sample-1){
			my ($geno,$reads)=(" "," ");
			if(exists $mut{$sample[$i]}){
				my @splits=split /:/,$mut{$sample[$i]};
				$reads="$splits[2],".$library{$title}{$alt}{'Ref Allele'}.":$splits[4],$alt:$splits[5]";
			}
			$SNV->write($row,@title+@sample+$i,$reads,$format{"left"});
		}
		for my $i(0..@title2-1){
			my $v=" ";
			$v=$library{$title}{$alt}{$title2[$i]} if(exists $library{$title}{$alt}{$title2[$i]});
			$v=$annotation{$title}{$alt}{$title2[$i]} if(exists $annotation{$title}{$alt}{$title2[$i]});
			$v=$analysis{$title}{$alt}{$title2[$i]} if(exists $analysis{$title}{$alt}{$title2[$i]});
			$SNV->write($row,@title+@sample+@sample+$i,$v,$format{"normal"});
		}
		$row++;
	}
}
package::readme::run($Readme,\%format,"readme_snvallMT.txt");