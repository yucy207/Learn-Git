package package::main;
use strict;
use warnings;

sub readPara{
	my $file=shift @_;
	my %hash=();
	open FILE,$file or die "[ERR] Loss Para File,$file\n";
	while(my $line=<FILE>){
		$line=~ s/[\s\t\r\n]//g;
		my @split_line=split /\=/,$line;
		next if(@split_line!=3);
		$hash{$split_line[0]}{$split_line[1]}=$split_line[2];
	}
	close FILE;
	return %hash;
}

sub readConfig{
	my $file=shift @_;
	my %hash=();
	open FILE,$file or die "[ERR] Loss Config File,$file\n";
	while(my $line=<FILE>){
		$line=~ s/[\s\t\r\n]//g;
		my ($name,$value)=split /\=/,$line;
		next if(!defined($value));
		if(exists $hash{$name} and ($name eq "case" or $name eq "control") ){
			$hash{$name}.=",".$value;
		}elsif(exists $hash{$name} and $name eq 'SomaticPair'){
			$hash{$name}.=";".$value;
		}else{
			$hash{$name}=$value;
		}
	}
	close FILE;
	return %hash;
}

sub showConfig{
	my $config=shift @_;
	my @case=Sample($config,"case");
	my @control=Sample($config,"control");
	my $nCase=@case;
	my $nControl=@control;
	my $GenderFile="[NA]";$GenderFile=$config->{"GenderFile"} if(exists $config->{"GenderFile"});
	my $Fastq="[NA]";$Fastq=$config->{"Fastq"} if(exists $config->{"Fastq"});
	my $RNAOutput="[NA]";$RNAOutput=$config->{"RNAOutput"} if(exists $config->{"RNAOutput"});
	my $Output="[NA]";$Output=$config->{"Output"} if(exists $config->{"Output"});
	my $Report="[NA]";$Report=$config->{"Report"} if(exists $config->{"Report"});
	my $Species="[NA]";$Species=$config->{"Species"} if(exists $config->{"Species"});
	print "case:       $nCase\n";
	print "control:    $nControl\n";
	print "GenderFile: $GenderFile\n";
	print "Fastq:      $Fastq\n";
	print "RNAOutput:  $RNAOutput\n";
	print "Output:     $Output\n";
	print "Report:     $Report\n";
	print "Species:    $Species\n";
}

# 获取当前项目分析的样本
sub get_sample{
    my $hashConfig = shift @_;
    my @need_types = @_;
    my @samples;

    foreach my $need_type(@need_types)
    {
        next if(!exists($hashConfig->{$need_type}));
        foreach my $sample(split /,/, $hashConfig->{$need_type})
        {
            push @samples, $sample if($sample=~/\w/);
        }
    }

    return @samples;
}

# 创建目录
sub make_dir{
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
}

# 检验目录是否存在
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# 检验文件是否为空
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}
sub Sample{
	my $config=shift @_;
	my $key=shift @_;
	my @array=();
	if(exists $config->{$key}){
		foreach(split /,/,$config->{$key}){
			push @array,$_ if($_=~ /\w/);
		}
	}
	return @array;
}

sub SamplePair{
	my $config = shift @_;
	my %hashPairSample;
    my $pairInfo=(exists($config->{'SomaticPair'})) ? $config->{'SomaticPair'} : "";
	foreach my $pairSample((split /;/,$pairInfo)){
		next if($pairSample!~/\w/ or $pairSample!~/,/);
		my ($tumorSample,$normalSample)=split /,/,$pairSample;
		my $pairSampleName=$tumorSample.".VS.".$normalSample;
		$hashPairSample{$pairSampleName}{"Case"}    =$tumorSample;
		$hashPairSample{$pairSampleName}{"Control"} =$normalSample;
	}
	return %hashPairSample;
}

sub checkFileSizeOK{
	my $file=shift @_;
	if(-e $file and -s $file>0){
		return 1;
	}else{
		return 0;
	}
}

sub checkBedFileOK{
	my $file=shift @_;
	if(-e $file){
		my %tmp=();
		my $chrNow="";
		open FILE,$file;
		while(my $line=<FILE>){
			next if($line=~ /^\@/);
			$line=~ s/[\r\n]//g;
			$line=~ s/\s+$//;
			my @split_line=split /\t/,$line;
			# 数据信息不正确，返回错误
			return 0 if(@split_line<5 and $line=~ /\w/);
			# 不允许染色体不连续
			$chrNow=$split_line[0] if(not exists $tmp{"chr"}{$split_line[0]});
			return 0 if($chrNow ne $split_line[0]);
			# 不允许染色体位置不排序
			return 0 if(exists $tmp{"chr"}{$split_line[0]} and $split_line[1]<$tmp{"chr"}{$split_line[0]});
			$tmp{"chr"}{$split_line[0]}=$split_line[1];
		}
		close FILE;
	}else{
		return 0;
	}
}

sub readLibrary{
	my $file=shift @_;
	my %library=();
	print "Reading $file...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$q,$rs,$freq,$ref,$alt,$chr,$pos,$gene,$region,$function,$hgvs,$seq5,$seq3,$blast)=split /\t/,$line;
		$library{$title}{$alt}{"SNP Calling Quality"}=$q;
		$library{$title}{$alt}{"SNP ID"}=$rs;
		$library{$title}{$alt}{"Freq_Alt (1000g)"}=$freq;
		$library{$title}{$alt}{"Ref Allele"}=$ref;
		$library{$title}{$alt}{"Alt Allele"}=$alt;
		$library{$title}{$alt}{"Chrs"}=$chr;
		$library{$title}{$alt}{"Position"}=$pos;
		$library{$title}{$alt}{"Position2"}=(split /\-/,$pos)[0];
		$library{$title}{$alt}{"Gene"}=$gene;
		$library{$title}{$alt}{"Gene Region"}=$region;
		$library{$title}{$alt}{"Function"}=$function;
		$library{$title}{$alt}{"Predicted Protein Variants"}=$hgvs;
		$library{$title}{$alt}{"5' FLANKING SEQUENCE"}=$seq5;
		$library{$title}{$alt}{"3' FLANKING SEQUENCE"}=$seq3;
		$library{$title}{$alt}{"HOMOLOGY HITS"}=$blast;
	}
	close FILE;
	print "OK\n";
	return %library;
}
sub readLibrarySomatic{
	my $file=shift @_;
	my %library=();
	print "Reading $file ...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$quality,@somaticInfos)=split /\t/,$line;
		foreach my $somaticInfo(@somaticInfos){
			my ($pairSample,$tmp)=split /;/,$somaticInfo,2;
		    $library{$title}{$alt}{$pairSample}=$tmp;
		}
	}
	close FILE;
	print " OK\n";
	return %library;
}

sub readAnnotation{
	my $file=shift @_;
	my %annotation=();
	print "Reading $file...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$type,$content)=split /\t/,$line,4;
		if($type eq "HGMD_site" or $type eq "HGMD_site_class"){
		   next if(exists($annotation{$title}{$alt}{'HGMD_site_class'}) and $annotation{$title}{$alt}{'HGMD_site_class'}=~/DM/);
		}elsif($type eq 'Strand Orientation'){# 表头名称发生变动
		   $annotation{$title}{$alt}{"Gene Strand Orientation"}=$content;
		   next;
		}
		$annotation{$title}{$alt}{$type}=$content;
	}
	close FILE;
	print "OK\n";
	return %annotation;
}
sub readSomaticAnnotation{
	my $file=shift @_;
	my $hashSomatic = shift @_;
	print "Reading Somatic Mutations From $file ...";
	my %annotation=();
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$type,$content)=split /\t/,$line,4;
		next if(!exists($hashSomatic->{$title}{$alt}));
		$annotation{$title}{$alt}{$type}=$content;
	}
	close FILE;
	print " OK\n";
	return %annotation;
}
sub readAnalysis{
	my $file=shift @_;
	my %analysis=();
	print "Reading $file...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$genoQ,$s0,$s1,$s2)=split /\t/,$line;
		$analysis{$title}{$alt}{"Genotyping Quality"}=$genoQ;
		$analysis{$title}{$alt}{"Mutation 0"}=$s0;
		$analysis{$title}{$alt}{"Mutation 1"}=$s1;
		$analysis{$title}{$alt}{"Mutation 2"}=$s2;
	}
	close FILE;
	print "OK\n";
	return %analysis;
}
sub readSomaticAnalysis{
	my $file=shift @_;
	my $hashSomatic = shift @_;
	my %analysis=();
	print "Reading $file...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$genoQ,$s0,$s1,$s2)=split /\t/,$line;
		next if(!exists($hashSomatic->{$title}{$alt}));
		$analysis{$title}{$alt}{"Genotyping Quality"}=$genoQ;
		$analysis{$title}{$alt}{"Mutation 0"}=$s0;
		$analysis{$title}{$alt}{"Mutation 1"}=$s1;
		$analysis{$title}{$alt}{"Mutation 2"}=$s2;
	}
	close FILE;
	print "OK\n";
	return %analysis;
}
sub readRunConfig{
	my $file=shift @_;
	my %run=();
	my $runID="";
	my %weiyi=();
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		$line=~ s/^\s+//;
		$line=~ s/\s+$//;
		next if($line=~ /^\#/);
		if($line=~ /::RUN/){
			$runID=$line;
			$runID=~ s/::RUN//;
			$runID=~ s/[\s\t]//g;
			%weiyi=();
			next;
		}
		if($runID ne ""){
			my @split_line=split /:/,$line;
			next if(@split_line!=2);
			$weiyi{$split_line[0]}++;
			$run{$runID}{$split_line[0]}{$weiyi{$split_line[0]}}=$split_line[1];
		}
	}
	close FILE;
	return %run;
}

sub realMut{
	my ($Zygotes,$title,$Gender)=@_;
	if($Zygotes eq "HET"){
		return 1;
	}elsif($Zygotes eq "HOMA"){
		if($title=~ /[XY]/i and $Gender eq "male"){
			return 1;
		}else{
			return 2;
		}
	}
	return 0;
}


1