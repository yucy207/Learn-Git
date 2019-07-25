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
	my $groupcount=0;
	while(my $line=<FILE>){
		$line=~ s/[\s\t\r\n]//g;
		next if($line!~/\w/ or $line=~/^#/);
		my @datas=split /=/,$line;
		if($datas[0] eq 'Group' or $datas[0] eq 'Show'){
		   $groupcount++;
		   $hash{$datas[0]}{$groupcount}=$datas[1];
		   $hash{$datas[0]}{$groupcount}=$datas[1];
		   next;
		} 
		if($datas[0] eq 'Sample'){
		   $hash{'Sample'}{$datas[1]}{$datas[2]}=$datas[3];
		   next;
		}
		$hash{$datas[0]}=$datas[1];
	}
	close FILE;
	return %hash;
}

sub readGroup{
	my $file=shift @_;
	my %hash=();
	my $count=0;
	open FILE,$file or die "[ERR] Loss Config File,$file\n";
	while(my $line=<FILE>){
	    $count++;
		$line=~ s/[\s\t\r\n]//g;
		my ($case,$control)=split /;/,$line;
		$hash{$count}{'case'}=$case;
		$hash{$count}{'control'}=$control;
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
	print "RNAOutput:  $RNAOutput\n" if(exists $config->{"RNAOutput"});
	print "Output:     $Output\n";
	print "Report:     $Report\n";
	print "TargetBed:  $config->{'TargetBed'}\n" if(exists($config->{'TargetBed'}));
	print "RealignBed: $config->{'RealignBed'}\n" if(exists($config->{'RealignBed'}));
	print "Species:    $Species\n";
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
	return %library;
}

sub readAnnotation{
	my $file=shift @_;
	my %annotation=();
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$type,$content)=split /\t/,$line,4;
		$annotation{$title}{$alt}{$type}=$content;
	}
	close FILE;
	return %annotation;
}

sub readAnalysis{
	my $file=shift @_;
	my %analysis=();
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


sub checktag{
    my $config=shift @_;
	my @samples;
	my @groups=keys %{$config->{'Group'}};
	my %error;
	foreach my $group(@groups){
	    die "Lost Control\n" if($config->{'Group'}{$group}!~/;/);
	    my @tmps=split /[,;]/,$config->{'Group'}{$group};
		map{push @samples,$_ if($_=~/\w/);}@tmps;
	}
    my $flat=0;# 检错
	foreach my $sample(@samples){
	    my $tag="";$tag=$config->{'Sample'}{$sample}{'tag'} if(exists($config->{'Sample'}{$sample}{'tag'}));
		my $sex="";$sex=$config->{'Sample'}{$sample}{'sex'} if(exists($config->{'Sample'}{$sample}{'sex'}));
		if(not -e $tag or -s $tag==0){
		   $error{$sample}{'tag'}=$tag;$flat=1;
		}
		if($sex ne 'male' and $sex ne 'female'){
		   $error{$sample}{'sex'}=$sex;$flat=1;
		}
	}
	if($flat==1){
	   my @samples=keys %error;
	   print "Sample Error \n";
	   foreach my $sample(@samples){
	       my @infos=keys %{$error{$sample}};
		   map{print "$sample : $_ : $error{$sample}{$_}\n"}@infos;
	   }
	   die;
	}
}

1