package package::nipt::main;
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
	my %exists=();
	my $weiyi=0;
	open FILE,$file or die "[ERR] Loss Config File,$file\n";
	while(my $line=<FILE>){
		next if($line=~ /\#/);
		if($line=~ /\=/){
			$line=~ s/[\s\t\r\n]//g;
			my @split_line=split /\=/,$line;
			next if(@split_line!=2);
			$hash{$split_line[0]}=$split_line[1];
		}else{
			$line=~ s/[\r\n]//g;
			$line=~ s/^[\s\t]+//;
			$line=~ s/[\s\t]+$//;
			my @split_line=split /[\s\t]+/,$line;
			next if(@split_line!=2);
			print "[WARN] $split_line[0] is in more than 1 Group\n" if(exists $exists{$split_line[0]});
			$exists{$split_line[0]}=1;
			$hash{"Sample"}{$split_line[1]}{$weiyi}=$split_line[0];
			$weiyi++;
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
	my $Output="[NA]";$Output=$config->{"Output"} if(exists $config->{"Output"});
	my $Report="[NA]";$Report=$config->{"Report"} if(exists $config->{"Report"});
	my $Species="[NA]";$Species=$config->{"Species"} if(exists $config->{"Species"});
	print "case:       $nCase\n";
	print "control:    $nControl\n";
	print "GenderFile: $GenderFile\n";
	print "Fastq:      $Fastq\n";
	print "Output:     $Output\n";
	print "Report:     $Report\n";
	print "Species:    $Species\n";
	die "[ERR] must have case and control\n" if($nCase==0 or $nControl==0);
}

sub Sample{
	my $config=shift @_;
	my $key=shift @_;
	my %hash=();
	my @array=();
	my %exists=();
	if(exists $config->{"Sample"}){
		foreach my $group(keys %{$config->{"Sample"}}){
			foreach my $weiyi(sort{$a<=>$b}keys %{$config->{"Sample"}{$group}}){
				my $sample=$config->{"Sample"}{$group}{$weiyi};
				next if(exists $exists{$sample});
				$exists{$sample}=1;
				if($key eq "case"){
					push @array,$sample if($group ne "Control");
				}elsif($key eq "control"){
					push @array,$sample if($group eq "Control");
				}
			}
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



1