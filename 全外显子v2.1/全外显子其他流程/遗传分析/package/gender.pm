package package::gender;
use strict;
use warnings;

sub read{
	my $file=shift @_;
	my %gender=();
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		$line=~ s/^\s+//;
		$line=~ s/\s+$//g;
		my @split_line=split /\s+/,$line,2;
		next if(@split_line!=2);
		$gender{$split_line[0]}="male" if($split_line[1] eq "1" or $split_line[1] eq "male");
		$gender{$split_line[0]}="female" if($split_line[1] eq "0" or $split_line[1] eq "female");
	}
	close FILE;
	return %gender;
}

sub read2{
	my $config=shift @_;
	my %gender=();
	my $genderfile=(exists($config->{'GenderFile'})) ? $config->{'GenderFile'}:"/";
	my $genderfilelog=$config->{'Report'}."/gender.log";
	if(-f $genderfile){
	   open IN,$genderfile;
	   while(<IN>){
	      $_=~s/[\r\n]//g;
		  next if($_!~/\w/ or $_=~/^#/);
		  my @split_line=split /\s+/,$_,2;
		  next if(@split_line!=2);
		  $split_line[0]=~s/[\r\n\s]//g;
		  $split_line[1]=~s/[\r\n\s]//g;
		  $gender{$split_line[0]}="male" if($split_line[1] eq "1" or $split_line[1] eq "male");
		  $gender{$split_line[0]}="female" if($split_line[1] eq "0" or $split_line[1] eq "female");
	   }
	   close IN;
	}
	if(-f $genderfilelog){
	   open IN,$genderfilelog;
	   while(<IN>){
	      $_=~s/[\r\n]//g;
		  next if($_!~/\w/ or $_=~/^#/);
		  my @split_line=split /\s+/,$_,2;
		  next if(@split_line!=2);
		  $split_line[0]=~s/[\r\n\s]//g;
		  $split_line[1]=~s/[\r\n\s]//g;
          if($split_line[1]=~/\// and $split_line[1]=~/=/){
		     my ($y,$x,$value)=split /[\/=]/,$split_line[1];
			 $gender{$split_line[0]}="male" if($value>0.4);
			 $gender{$split_line[0]}="female" if($value<0.1);
		  }
	   }
	   close IN;
	}	
	return %gender;
}

sub get{
	my $config=shift @_;
	my @case=package::main::Sample($config,"case");
	my @control=package::main::Sample($config,"control");
	my @sample=(@case,@control);
	open LOG,">".$config->{"Report"}."/gender.log";
	if(exists $config->{"GenderFile"} and -e $config->{"GenderFile"}){
		my %gender=package::gender::read($config->{"GenderFile"});
		my $isok=1;
		foreach(@sample){
			$isok=0 if(not exists $gender{$_});
		}
		if($isok==1){
			print LOG "# User Defined ".$config->{"GenderFile"}."\n";
			foreach(@sample){
				print LOG "$_\t$gender{$_}\n";
			}
		}else{
			subGet($config);
		}
	}else{
		subGet($config);
	}
	close LOG;
}

sub subGet{
	print LOG "# Auto Analysis\n";
	my $config=shift @_;
	my @case=package::main::Sample($config,"case");
	my @control=package::main::Sample($config,"control");
	my @sample=(@case,@control);
	print "\nProcess...\n";
	my $out="";
	my ($good,$bad)=(0,0);
	my %gender=();
	foreach(@sample){
		if(-e $config->{"Report"}."/status/$_.tag"){
			my $v=subRun($_,$config->{"Report"}."/status/$_.tag");
			if($v eq "male" or $v eq "female"){
				$gender{$_}=$v;
				$good++;
			}else{
				$out.="$_,";
				$bad++;
			}
		}else{
			$out.="$_,";
			$bad++;
		}
	}
	if($out ne ""){
		print "[WARN] Good: $good, Bad:$bad, Sample:$out\n";
		print "\n[Option] Set All Bad to Female?[y/n]";
		my $input=package::input::run();
		if($input eq "y"){
			print LOG "# User Defined All Female\n";
			foreach(@sample){
				next if(exists $gender{$_});
				$gender{$_}="female";
			}
		}
	}
	print LOG "\n";
	foreach(@sample){
		my $v="Unknown";$v=$gender{$_} if(exists $gender{$_});
		print LOG "$_\t$v\n";
	}
}

sub subRun{
	my $name=shift @_;
	my $file=shift @_;
	my (%hash,%hash_num)=();
	open FILE,$file;
	while(my $line=<FILE>){
		my ($chr,$tp)=split /\t/,$line,2;
		next if($chr ne "X" and $chr ne "Y");
		$line=~ s/[\r\n]//g;
		next if($line!~ /\d/);
		my @split_line=split /\t/,$line;
		$hash{$split_line[0]}+=$split_line[7];
		$hash_num{$split_line[0]}++;
	}
	close FILE;
	return -1 if(not exists $hash{"X"} or not exists $hash{"Y"} or $hash{"X"}==0 or $hash{"Y"}==0);
	my $value=($hash{"Y"}/$hash_num{"Y"})/($hash{"X"}/$hash_num{"X"});
	my $type="";
	print LOG "$name\t".($hash{"Y"}/$hash_num{"Y"})."/".($hash{"X"}/$hash_num{"X"})."=".$value."\n";
	if($value<0.1){$type="female";}
	elsif($value>0.4){$type="male";}
	else{
		print "[WARN] $name Y/X=$value, cannot be confirm!please choose\n";
		print "1:male\n";
		print "2:female\n";
		print "3:Unknown\n";
		print "you choose[1/2/3]:";
		my $input=package::input::run();
		if($input eq 1){$type="male";print LOG "[WARN] confirm $name as male\n";}
		elsif($input eq 2){$type="female";print LOG "[WARN] confirm $name as female\n";}
	}
	if($type=~ /\w/){return $type;}
	else{return -1;}
}


1