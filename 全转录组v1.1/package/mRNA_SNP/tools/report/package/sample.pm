package package::sample;
use strict;
use warnings;

sub run{
	my $SAMPLE=shift @_;
	my $format=shift @_;
	my $library=shift @_;
	my $annotation=shift @_;
	my $analysis=shift @_;
	my $filter=shift @_;
	my $gender=shift @_;
	my $run=shift @_;
	my $hashassigngender=shift @_;
	$SAMPLE->set_row(0, 60);
	my @case=();@case=split /,/,$run->{"Case"}{1} if(exists $run->{"Case"}{1});
	my @control=();@control=split /,/,$run->{"Control"}{1} if(exists $run->{"Control"}{1});
	my @sample=(@case,@control);
	my %type=();
	foreach(@case){$type{$_}="case";}
	foreach(@control){$type{$_}="control";}
	my %tmp=();
	foreach my $id(sort{$a <=> $b}keys %{$filter}){
		foreach my $title(keys %{$filter->{$id}}){
			foreach my $alt(keys %{$filter->{$id}{$title}}){
			    next if(!exists($library->{$title}{$alt}->{"Predicted Protein Variants"}) or $library->{$title}{$alt}->{"Predicted Protein Variants"}!~/\w/);
				my $Priority=$filter->{$id}{$title}{$alt}{"First Priority"};
				my $Gene=$library->{$title}{$alt}{"Gene"};
				my %mut=();
				foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 0"}){$mut{(split /:/,$_)[0]}=0 if($_=~ /\w/);}
				foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 1"}){$mut{(split /:/,$_)[0]}=1 if($_=~ /\w/);}
				foreach(split /,/,$analysis->{$title}{$alt}{"Mutation 2"}){$mut{(split /:/,$_)[0]}=2 if($_=~ /\w/);}
				for my $i(0..@sample-1){
					next if(not exists $mut{$sample[$i]});
					$tmp{"priority"}{$sample[$i]}{$Gene}.=$Priority;
					$tmp{"hgmd"}{$sample[$i]}{$Gene}{$title}="";
					$tmp{"hgmd"}{$sample[$i]}{$Gene}{$title}=$annotation->{$title}{$alt}{"HGMD_site"} if(exists $annotation->{$title}{$alt}{"HGMD_site"});
					my ($ppv)=split /[,;]/,$library->{$title}{$alt}->{"Predicted Protein Variants"};$ppv=$Gene.":".$ppv if($ppv!~ /^$Gene:/);
					$tmp{"ppv"}{$sample[$i]}{$Gene}{$title}=$ppv;
					if($mut{$sample[$i]}==1){
						$tmp{"count"}{$sample[$i]}{$Gene}+=1;
						$tmp{"title"}{$sample[$i]}{$Gene}{$title}=1;
						my $ppv=$tmp{"ppv"}{$sample[$i]}{$Gene}{$title};
						$tmp{"ppvCase1"}{$Gene}{$sample[$i]}{$ppv}=1 if($type{$sample[$i]} eq "case");
						$tmp{"ppvControl1"}{$Gene}{$sample[$i]}{$ppv}=1 if($type{$sample[$i]} eq "control");
					}
					if($mut{$sample[$i]}==2){
						$tmp{"count"}{$sample[$i]}{$Gene}+=2;
						$tmp{"title"}{$sample[$i]}{$Gene}{$title}=2;
						my $ppv=$tmp{"ppv"}{$sample[$i]}{$Gene}{$title};
						$tmp{"ppvCase2"}{$Gene}{$sample[$i]}{$ppv}=1 if($type{$sample[$i]} eq "case");
						$tmp{"ppvControl2"}{$Gene}{$sample[$i]}{$ppv}=1 if($type{$sample[$i]} eq "control");
					}
				}
			}
		}
	}
	# ppv count
	my %ppvGroup=();
	foreach my $gene(keys %{$tmp{"ppvCase1"}}){
		foreach my $case(keys %{$tmp{"ppvCase1"}{$gene}}){
			foreach my $key1(keys %{$tmp{"ppvCase1"}{$gene}{$case}}){
				foreach my $key2(keys %{$tmp{"ppvCase1"}{$gene}{$case}}){
					next if($key1 eq $key2);
					foreach my $control(keys %{$tmp{"ppvControl1"}{$gene}}){
						if(exists $tmp{"ppvControl1"}{$gene}{$control}{$key1} and exists $tmp{"ppvControl1"}{$gene}{$control}{$key2}){
							$ppvGroup{$case}{$gene}{$key1}=1;
							$ppvGroup{$case}{$gene}{$key2}=1;
							$ppvGroup{$control}{$gene}{$key1}=1;
							$ppvGroup{$control}{$gene}{$key2}=1;
						}
					}
				}
			}
		}
	}
	foreach my $gene(keys %{$tmp{"ppvCase2"}}){
		foreach my $case(keys %{$tmp{"ppvCase2"}{$gene}}){
			foreach my $key1(keys %{$tmp{"ppvCase2"}{$gene}{$case}}){
				foreach my $control(keys %{$tmp{"ppvControl2"}{$gene}}){
					if(exists $tmp{"ppvControl2"}{$gene}{$control}{$key1}){
						$ppvGroup{$case}{$gene}{$key1}=1;
						$ppvGroup{$control}{$gene}{$key1}=1;
					}
				}
			}
		}
	}
	# output
	my %hash_count=();
	foreach my $key1(keys %{$tmp{"count"}}){
		foreach my $key2(keys %{$tmp{"count"}{$key1}}){
			next if($tmp{"count"}{$key1}{$key2}==1 and exists $run->{"Model"}{1} and $run->{"Model"}{1} eq "Recessive");
			$hash_count{"gene"}{"all"}{$key2}+=1;
			$hash_count{"sample"}{$key1}+=1;
			if($type{$key1} eq "case"){$hash_count{"gene"}{"case"}{$key2}+=1;}
			elsif($type{$key1} eq "control"){$hash_count{"gene"}{"control"}{$key2}+=1;}
			$hash_count{"exists"}{$key1}{$key2}=1;
		}
	}
	$SAMPLE->write(0, 0,"Sample Name",$format->{'title'});
	$SAMPLE->write(0, 1,"Gender",$format->{'title'});
	$SAMPLE->write(0, 2,"Case OR Control",$format->{'title'});
	$SAMPLE->write(0, 3,"Mutation Gene",$format->{'title'});
	$SAMPLE->write(0, 4,"Number\n(Mutation Gene)",$format->{'title'});
	$SAMPLE->write(1, 4,"Summary(Case)",$format->{'normal'});
	$SAMPLE->write(2, 4,"Summary(Control)",$format->{'normal'});
	for my $i(0..3){
		$SAMPLE->write(1, $i," ",$format->{'normal'});
		$SAMPLE->write(2, $i," ",$format->{'normal'});
	}
	my $col=5;
	foreach my $key1(sort{$a cmp $b}keys %{$hash_count{"gene"}{"all"}}){
		$SAMPLE->write_string(0, $col,$key1,$format->{'title'});
		$SAMPLE->write(1, $col,$hash_count{"gene"}{"case"}{$key1},$format->{'normal'});
		$SAMPLE->write(2, $col,$hash_count{"gene"}{"control"}{$key1},$format->{'normal'});
		$col++;
	}
	my $row=3;
	foreach my $key1(@sample){
		next if(not exists $hash_count{"sample"}{$key1});
		$col=0;
		$SAMPLE->write_string($row, $col,$key1,$format->{'normal'});$col++;
		my $getgender=(exists($hashassigngender->{$key1})) ? $hashassigngender->{$key1}:"Unknown";
		$SAMPLE->write($row, $col,$getgender,$format->{'normal'});$col++;
		$SAMPLE->write($row, $col,$type{$key1},$format->{'normal'});$col++;
		my $mutgene="";
		foreach my $key2(sort{$a cmp $b}keys %{$hash_count{"gene"}{"all"}}){
			if(exists $hash_count{"exists"}{$key1}{$key2}){
				$mutgene.=$key2.",";
			}
		}
		$SAMPLE->write($row, $col,$mutgene,$format->{'normal'});$col++;
		$SAMPLE->write($row, $col,$hash_count{"sample"}{$key1},$format->{'normal'});$col++;
		foreach my $key2(sort{$a cmp $b}keys %{$hash_count{"gene"}{"all"}}){
			if(exists $hash_count{"exists"}{$key1}{$key2}){
				my @array=();
				foreach my $title(keys %{$tmp{"title"}{$key1}{$key2}}){
					my $ppv=$tmp{"ppv"}{$key1}{$key2}{$title};
					my $style="";
					$style.="bold" if($tmp{"hgmd"}{$key1}{$key2}{$title} ne "");
					$style.="blue" if($tmp{"title"}{$key1}{$key2}{$title}>=2);
					$style.="italic" if(exists $ppvGroup{$key1}{$key2}{$ppv});
					$style="left" if($style eq "");
					push @array,$format->{$style};
					push @array,"$ppv,";
				}
				my $firstelse=$tmp{"priority"}{$key1}{$key2};
				$firstelse=~ s/First1//g;
				my $numfirst=()=$tmp{"priority"}{$key1}{$key2}=~ /First1/g;
				my $redfirst=1;$redfirst=2 if(exists $run->{"Model"}{1} and $run->{"Model"}{1} eq "Recessive");
				push @array,$format->{'orange'} if($numfirst>=$redfirst or $firstelse!~ /\w/);
				$SAMPLE->write_rich_string($row, $col,@array);$col++;
			}
			else{
				$SAMPLE->write($row, $col,"",$format->{'normal'});$col++;
			}
		}
		$row++;
	}

	
	
}


1