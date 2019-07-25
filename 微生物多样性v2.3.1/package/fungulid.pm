package fungulid;

use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata    = shift;
	my $base        = shift;

	my @groups    = @{$metadata->{'group'}};
	my $group_num = @groups;
	my $manager   = new Parallel::ForkManager($group_num);
	
	for my $group (@groups){
		
		my $groupfile = (split /:/,$group)[0];
		my $groupbase = (split /:/,$group)[1];
		
		die "ERROR: Can not find $groupfile, Please Check It!\n" if not -e $groupfile;

		my @groupname     = split /\//,$groupfile;
		my ($groupname)   = $groupname[-1] =~ /(.*).xls/;
		my $intm_result   = qq{$metadata->{intm_result}/$groupname};
		my $projectpath   = qq{$metadata->{result}/$groupname};

		if (-e qq{$intm_result/tmp/FG/fungulid.finish}){

			print qq{$groupname fungulid分析已经运行完成!\n 是否重新运行? [y/n]};
			my $option = <STDIN>;
			$option =~ s/[\r\n]//g;
			next if $option ne "y";
			system qq{rm $intm_result/tmp/FG/fungulid.finish};
			system qq{rm $intm_result/tmp/FG/fungulid.log};
		}

		# Forks and returns the pid for the child:
		$manager->start and next;

		### fungulid分析
		fungulid($intm_result,$projectpath,$groupfile,$groupbase,$metadata,$base);
		
		system qq{touch $intm_result/tmp/FG/fungulid.finish};
		print qq{$groupname fungulid分析已经运行完成!\n};

		$manager->finish;
	}
	$manager->wait_all_children;
	print qq{fungulid分析已经运行完成!\n};
}

sub fungulid
{
	my $intm_result = shift;
	my $projectpath = shift;
	my $groupfile   = shift;
	my $groupbase   = shift;
	my $metadata    = shift;
	my $base        = shift;

	my $dissi       = qq{$metadata->{dissi}};
	my $qiime		= qq{$base->{qiime_bin}};
	my $Rscript     = qq{$base->{Rscript_bin}};
	my $Rscript_3_5_3     = qq{$base->{Rscript_3_5_3}};
	my $scriptdir   = qq{$base->{scriptdir}};

	my $FUNGuild = qq{$projectpath/FUNGuild};
	system qq{mkdir -p $FUNGuild} if not -d $FUNGuild;

	my $tmpPath = qq{$intm_result/tmp/FG};
	system qq{mkdir -p $tmpPath} if not -d $tmpPath;

	#my $otufile = qq{$intm_result/tmp/subsample_otu.tax.$dissi.xls};
	my $otufile = -e qq{$intm_result/tmp/subsample_otu.tax.$dissi.xls} ? qq{$intm_result/tmp/subsample_otu.tax.$dissi.xls} : qq{$intm_result/tmp/otu.tax.$dissi.xls};

	die "ERROR: Can not find $otufile , Please Check It!\n" if not -e $otufile;

	#system qq{touch $tmpPath/fungulid.finish};
	return if -e qq{$tmpPath/fungulid.finish};

	### Run 
	### Format ###
	otu_format($otufile, "$tmpPath/otu.format.xls");
	my $id = `id`;
	my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;
	system qq{$qiime bash -c 'python $scriptdir/Guilds_v1.1.py -otu $tmpPath/otu.format.xls -db fungi;chown -R $uid:$gid $tmpPath/otu.format.guilds.txt' &>> $tmpPath/fungulid.log};
	### Format ###
	format_guilds("$tmpPath/otu.format.guilds.txt", "$FUNGuild/guilds.xls");
	#system qq{rm $tmpPath/otu.format.guilds.txt};

	### plot barplot ###
	system qq{$Rscript $scriptdir/FUNGuild.barplot.r $FUNGuild/guilds.xls $groupfile $FUNGuild &>> $tmpPath/fungulid.log};

	###
	system qq{mkdir -p $FUNGuild/Heatmap};
	system qq{$Rscript_3_5_3 $scriptdir/16S.PI.heatmap.r $FUNGuild/guilds.xls  $groupfile 30 "$FUNGuild/Heatmap" &>> $tmpPath/fungulid.log};

	
	###差异分析
	my @samplegroup = samplegroup($groupfile,$groupbase,$intm_result);
	
	for my $samplegroup (@samplegroup){
		
		my $path = "";
		
		($path) = (split /\//,$samplegroup)[-1] =~ /(.*).xls/ if ((split /\//,$samplegroup)[-1] =~ /_vs_/);
	
		### t_test ###
		system qq{$Rscript $scriptdir/FUNGuild_pair_diff_analysis.r $FUNGuild/guilds.xls $samplegroup $FUNGuild &>> $tmpPath/fungulid.log};
		### anova ###
		my $num = `less $samplegroup | cut -f 2 | sort |  uniq | wc -l`;
		if ($num >2 ){
			system qq{$Rscript $scriptdir/FUNGuild.ANOVA.DEA.r $FUNGuild/guilds.xls $samplegroup P 0.05 $FUNGuild/$path &>> $tmpPath/fungulid.log};
		}
	}
}

sub otu_format
{
	my $in  = shift;
	my $out = shift;

	open(IN, $in) || die "cannot open inputfile\n";
	open(OUT,">$out") || die "cannot make format_file\n";
	my @head = split/\t/,<IN>;
	$head[0] =~s/OTUId/OTU_ID/;
	my ($tax_index) = grep{$head[$_] eq 'Taxonomy'} 0..$#head;
	my $head = join("\t",@head[0..$tax_index-2],$head[$tax_index]);
	print OUT "$head\n";
	while(<IN>){

		$_ =~s/[\r\n]//g;
		my @a = split/\t/,$_;
		my @line = grep(/.+/,@a);#去除空的元素
		my $tax = join(";",@line[$tax_index+1..$#line]);
		my $res = join("\t",@line[0..$tax_index-2],$tax);
		print OUT "$res\n";
	}
	close IN;
	close OUT;
}

sub format_guilds
{
	my $in  = shift;
	my $out = shift;

	open(IN, $in) || die "cannot open inputfile\n";
	open(OUT,">$out") || die "cannot make outfile\n";
	my @head1 	= split/\t/,<IN>;
	my ($tax) 	= grep{$head1[$_] eq 'Taxonomy'} 0..$#head1;
	my ($guild) = grep{$head1[$_] eq 'Guild'} 0..$#head1;
	my ($confidence) = grep{$head1[$_] eq 'Confidence Ranking'} 0..$#head1;
	my $head1 = join("\t","Guild",@head1[1..$tax-1]);
	print OUT "$head1\n";

	my %hash = ();
	while(<IN>){

		my @line = split/\t/,$_;
		next if $line[$confidence] eq "-" || $line[$confidence] eq "Possible";
		$line[$guild] =~s/\-/ \/ /g;
		foreach my $i (1..$tax-1){
			$hash{$line[$guild]}{$i} += $line[$i];
		}

	}

	foreach my $x(sort keys %hash){
		my $v = "$x";
		foreach my $y(sort {$a<=>$b} keys %{$hash{$x}}){
			$v.="\t$hash{$x}{$y}";
		}
		print OUT"$v\n";
	}
	close IN;
	close OUT;
}

sub samplegroup	{
	
	my $groupfile = shift;
	my $groupbase = shift;
	my $grouppath = shift;
	
	my @samplegroup;
	
	push @samplegroup, $groupfile;
	
	if (defined $groupbase){
		
		my @groupname = split /;/,$groupbase;
		
		for my $groupname (@groupname){

			my @group = split /,/, $groupname;
			
			my $group = join "_vs_", @group; 
			
			my $samplegroup = qq{$grouppath/$group\.xls};
			
			open FILE ,$groupfile;
			open OUT,">$samplegroup";
			while (<FILE>) {
				
				s/[\r\n]//g;
				
				my @arr = split /\t/, $_;
				
				print OUT qq{$_\n} if ($arr[1] ~~ @group);

			}
			close FILE;
			close OUT;
			push @samplegroup, $samplegroup;
			
		}
		
	}
	# else{
		
		# push @samplegroup, $groupfile;
	# }
	
	return @samplegroup;
}

1;
