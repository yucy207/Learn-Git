
use strict;
use Cwd 'abs_path';

die "Usage: 
	perl 16S.BetaDiversity.Lefse.pl otustatRAfile groupfile outputpath 
		parameters ->
		    otustatRAfile: [file -> always otu.tax.0.03.stat.RA.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n" if (@ARGV!=3);
	
my $otustatRAfile = $ARGV[0];	
my $groupfile = $ARGV[1];
my $outputpath = $ARGV[2];
$outputpath = abs_path($outputpath);

my $id = `id`;
my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;

# my $LefseDIR = "/home/panrf/Softwares/nsegata-lefse/";	
# die "ERROR:*** Can not find Lefse sofeware in $LefseDIR, please check it!***\n" if not -e $LefseDIR;
	
#my ($dissi) = $otustatRAfile =~ /\.([\d\.]+)\./;

open GROUP, $groupfile or die "Can not open $groupfile, Please Check It!\n";
open STAT, $otustatRAfile or die "Can not open $otustatRAfile, Please Check It!\n";
open FILE, qq{>$outputpath/taxon.stat.for.lefse.xls} or die "Can not write to file!\n";

############ prepare data for lefse
my %sampleinfo;
while(<GROUP>){

        s/\s+$//g;
	chomp;
	next if /^\s+$/;
	
	my @data = split /\s+/;
			
	die "Sample name duplication, please check sample.groups!\n" if exists $sampleinfo{$data[0]};
	
	$sampleinfo{$data[0]} = $data[1];

}
close GROUP;

while(<STAT>){

	chomp;
	s/\s+$//;
	my @names = split /\t/;
	
	if (/Taxlevel/){
	
		for my $i (1..$#names){
		
			if (exists $sampleinfo{$names[$i]}){
			
				$names[$i] = $sampleinfo{$names[$i]};       # change sample name to groups name 
			
			}else{
			
				die "Sample names in $groupfile and $otustatRAfile do not match, Check It !!!\n";
			
			}
		
		}
		
	}

	my $out = join "\t",@names[0..$#names];    # do not export Size column
 	print FILE qq{$out\n};

}

close STAT;
close FILE;

############ lefse
chdir $outputpath;
system qq{docker run --rm -u $uid:$gid -v /home:/home lefse format_input.py $outputpath/taxon.stat.for.lefse.xls $outputpath/taxon.stat.format.xls -c 1 -o 1000000};
system qq{docker run --rm -u $uid:$gid -v /home:/home lefse run_lefse.py $outputpath/taxon.stat.format.xls $outputpath/taxon.stat.res.xls};

open(STAT,"$outputpath/taxon.stat.res.xls") || die "cannot open taxon.stat.res.xls\n";
open(OUT,">$outputpath/taxon.stat.res1.xls") || die "cannot make taxon.stat.res1.xls\n";
while(my $line = <STAT>){
	my $taxon = (split/\t/,$line)[0];
	my @arr = split/\./,$taxon;
	next if $arr[$#arr] =~/Unassigned/;
	print OUT "$line";
}
close STAT;
close OUT;

open(STAT,"$outputpath/taxon.stat.res.xls") || die "cannot open taxon.stat.res.xls\n";
	open(OUT,">$outputpath/taxon.stat.res2.xls") || die "cannot make taxon.stat.res2.xls\n";
	while(my $line = <STAT>){
		my $taxon = (split/\t/,$line)[0];
		my @arr = split/\./,$taxon;
		next if $arr[$#arr] =~/Unassigned/;
		next if $arr[$#arr] =~/s__/;
		print OUT "$line";
	}
close STAT;
close OUT;

system qq{docker run --rm -u $uid:$gid -v /home:/home zhengy/r-3.4.3:v1 Rscript /home/panrf/LearningPipeline/Metagenomics/lefse4plot.r $outputpath/taxon.stat.res1.xls 50 $outputpath};
system qq{docker run --rm -u $uid:$gid -v /home:/home lefse plot_cladogram.py $outputpath/taxon.stat.res4plot.xls $outputpath/taxon.stat.cladogram.pdf --format pdf --dpi 600 --title_font_size 8 --label_font_size 3 --class_legend_font_size 5 --right_space_prop 0.25 --left_space_prop 0.1};

system qq{docker run --rm -u $uid:$gid -v /home:/home zhengy/r-3.4.3:v1 Rscript /home/panrf/LearningPipeline/Metagenomics/lefse4plot.r $outputpath/taxon.stat.res2.xls 10 $outputpath};
system qq{docker run --rm -u $uid:$gid -v /home:/home lefse plot_res.py $outputpath/taxon.stat.res4plot.xls $outputpath/taxon.stat.res.pdf --width 10 --left_space 0.35 --format pdf --dpi 600};

#system qq{rm -rf taxon.stat.for.lefse.xls};
#system qq{rm -rf taxon.stat.res4plot.xls};
#system qq{rm -rf taxon.stat.format.xls};
#system qq{rm -rf taxon.stat.res.xls};
#system qq{rm -rf taxon.stat.res1.xls};
