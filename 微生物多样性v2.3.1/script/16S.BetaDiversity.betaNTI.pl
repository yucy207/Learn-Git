#!/usr/bin/perl
use Cwd 'abs_path';
use Parallel::ForkManager;

my ($otufile, $groupfile, $trefile, $outputPath) = @ARGV;
die "Usage: perl $0 otufile groupfile trefile outputPath\n" if scalar @ARGV != 4;

my $num = `less $groupfile | cut -f 2 | sort -u | wc -l`;
my $num_sample =  `less $groupfile | cut -f 1 | wc -l`;
my $rep = $num_sample/$num;
die "BetaNTI分析组内重复少于3\n" if $rep < 3;

$otufile = abs_path($otufile);
$groupfile = abs_path($groupfile);
$trefile = abs_path($trefile);
$outputPath = abs_path($outputPath);

my $tmpPath = qq{$outputPath/tmp};
system qq{mkdir -p $outputPath} if not -d $outputPath;
system qq{mkdir -p $tmpPath} if not -d $tmpPath;

my $id = `id`;
my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;

system qq{docker run --rm -v /home:/home -u $uid:$gid zhengy/r-3.4.3:v1 Rscript /home/zhengy/bin/modules/script/16S.BetaDiversity.format.r $otufile $groupfile $tmpPath};

system qq{cp $trefile $tmpPath/tre} if not -e qq{$tmpPath/tre};

chdir $tmpPath;

my $MAX_PROCESSES = @files;
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

my @files = `find $tmpPath -name '*.format.xls'`;

for my $file (@files){
	chomp $file;
	my $pid = $pm->start and next;
	my ($name) = $file =~ /\/([^\.\/]*)\.format.xls/;
	#system qq{/home/zhengy/panrf/software/phylocom/phylocom-4.2/phylocom comdistnt -m 0 -r 999 -n -f tre -s $file > $tmpPath/$name\.NTI.xls};
	system qq{docker run --rm -v $tmpPath:/data -u $uid:$gid chengsy_16s/phylocom phylocom comdistnt -m 0 -r 999 -n -f /data/tre -s /data/$name.format.xls > $tmpPath/$name.NTI.xls};
	
	$pm->finish;
}
$pm->wait_all_children;

system qq{docker run --rm -v /home:/home -u $uid:$gid zhengy/r-3.4.3:v1 Rscript /home/zhengy/bin/modules/script/16S.BetaDiversity.NTI.r $groupfile $tmpPath $outputPath};
system qq{rm -rf $tmpPath};