#!/usr/bin/perl

BEGIN{
	push @INC,"/home/genesky/pipeline/metagenome_16s_18s_its_sequencing/v2.3.1/package";
}

$|=1;
use strict;
use warnings;
use Getopt::Long;
use config;
use quality_control;
use otu_fungene_noresample;
use otu_fungene_resample;
use community;
use alphaDiversity;
use betaDiversity;
use environmentFactors;

my ($skip,
	$only,
	$help);

GetOptions(
	'skip=s'          => \$skip,
	'only=s'          => \$only,
	'help|h!'         => \$help
);

Usage() if ( $help or (scalar @ARGV != 1) );

sub Usage{
my $help ="
Usage: perl $0  config.txt

    [ Options ]
    --------------- steps  ----------------
    -skip          skip steps, eg [4,5,6,7,8]
    -only          only steps, eg [1,2,3]
    -help|-h       print help message
	
    [ Arguments ]
    --------------- argus  ----------------
    config.txt

    [ step list ]
    ---------------  list  ----------------
    1 质量控制
    2 OTU聚类注释
    3 Resample抽平分析
    4 Community群落分析
    5 Alpha多样性分析
    6 Beta多样性分析
    7 功能预测分析
    8 环境因子分析
\n";
print $help;
exit;
}

if ($only and $skip) {
	die "you can just pick only one of options [only / skip]\n";
}

my @run_steps = ();
if ($only) {
	foreach my $x (split /,/, $only) {
		push @run_steps, $x;
	}
} else {
	my @skips = split /,/, $skip if defined ($skip);
	foreach my $x (1..8) {
		next if $x ~~ @skips;
		push @run_steps, $x;
	}
}

$SIG{INT} = sub {
    print "you has enter the CTRL+C keys,now exit\n";
    exit 0;
};

my ($config) = @ARGV;
#die "perl $0 config.txt \n" if scalar @ARGV != 1;

my $metadata = config::read_config($config);
my $base     = config::base_config();
my $type     = qq{$metadata->{type}};
			 
my $suffix = "="  x 30;
sub get_time{
	my $time = `date`;
	chomp $time;
	return $time;
}

if (exists $metadata->{'raw_data'} and (1 ~~ @run_steps)){
	print qq{$suffix}.get_time().qq{$suffix运行质量控制\n};
	quality_control::run($metadata, $base);
}

if (2 ~~ @run_steps){ 
	print qq{$suffix}.get_time().qq{$suffix运行fungene聚类注释分析\n};
	otu_fungene_noresample::run($metadata, $base);
}

if (exists $metadata->{'group'}){

	if (3 ~~ @run_steps){
		print qq{$suffix}.get_time().qq{$suffix运行Resample抽平分析\n};
		otu_fungene_resample::run($metadata, $base);
	}
	
	if (4 ~~ @run_steps){
		print qq{$suffix}.get_time().qq{$suffix运行Community群落分析\n};
		community::run($metadata, $base);
	}
	
	if (5 ~~ @run_steps){
		print qq{$suffix}.get_time().qq{$suffix运行Alpha多样性分析\n};
		alphaDiversity::run($metadata, $base);
	}
	
	if (6 ~~ @run_steps){
		print qq{$suffix}.get_time().qq{$suffix运行Beta多样性分析\n};
		betaDiversity::run($metadata, $base);
	}
	
	if (7~~ @run_steps){
		if (uc $type eq "16S"){
			print qq{$suffix}.get_time().qq{$suffix运行PICRUSt功能分析\n};
			picrust::run($metadata, $base);
		}
		if (uc $type =~ /ITS/){
			print qq{$suffix}.get_time().qq{$suffix运行fungulid分析\n};
			fungulid::run($metadata, $base);
		}
	}
}

if (exists $metadata->{'env'} and (8 ~~ @run_steps)){
	print qq{$suffix}.get_time().qq{$suffix运行ENV分析\n};
	environmentFactors::run($metadata, $base);
}