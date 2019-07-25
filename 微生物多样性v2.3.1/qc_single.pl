#!/usr/bin/perl

BEGIN{
	push @INC,"/home/genesky/pipeline/metagenome_16s_18s_its_sequencing/v2.3.1/package";
}

$|=1;
use strict;
use warnings;
use config;
use quality_control_single;

my ($config) = @ARGV;
die "perl $0 config.txt \n" if scalar @ARGV != 1;

my $metadata = config::read_config($config);
my $base     = config::base_config();

my $suffix = "="  x 30;
sub get_time{
	my $time = `date`;
	chomp $time;
	return $time;
}

print qq{$suffix}.get_time().qq{$suffix运行质量控制\n};

quality_control_single::run($metadata, $base);
