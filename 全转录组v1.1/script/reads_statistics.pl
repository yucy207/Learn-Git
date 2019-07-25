#!/usr/bin/env perl

BEGIN{
	push @INC,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/";
}

$|=1;
use strict;
use warnings;
use package::config;
use package::get_time;
use package::format;
use package::parallel;

use package::reads_evaluate;



my ($config) = @ARGV;
die "Usage:perl $0 config.txt!\n" if scalar @ARGV != 1;

my $metadata = package::config::read_config($config);
my $base     = package::config::base_config();

my $suffix = "="  x 30;

print package::get_time::run().qq{${suffix}运行数据统计${suffix}\n};
package::reads_evaluate::run($metadata, $base);












