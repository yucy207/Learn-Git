#!/usr/bin/env perl


BEGIN{
	push @INC,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/";
}

$|=1;
use strict;
use warnings;
use package::config;
use package::get_time;
use package::parallel;

use package::assemble_isoform;
use package::merge_transcript;

my ($config) = @ARGV;
die "Usage:perl $0 config.txt!\n" if scalar @ARGV != 1;

my $metadata = package::config::read_config($config);
my $base     = package::config::base_config();

my $suffix = "="  x 30;

print package::get_time::run().qq{${suffix}运行转录本组装${suffix}\n};
package::assemble_isoform::run($metadata, $base);

print package::get_time::run().qq{${suffix}运行合并转录本${suffix}\n};
package::merge_transcript::run($metadata, $base);

