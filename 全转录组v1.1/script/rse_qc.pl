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

use package::rseqc;
use package::RNA_evaluate;

my ($config) = @ARGV;
die "Usage:perl $0 config.txt!\n" if scalar @ARGV != 1;

my $metadata = package::config::read_config($config);
my $base     = package::config::base_config();

my $suffix = "="  x 30;

my $rseqc_proc = sub {

	print package::get_time::run().qq{${suffix}运行rseqc${suffix}\n};
	package::rseqc::run($metadata, $base);

};

my $rna_proc = sub {

	print package::get_time::run().qq{${suffix}运行RNA_evaluate${suffix}\n};
	package::RNA_evaluate::run($metadata, $base);

};

package::parallel::run_sub($rseqc_proc, $rna_proc);













