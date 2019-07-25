#!/usr/bin/env perl

BEGIN{
	push @INC,"/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/";
}

$|=1;
use strict;
use warnings;

use Getopt::Long;
use package::config;
use package::get_time;
use package::format;
use package::parallel;

use package::quality_control;

my ($force_sample,    
	$thread_qc,  
	$adapter,
	$help);

GetOptions(
	'force_sample=s'  => \$force_sample,
	'thread_qc=i'     => \$thread_qc,
	'adapter=s'       => \$adapter,
 	'help|h!'         => \$help
);

if ($help) {
	usage();
	exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl $0 config.txt
    [ Options ]
    -force_sample  force run sample, eg [sampleA, sampleB]          
    -thread_qc     thread of qc  [option]
    -adapter       adapter sequence [option:truseq, nextera]
    -help|-h       print help message

    [ Arguments ]
    config.txt
EOF
	print $help;
	exit;
}


my ($config) = @ARGV;
die "Usage:perl $0 config.txt!\n" if scalar @ARGV != 1;

my $metadata = package::config::read_config($config);
my $base     = package::config::base_config();

$base->{'thread_qc'}    = $thread_qc    if $thread_qc;
$base->{'force_sample'} = $force_sample if $force_sample;

my $suffix = "="  x 30;

print package::get_time::run().qq{${suffix}运行质量控制${suffix}\n};
package::quality_control::run($metadata, $base, $adapter);












