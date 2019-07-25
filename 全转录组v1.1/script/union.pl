#!/usr/bin/env/perl

BEGIN{
	push @INC, "/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/";
}

use Getopt::Long;
use package::config;
use package::get_time;
use package::parallel;
# association
use package::association::mRNA_miRNA;
use package::association::lncRNA_miRNA;
use package::association::circRNA_miRNA;
use package::association::lncRNA_cis;
use package::association::lncRNA_trans;
# ceRNA
use package::ceRNA::mRNA_miRNA_lncRNA;
use package::ceRNA::mRNA_miRNA_circRNA;

my ($skip, 
	$only,      
	$help, 
	$list);

GetOptions(
	'skip=s'          => \$skip,
	'only=s'          => \$only,
 	'help|h!'         => \$help,
 	'list!'           => \$list
);

if ($list) {
my $help =<<EOF;
1 mRNA和miRNA相互作用分析
2 lncRNA和miRNA相互作用分析
3 circRNA和miRNA相互作用分析
4 lncRNA和mRNA相互作用分析(Cis和Trans分析)
5 mRNA_miRNA_lncRNA ceRNA分析
6 mRNA_miRNA_circRNA ceRNA分析
EOF
print qq{$help\n};
exit;	
}

if ($help) {
	usage();
	exit();
}

sub usage
{
my $help =<<EOF;
Usage: perl $0  config.txt
    [ Options ]
    --------------- steps  ----------------
    -skip          skip steps, eg [3,4,5]
    -only          only steps, eg [1,2,3]
    -help|-h       print help message
    -list          list steps 

    [ Arguments ]
    config.txt
EOF
	print $help;
	exit;

}

if ($only and $skip) {
	die "you can just pick only one of options [only / skip]\n";
}

my ($config) = $ARGV[$#ARGV];
my $metadata = package::config::read_config($config);
my $base     = package::config::base_config();

my @run_steps = ();
if ($only) {
	foreach my $x (split /,/, $only) {
		push @run_steps, $x;
	}
} else {
	my @skips = split /,/, $skip;
	foreach my $x (1..6) {
		next if $x ~~ @skips;
		push @run_steps, $x;
	}
}

$SIG{INT} = sub {
    print "you has enter the CTRL+C keys,now exit\n";
    exit 0;
};


my $suffix = "="  x 30;

if (1 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA和miRNA相互作用分析${suffix}\n};
	package::association::mRNA_miRNA::run($metadata, $base);
}
if (2 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行lncRNA和miRNA相互作用分析${suffix}\n};
	package::association::lncRNA_miRNA::run($metadata, $base);
}
if (3 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行circRNA和miRNA相互作用分析${suffix}\n};
	package::association::circRNA_miRNA::run($metadata, $base);
}
if (4 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行lncRNA Cis分析${suffix}\n};
	package::association::lncRNA_cis::run($metadata, $base);
	print package::get_time::run().qq{${suffix}运行lncRNA Trans分析${suffix}\n};
	package::association::lncRNA_trans::run($metadata, $base);
}
if (5 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA_miRNA_lncRNA ceRNA分析${suffix}\n};
	package::ceRNA::mRNA_miRNA_lncRNA::run($metadata, $base);
}
if (6 ~~ @run_steps) {
	print package::get_time::run().qq{${suffix}运行mRNA_miRNA_circRNA ceRNA分析${suffix}\n};
	package::ceRNA::mRNA_miRNA_circRNA::run($metadata, $base);
}






