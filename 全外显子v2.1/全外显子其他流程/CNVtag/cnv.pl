use strict;
use warnings;
use Excel::Writer::XLSX;
use package::main;
use package::readtag;
use package::CNV;
use package::annotation_db;
use package::annotation;
use package::output;
use package::format;
use package::readme;
$|=1;
die "Usage perl $0 ConfigFile\n" if(@ARGV!=1);
my $ConfigFile=shift @ARGV;
my %config=package::main::readConfig($ConfigFile);
package::main::checktag(\%config);###tag���
my %hashtag=package::readtag::run(\%config);###��tag
my %hashanno=package::annotation::run(\%config,\%hashtag);
my @results=package::CNV::run(\%config,\%hashtag);###������Կ�����
my %hashCNV=%{$results[0]};
my %hashok=%{$results[1]};
package::output::run(\%config,\%hashtag,\%hashCNV,\%hashok,\%hashanno);###���
