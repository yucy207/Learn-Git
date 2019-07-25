package package::time;
use strict;
use warnings;

sub ymd{
	my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
	$year+=1900;
	$mon+=1;
	return "$year-$mon-$day";
}

sub hms{
	my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
	return "$hour:$min:$sec";
}


1