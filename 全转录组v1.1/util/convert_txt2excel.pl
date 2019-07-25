

use Excel::Writer::XLSX;
use Encode qw/decode/;

my ($txt, $excel) = @ARGV;

open TXT, $txt or die "Can't open $txt!\n";

my $workbook  = Excel::Writer::XLSX->new($excel);
my $worksheet  = $workbook->add_worksheet();
$worksheet->set_column( 'A:Z', 20);

my %format=();
$format{'title'} = $workbook->add_format();
$format{'title'}->set_align('center');
$format{'title'}->set_align('vcenter');
$format{'title'}->set_size(12);
$format{'title'}->set_font("Times New Roman");
$format{'title'}->set_border();
$format{'title'}->set_bg_color("yellow");
$format{'title'}->set_color("black");
	
$format{'normal'} = $workbook->add_format();
$format{'normal'}->set_align('center');
$format{'normal'}->set_align('vcenter');
$format{'normal'}->set_size(12);
$format{'normal'}->set_font("Times New Roman");
$format{'normal'}->set_border();


my $row = 0;
while (<TXT>) {
	chomp;
	my @arr    = split /\t/;

	if ($row == 0) {
		$worksheet->write_row($row, 0, \@arr, $format{'title'});
	} else {
		$worksheet->write_row($row, 0, \@arr, $format{'normal'});
	}
	$row++;
}
close TXT;
$workbook->close();
