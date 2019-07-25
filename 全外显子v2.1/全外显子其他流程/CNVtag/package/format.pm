package package::format;
use strict;
use warnings;


sub addformat{
	my ($workbook)=shift @_;
	my %format=();
	$format{'title'} = $workbook->add_format();
	$format{'title'} ->set_align('center');
	$format{'title'} ->set_align('vcenter');
	$format{'title'} ->set_size(12);
	$format{'title'} ->set_font("Times New Roman");
	$format{'title'} ->set_border();
	$format{'title'} ->set_bg_color("yellow");
	$format{'title'} ->set_color("black");
	
	$format{'normal'} = $workbook->add_format();
	$format{'normal'} ->set_align('center');
	$format{'normal'} ->set_align('vcenter');
	$format{'normal'} ->set_size(12);
	$format{'normal'} ->set_font("Times New Roman");
	$format{'normal'} ->set_border();
	
	$format{'normalgrey'} = $workbook->add_format(color => "#C0C0C0");
	$format{'normalgrey'} ->set_align('center');
	$format{'normalgrey'} ->set_align('vcenter');
	$format{'normalgrey'} ->set_size(12);
	$format{'normalgrey'} ->set_font("Times New Roman");
	$format{'normalgrey'} ->set_border();

	$format{'redgrey'} = $workbook->add_format(color => "#C0C0C0");
	$format{'redgrey'} ->set_align('center');
	$format{'redgrey'} ->set_align('vcenter');
	$format{'redgrey'} ->set_size(12);
	$format{'redgrey'} ->set_font("Times New Roman");
	$format{'redgrey'} ->set_border();		
	$format{'redgrey'} ->set_bg_color("red");

	$format{'red'} = $workbook->add_format();
	$format{'red'} ->set_align('center');
	$format{'red'} ->set_align('vcenter');
	$format{'red'} ->set_size(12);
	$format{'red'} ->set_font("Times New Roman");
	$format{'red'} ->set_border();		
	$format{'red'} ->set_bg_color("red");

	$format{'greengrey'} = $workbook->add_format(color => "#C0C0C0");
	$format{'greengrey'} ->set_align('center');
	$format{'greengrey'} ->set_align('vcenter');
	$format{'greengrey'} ->set_size(12);
	$format{'greengrey'} ->set_font("Times New Roman");
	$format{'greengrey'} ->set_border();		
	$format{'greengrey'} ->set_bg_color("red");

	$format{'green'} = $workbook->add_format();
	$format{'green'} ->set_align('center');
	$format{'green'} ->set_align('vcenter');
	$format{'green'} ->set_size(12);
	$format{'green'} ->set_font("Times New Roman");
	$format{'green'} ->set_border();		
	$format{'green'} ->set_bg_color("green");

	
	$format{'small'} = $workbook->add_format();
	$format{'small'} ->set_align('vcenter');
	$format{'small'} ->set_size(10);
	$format{'small'} ->set_font("Times New Roman");
	$format{'small'} ->set_border();
	
	$format{'seq'} = $workbook->add_format();
	$format{'seq'} ->set_align('vcenter');
	$format{'seq'} ->set_size(11);
	$format{'seq'} ->set_font("Courier New");
	$format{'seq'} ->set_border();
	
	$format{'left'} = $workbook->add_format();
	$format{'left'} ->set_align('vcenter');
	$format{'left'} ->set_size(12);
	$format{'left'} ->set_font("Times New Roman");
	$format{'left'} ->set_border();
	
	$format{'orange'} = $workbook->add_format();
	$format{'orange'} ->set_align('vcenter');
	$format{'orange'} ->set_size(12);
	$format{'orange'} ->set_font("Times New Roman");
	$format{'orange'} ->set_bg_color("#fac090");
	$format{'orange'} ->set_border();

	$format{'bold'} = $workbook->add_format( bold => 1 );
	$format{'blue'} = $workbook->add_format( color => "#538ed5" );
	$format{'boldblue'} = $workbook->add_format( bold => 1, color => "#538ed5" );
	
	$format{'readme1'} = $workbook->add_format();
	$format{'readme1'}->set_align('center');
	$format{'readme1'}->set_align('vcenter');
	$format{'readme1'}->set_bold();
	$format{'readme1'}->set_size(14);
	$format{'readme1'}->set_font("Times New Roman");
	$format{'readme1'}->set_border();

	$format{'readme2'} = $workbook->add_format();
	$format{'readme2'}->set_align('vcenter');
	$format{'readme2'}->set_bold();
	$format{'readme2'}->set_size(14);
	$format{'readme2'}->set_font("Times New Roman");

	$format{'readme2tmp'} = $workbook->add_format();
	$format{'readme2tmp'}->set_right();

	$format{'readme3'} = $workbook->add_format();
	$format{'readme3'}->set_align('center');
	$format{'readme3'}->set_align('vcenter');
	$format{'readme3'}->set_bold();
	$format{'readme3'}->set_size(11);
	$format{'readme3'}->set_font("Times New Roman");
	$format{'readme3'}->set_border();

	$format{'readme4'} = $workbook->add_format();
	$format{'readme4'}->set_align('vcenter');
	$format{'readme4'}->set_bold();
	$format{'readme4'}->set_size(11);
	$format{'readme4'}->set_font("Times New Roman");
	$format{'readme4'}->set_border();

	$format{'readme5'} = $workbook->add_format();
	$format{'readme5'}->set_align('vcenter');
	$format{'readme5'}->set_size(11);
	$format{'readme5'}->set_font("Times New Roman");
	$format{'readme5'}->set_border();
	$format{'readme5'}->set_text_wrap();

	return %format;
}

1