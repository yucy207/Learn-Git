package package::readme;
use strict;
use warnings;
use Encode;

sub run{
	my ($sheet,$format,$file)=@_;
	$sheet->set_row(0, 65);
	$sheet->set_column('A:A', 35);
	$sheet->set_column('B:B', 25);
	$sheet->set_column('C:C', 110);
	my $row=0;
	open FILE,$file;
	while(<FILE>){
		$_=~ s/\^/\n/g;
		my @split_line=split /\t/,$_;
		my $col=0;
		foreach(@split_line)
		{
			my $text=decode("gb2312",$_);
			if($row==0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme1'});}
			if($row==0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme2'});}
			if($row==0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme2tmp'});}
			if($row>0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme3'});}
			if($row>0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme4'});}
			if($row>0 and $col==2){$sheet->write($row,$col,$text,$format->{'readme5'});}
			$col++;
		}
		$row++;
	}
	close FILE;
}


1