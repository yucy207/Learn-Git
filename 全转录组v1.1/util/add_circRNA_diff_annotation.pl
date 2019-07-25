#！usr/bin/perl
use strict;
use warnings;

die "======= Add_CircRNA_diff_annotation =========
Usage: 
    perl add_circRNA_diff_annotation.pl summary_file exp_file diff_file organ outfile
        parameters ->
            summary_file: [file -> always project/circRNA/circRNA_predict/result/Summary.xls];
                exp_file: [file -> always circRNA.expression.add.annotation.xls];
               diff_file: [file -> always circRNA.different.expression.xls];   
                   organ: [organ-> always like human,mouse,fly,zebrafish];       
                 outfile: [file -> output file circRNA.different.expression.add.annotation.xls];\n" if (@ARGV != 5);

my($sum, $exp, $diff, $organ, $out) = @ARGV;

### about Summary.xls ###
my %sum_hash = ();
open(SUM, $sum) or die "cannot open Summary.xls\n";
my $head = <SUM>;
chomp $head;
my @head = split/\t/,$head;
pop @head;
shift @head; 
my $sum_head = join("\t",@head);
while(my $line = <SUM>){
	chomp $line;
	my @arr = split/\t/,$line;
	my $id  = $arr[0];
	$id =~s/\|/\-/;
	pop @arr;
	shift @arr;
	#print "$line\n";
	$sum_hash{$id} = join("\t",@arr);
}
close SUM;

### about circRNA.expression.add.annotation.xls ###
my %exp_hash = ();
open(EXP, $exp) or die "cannot open circRNA.expression.add.annotation.xls\n";
my $head = <EXP>;
chomp $head;
my @head = split/\t/,$head;
my @Organ = ("human","mouse","fly","zebrafish");
my $index;
if(grep {$_ eq $organ} @Organ){
	$index = $#head - 5;
}else{
	$index = $#head - 1;
}
my $exp_head = join("\t",@head[$index..$#head]);
while(my $line = <EXP>){
	chomp $line;
	my @arr = split/\t/,$line;
	my $id  = $arr[0]; 
	$exp_hash{$id} = join("\t",@arr[$index..$#arr]);
}
close EXP;

#### add ####
open(IN, $diff) or die "cannot open it\n";
open(OUT, ">$out") or die "cannot make it\n";
my $title = <IN>;
chomp $title;
my $title1 = $title."\t".$sum_head."\t".$exp_head;
print OUT "$title1\n";
while(my $line = <IN>){
	chomp $line;
	my $id = (split/\t/,$line)[0];
	my $res;
	if(exists $sum_hash{$id}){$res = $line."\t".$sum_hash{$id};}
	if(exists $exp_hash{$id}){$res = $res."\t".$exp_hash{$id};}
	print OUT "$res\n";
}
close IN;
close OUT;