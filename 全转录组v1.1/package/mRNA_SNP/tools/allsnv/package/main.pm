package package::main;
use strict;
use warnings;

sub readLibrary
{
	my $file=shift @_;
	my %library=();
	print "Reading $file...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$q,$rs,$freq,$ref,$alt,$chr,$pos,$gene,$region,$function,$hgvs,$seq5,$seq3,$blast)=split /\t/,$line;
		$library{$title}{$alt}{"SNP Calling Quality"}=$q;
		$library{$title}{$alt}{"SNP ID"}=$rs;
		$library{$title}{$alt}{"Freq_Alt (1000g)"}=$freq;
		$library{$title}{$alt}{"Ref Allele"}=$ref;
		$library{$title}{$alt}{"Alt Allele"}=$alt;
		$library{$title}{$alt}{"Chrs"}=$chr;
		$library{$title}{$alt}{"Position"}=$pos;
		$library{$title}{$alt}{"Position2"}=(split /\-/,$pos)[0];
		$library{$title}{$alt}{"Gene"}=$gene;
		$library{$title}{$alt}{"Gene Region"}=$region;
		$library{$title}{$alt}{"Function"}=$function;
		$library{$title}{$alt}{"Predicted Protein Variants"}=$hgvs;
		$library{$title}{$alt}{"5' FLANKING SEQUENCE"}=$seq5;
		$library{$title}{$alt}{"3' FLANKING SEQUENCE"}=$seq3;
		$library{$title}{$alt}{"HOMOLOGY HITS"}=$blast;
	}
	close FILE;
	print "OK\n";
	return %library;
}



sub readAnnotation
{
	my $file=shift @_;
	my %annotation=();
	print "Reading $file...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$type,$content)=split /\t/,$line,4;
		if($type eq "HGMD_site" or $type eq "HGMD_site_class"){
		   next if(exists($annotation{$title}{$alt}{'HGMD_site_class'}) and $annotation{$title}{$alt}{'HGMD_site_class'}=~/DM/);
		}elsif($type eq 'Strand Orientation'){# 表头名称发生变动
		   $annotation{$title}{$alt}{"Gene Strand Orientation"}=$content;
		   next;
		}
		$annotation{$title}{$alt}{$type}=$content;
	}
	close FILE;
	print "OK\n";
	return %annotation;
}


sub readAnalysis
{
	my $file=shift @_;
	my %analysis=();
	print "Reading $file...";
	open FILE,$file;
	while(my $line=<FILE>){
		$line=~ s/[\r\n]//g;
		my ($title,$alt,$genoQ,$s0,$s1,$s2)=split /\t/,$line;
		$analysis{$title}{$alt}{"Genotyping Quality"}=$genoQ;
		$analysis{$title}{$alt}{"Mutation 0"}=$s0;
		$analysis{$title}{$alt}{"Mutation 1"}=$s1;
		$analysis{$title}{$alt}{"Mutation 2"}=$s2;
	}
	close FILE;
	print "OK\n";
	return %analysis;
}


1