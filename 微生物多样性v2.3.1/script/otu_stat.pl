#!/usr/bin/perl
#use strict;
#use warnings;

my ($otufile,$outputpath) = @ARGV;

die "perl $0 otufile outputpath\n" if scalar @ARGV != 2;

my ($dissa) = $otufile =~ /otu.tax.(.*).xls/;
my @taxonclass = ("superkingdom", "phylum", "class", "order", "family", "genus", "species");
my @tags = ("k", "p", "c", "o", "f", "g", "s");
my $index = 0;
my %hash;

open INPUT, "$otufile";
open OUTPUT, ">$outputpath/otu.tax.$dissa.for.stat.xls";

while (<INPUT>) {

	$_ =~ s/[\r\n]//g;
	if (/Taxonomy/) {
	
		my @title = split /\t/;
		($index) = grep { $title[$_] eq "Taxonomy" } 0..$#title;   #提取符合条件的下标，0..$#arr 从0到@arr的最后一个，等价于@arr；
		print OUTPUT "Taxonomy"."\t".(join "\t", @title[1..($index-2)])."\t"."Abundance"."\n";
		
	}
	else{
	
		my @arr   = split /\t/;
		my @data  = @arr[1..($index-1)];  #存储数据
		my @taxon = (@arr[($index+1)..$#arr]);  #存储分类信息

		my ($taxn_no) = grep { $taxon[$_] eq "" } 0..$#taxon;  #返回第一个为空值的下标
		
		my @taxonomy;
		
		# 注释层无中断,取全部注释内容
		if ($taxn_no eq "") {
			for my $s (0..$#tags) {
				$taxonomy[$s] = $tags[$s]."__".$taxon[$s];
			}
		} 
		# 注释层中断,取中断位置前的注释内容
		else {
			for my $s (0..($taxn_no-1)) {
				$taxonomy[$s] = $tags[$s]."__".$taxon[$s];
			}
		}
		my $str = join "|", @taxonomy;

		#如果存在就累加
		if (exists $hash{$str}) {
			for my $num (0..$#data) { 
				$hash{$str}[$num] += $data[$num];
			}
		}
		else {
			$hash{$str} =\@data;
		}
		
	}

}
close INPUT;

foreach my $key (sort keys %hash){
	print OUTPUT $key."\t".(join "\t", @{$hash{$key}})."\n";
	# print $key."\t".(join "\t", @{$hash{$key}})."\n";
}
close OUTPUT;


open FORSTAT, "$outputpath/otu.tax.$dissa.for.stat.xls";
open STAT, ">$outputpath/otu.tax.$dissa.stat.xls";

my %stats;
while(<FORSTAT>){

	chomp;
	next if /^\s+$/;
	
	if (/taxonomy/){
	
		print STAT qq{$_\n};
	
	}else{
	
		my @st = split /\t/;
		my @Tax = split /\|/, $st[0];
		
		if ($#Tax < 1){       # Only kingdom
			
			if (exists $stats{$Tax[0]}){
			
				for my $i (0..($#st-1)) {
				
					$stats{$Tax[0]}->[$i] = $stats{$Tax[0]}->[$i] + $st[$i+1];
					
				}
			
			}else{
			 
				@{$stats{$Tax[0]}} = @st[1..$#st];
			
			}
		
		}else{
		
			for my $t (0..$#Tax){
			
				my $test = join "|",@Tax[0..$t];      # 存在上层注释则数据合并
			
				if (exists $stats{$test}){    
			
					for my $i (0..($#st-1)) {
					
						$stats{$test}->[$i] = $stats{$test}->[$i] + $st[$i+1];
						
					}
			
				}else{
			
					@{$stats{$test}} = @st[1..$#st];
			
				}
				
			}
			
		}
		
	}

}
close FORSTAT;

foreach my $s (sort keys %stats){

	print STAT $s."\t".(join "\t",@{$stats{$s}})."\n";
	
}
close STAT;
