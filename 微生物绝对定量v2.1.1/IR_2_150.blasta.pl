use strict;
use Cwd 'abs_path';
use Parallel::ForkManager;

my ($inpath, $outpath) = @ARGV;
$inpath = abs_path($inpath);
$outpath = abs_path($outpath);

####IR比对结果目录
my $inter_ref_db = "/home/zhengy/bin/modules/database/SpikeIN_database/SpikeIN/16S_SPIKEIN.fasta";    #####若内参有变需要改动
my $map_dir = "$outpath/IR_2_150/IR_map";            ###存放比对结果
system qq{mkdir -p $map_dir} if not -d $map_dir;

my $fas_dir = "$outpath/IR_2_150/raw_data_fasta";    ###存放原始数据
system qq{mkdir -p $fas_dir} if not -d $fas_dir;

my @samples = `ls $inpath | grep "_R1.fastq.gz" | sed 's/_R1.fastq.gz//g'`;
my %raw_stat;
open STA, ">$map_dir/rowdata.stat.xls";
print STA "sample\traw_reads\n";

my $MAX_PROCESSES = 18;
my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

for my $sample (@samples){
	my $pid = $pm->start and next;
	chomp $sample;
	my $gz1 = "$inpath/$sample\_R1.fastq.gz";
	my $fasta1 = "$fas_dir/$sample\_R1.fastq";
	open GZ, qq{zcat $gz1 |};
	open FA, ">$fasta1";
	my $count = 0;
	while (my $id = <GZ>) {
		my $seq = <GZ>;
		my $jia = <GZ>;
		my $q   = <GZ>;
		chomp $id;
		my @tmp = split /\s/, $id;
		$tmp[0] =~ s/\@/>/g;
		$tmp[1] = "barcodelabel\=".$sample;
		my $idnew = "$tmp[0]\;$tmp[1]";
		print FA "$idnew\n$seq";
		$count++;

	}
	close GZ;
	close FA;
	$raw_stat{$sample} = $count;
	print STA "$sample\t$count\n";

	###IR比对
	my $IR_map = "$map_dir/$sample.IR.map.xls";
	system qq{/home/genesky/software/blast+/2.7.1/bin/blastn -query $fasta1 -db $inter_ref_db -evalue 1e-40 -outfmt 6 -max_target_seqs 100000000 -perc_identity 96 -out $IR_map 2> $map_dir/$sample.log};

	$pm->finish;
}
$pm->wait_all_children; 
close STA;



###统计结果
my $count_dir = "$outpath/IR_2_150/L-120_M-0.96";
system qq{mkdir -p $count_dir} if not -d qq{$count_dir};

####提取内参名字
open DBIR, $inter_ref_db;
my @dbir;
while (<DBIR>) {
	chomp;
	next if !/^>/;
	my ($irid) = $_ =~ /^>(\w+?)\s/;
	push @dbir, $irid;

}
close DBIR;
@dbir = sort @dbir;


my $stat = "$map_dir/rowdata.stat.xls";
open STAT, $stat;
open PER, ">$count_dir/IR_persent.L-120.M-0.96.xls";
open IRC, ">$count_dir/IR_count.L-120.M-0.96.xls";
open ALL, ">$count_dir/IR_and_Raw.xls";
print PER "sample\t".(join "\t", @dbir)."\n";   ###输入表头
print IRC "sample\t".(join "\t", @dbir)."\n";
print ALL "sample\traw_reads\tIR_reads\tIR(%)\n";

while (<STAT>) {
	chomp;
	next if /^sample/;
	my ($sample, $count) = split /\t/, $_;
	my $IR_map = "$map_dir/$sample.IR.map.xls";

	open IRM, $IR_map;
	###统计个数
	my %count_ir;
	while (<IRM>) {
		chomp;
		my @tmp = split /\t/,$_;
		$count_ir{$tmp[1]}++ if ($tmp[2] >= 96 && $tmp[3] >= 120);   ####匹配率大于96%，序列长度大于120
		
	}
	close IRM;
	

	###计算
	my $irall = 0;
	my @result_c;
	my @result_p;
	for my $key (0..$#dbir){
		my $ircount =  $count_ir{$dbir[$key]};
		$ircount = 0 if ($ircount eq "");
		my $irper   =  $ircount*100/$count;
		$irper      =  sprintf "%.2f",$irper;

		push @result_c, $ircount;
		push @result_p, $irper;
		$irall = $irall + $ircount;
	}
	my $irall_per = $irall*100/$count;
	$irall_per    = sprintf "%.2f", $irall_per;

	print PER $sample."\t".(join "\t", @result_p)."\n";   ##
	print IRC $sample."\t".(join "\t", @result_c)."\n";   ##
	print ALL "$sample\t$count\t$irall\t$irall_per\n";

}
close STAT;
close PER;
close IRC;
close ALL;

# system qq{rm -rf $map_dir};
system qq{rm -rf $fas_dir};

