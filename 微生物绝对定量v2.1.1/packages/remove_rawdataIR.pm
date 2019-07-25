package remove_rawdataIR;

# use strict;
use Parallel::ForkManager;

sub run
{
	my $metadata    = shift;
	my $base        = shift;

	my $inputpath   = qq{$metadata->{raw_data}};
	my $intm_result = qq{$metadata->{intm_result}};
	my $projectpath = qq{$metadata->{result}/NO_IR_rawdata};
	my $PCpath      = qq{$metadata->{result}/Positive_Control};
	my $renamelist  = qq{$metadata->{rename_sample}};
	my $IR_readID   = qq{$metadata->{intm_result}/tmp/IR_Detection/readID_IR.txt};

	my $scriptdir   = qq{$base->{scriptdir}/sub_script};

	system qq{mkdir -p $projectpath} if not -d $projectpath;


	if (-e $renamelist){
		Rename_id($IR_readID, $renamelist, $projectpath);
	}else{
		system qq{cp $IR_readID $projectpath};
	}	


	###开始删除
	Remove_IRseq($intm_result, "$inputpath", "$projectpath/readID_IR.txt", $projectpath, $PCpath);

	# system qq{touch $intm_result/remove_rawdataIR.finish};
	print qq{rawdata spike-in 序列删除已经运行完成\n};
}

#######============== 子程序 ==============

sub Remove_IRseq
{

	my $tmpPath  = shift;
	my $inpath   = shift;
	my $readID   = shift;
	my $NOIR_dir = shift;
	my $PC_dir   = shift;
	# print "$inpath\n$readID\n$NOIR_dir";
	# return if -e qq{$tmpPath/remove_rawdataIR.finish};
	
	###如果存在PC样本
	my %hash_check;
	if (-e "$PC_dir/sample_groups.xls") {
		open PC, "$PC_dir/sample_groups.xls";
		while (<PC>) {
			my ($pcname, $tmp2) = split /\t/, $_;
			$hash_check{$pcname} = "---";			
		}
		close PC;
	}

	###读入spike-in seq ID;
	my %Readid;
	open READ_ID, $readID;
	while (<READ_ID>) {
	  chomp;
	  my ($oldname) = $_ =~ /barcodelabel\=(.+)\;/;
	  next if exists $hash_check{$oldname};			##排除阳性对照样本

	  my ($readid, $tmp1) = split /\;/, $_;
	  $readid =~ s/_/:/g;
	  my @arr = split /:/, $readid;

	  $Readid{$oldname}{"1"}{$readid} = "---";
	  $Readid{$oldname}{"2"}{$readid} = "---";
	}
	close READ_ID;


	####统计
	open STAT, ">$NOIR_dir/IR_statistics.xls";
	print STAT "sample\traw_counts\tremain_counts\tremove_counts\tIR_percent\n";

	my $MAX_PROCESSES = 18;
	my $pm = Parallel::ForkManager->new($MAX_PROCESSES);

	foreach my $sample (keys %Readid){

	    my $pid = $pm->start and next;

	  	my $R1_filename = "$sample\_R1";
	  	my $R2_filename = "$sample\_R2"; 

	  	###删除数据并返回统计结果
	  	my $stat1 = remove_seqIR_sub("$inpath/$R1_filename\.fastq.gz", "$NOIR_dir/$R1_filename\.fastq", $R1_filename, %{$Readid{$sample}{"1"}});
	  	my $stat2 = remove_seqIR_sub("$inpath/$R2_filename\.fastq.gz", "$NOIR_dir/$R2_filename\.fastq", $R2_filename, %{$Readid{$sample}{"2"}});

	  	###输出统计结果
	  	print STAT "$stat1";
	  	print STAT "$stat2";

	  	###打包数据
	  	system qq{rm $NOIR_dir/$R1_filename\.fastq.gz} if -e "$NOIR_dir/$R1_filename\.fastq.gz";
	  	system qq{rm $NOIR_dir/$R2_filename\.fastq.gz} if -e "$NOIR_dir/$R2_filename\.fastq.gz";

	  	system qq{gzip $NOIR_dir/$R1_filename\.fastq};
	  	system qq{gzip $NOIR_dir/$R2_filename\.fastq};

	  	$pm->finish;

	}
	$pm->wait_all_children; 
	close STAT;

}


sub remove_seqIR_sub
{

  my ($raw_fastq, $out, $sample, %hash, $dirction) = @_;
  my $all = 0;       ##用于统计原始序列数
  my $remain = 0;    ##用于统计保留序列数
 
  # print "$raw_fastq";
  open FASTQ,  qq{zcat "$raw_fastq" |} or die "can not open $raw_fastq\n"; ##zcat,解压读取；zcat，获取压缩文件file内容；|,管道符；
  open OUTPUT, qq{>$out} or die "Can not write to file: \n";
  while (my $id = <FASTQ>) {
	  my $seq     = <FASTQ>;
	  my $comment = <FASTQ>;
	  my $quaity  = <FASTQ>;##一次性读取四行；

	  my ($readid)  = $id =~ /^\@(\S+)\s/;
	  $all++;
	  next if exists $hash{$readid};
	  print OUTPUT qq{$id$seq$comment$quaity};
	  $remain++;
  }
  close FASTQ;
  close OUTPUT;
  my $remove = $all - $remain;   ###统计删除的IR序列数
  my $IR_per = $remove/$all;
  my $stat = $sample."\t".$all."\t".$remain."\t".$remove."\t".$IR_per."\n";
  return $stat;

}


sub Rename_id
{

	my $readid  = shift;
	my $list    = shift;
	my $outpath = shift;


	my %hash = ();
	open NAME, $list;  ##两列，第一列原来的名字，第二列新名
	while (<NAME>) {
		$_ =~ s/[\r\n]//g;
		my ($name1,$name2) = split /\t/, $_;
		$hash{$name2} = $name1;
	}
	close NAME;

	my %tmp = ();
	open FASTA, $readid;
	open NEW, ">$outpath/readID_IR.txt";
	while (my $sample = <FASTA>) {
		my @arry = split /=/, $sample;
		my ($samplename) = $sample =~ /barcodelabel\=(.+)\;/;
		my ($tmp1, $otuid) = split /\t/, $arry[1];
		
		if (exists $hash{$samplename}){
			print NEW qq{$arry[0]=$hash{$samplename};\t$otuid};
			$tmp{$samplename} = $hash{$samplename};
		}

	}
	close FASTA;
	close NEW;

}

1;
