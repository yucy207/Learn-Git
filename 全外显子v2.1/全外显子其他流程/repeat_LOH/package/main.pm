package package::main;
use strict;
use warnings;

# ��ȡĬ�ϲ���
sub get_para{
    my %hashPara;

    # �����߳�������
    $hashPara{"Process"}{"repeat_LOH"}  = 20;
	
	# ��������
	$hashPara{'Human_hg19'}{'Genome'}       = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.fa";
	$hashPara{'Human_hg19'}{'Dict'}         = "/home/genesky/database/ucsc/hg19_modify/genome/hg19_modify.dict";
	$hashPara{'Human_hg19'}{'AnnovarBuild'} = "hg19";
	
	$hashPara{'Human_hg38'}{'Genome'}       = "/home/genesky/database/ucsc/hg38/genome/hg38.fa";
	$hashPara{'Human_hg38'}{'Dict'}         = "/home/genesky/database/ucsc/hg38/genome/hg38.dict";
	$hashPara{'Human_hg38'}{'AnnovarBuild'} = "hg38";

	# �������
	$hashPara{'Soft'}{'AnnovarDIR'} = "/home/genesky/software/annovar/2018Apr16";
	$hashPara{'Soft'}{'bedtools'}   = "/home/genesky/software/bedtools/2.28.0/bin/bedtools";
	$hashPara{'Soft'}{'R'}          = "/home/genesky/software/r/3.5.1/bin/R";
	$hashPara{'Soft'}{'RLib'}       = "/home/genesky/software/r/3.5.1/lib64/R/library";

    return %hashPara;
}

# ��ȡ�����ļ�
sub read_config{
	my $config_file = shift @_;
	my %hashConfig  = ();
	open CONFIG, $config_file or die "[ERR] Loss Config File, $config_file\n";
	while(my $line = <CONFIG>) {
		next if($line !~ /\w/ or $line =~ /^#/);
		$line             =~ s/[\s\t\r\n]//g;
		my ($key, $value) = split /\=/, $line;
		if($key eq 'case' or $key eq 'control') {
		    $hashConfig{$key} .= ",$value";
		}
		else {
		    $hashConfig{$key}  = $value;
		}
	}
	close CONFIG;
	return %hashConfig;
}

sub get_sample{
    my $hashConfig = shift @_;
    my @need_types = @_;
    my @samples;
    foreach my $need_type(@need_types) {
        next if(!exists($hashConfig->{$need_type}));
        foreach my $sample(split /,/, $hashConfig->{$need_type}) {
            push @samples, $sample if($sample =~ /\w/);
        }
    }
    return @samples;
}

# ����Ļ�������
sub get_input{
    my $input = <STDIN>;
    $input =~ s/[\r\n]//g;
    return $input;
}

# ��ȡʱ��
sub get_time{
    my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
    $year+=1900;
    $mon+=1;
    return "$year-$mon-$day $hour:$min:$sec";
}

# ����Ŀ¼
sub make_dir{
    my $dir = shift @_;
    mkdir $dir if(not -e $dir);
}

# �����ļ��Ƿ�Ϊ��
sub is_file_ok{
    my @files = @_;
    my $isOK = 1;
    foreach my $file(@files)
    {
        $isOK = 0 if (not -e $file or not -f $file or -s $file == 0);
    }    
    return $isOK;
}

# ����Ŀ¼�Ƿ����
sub is_dir_ok{
    my $dir = shift @_;
    my $isOK = 0;
    $isOK = 1 if(-e $dir and -d $dir);
    return $isOK;
}

# д������־
sub write_log{
    my $log_file = shift @_;
    my $title    = shift @_;
    my $hashCondition = shift @_;

    my $tab="##########";
    open LOG,">>$log_file";
    print LOG "$tab Start $title ".get_time()." $tab\n" if($title=~/\w/);
    # Finish
    my @finish_samples = exists $hashCondition->{'Finish'} ? keys %{$hashCondition->{'Finish'}} : ();
    print LOG 'Finish   : '. (join ",", @finish_samples) . "\n"; 
    # Good2Run
    my @good2run_samples = exists $hashCondition->{'Good2Run'} ? keys %{$hashCondition->{'Good2Run'}} : (); 
    print LOG 'Good2Run : '. (join ",", @good2run_samples) . "\n"; 
    # Error
    my @error_samples = exists $hashCondition->{'Error'} ? keys %{$hashCondition->{'Error'}} : ();
    print LOG 'Error    : '. (join ",", @error_samples) . "\n"; 

    print LOG "\n";
    close LOG;
}

# �������� ����������������
sub process_bar_array{
    my $process_name   = shift @_; # ��Ҫ�����Ķ���
    my $process_arrays = shift @_; # ���б�
    my $process_all_count = @$process_arrays;# ��Ҫ���������б�����
    my ($pos) = grep $$process_arrays[$_] eq $process_name, 0..$process_all_count-1;
    my $process_count = $pos+ 1; # ��ǰ�����λ��
    my $process_perc = sprintf "%0.2f", 100 * $process_count/$process_all_count; # ����
    my $windows = 100; # ���ڿ��
    my $finished = $windows * int($process_perc) / 100; # ���
    my $unfinished = $windows - $finished; # δ���
    print "\r[". '>' x $finished . ' ' x $unfinished . "]  [$process_count/$process_all_count]  $process_perc% ";
    print "\n" if($process_count == $process_all_count);   
}

1