package package::evidence;
use strict;
use warnings;
use Encode;
use Excel::Writer::XLSX;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Get Bam Evidence ".package::utils::get_time()." ##########\n";
    my $report_dir  = $hashConfig->{'Report'}; 
    my $output_dir  = $hashConfig->{'Output'}; 
    my $SVdbtype    = (exists $hashConfig->{'GenomeBed'}) ? "SVdbWGS" : "SVdbWES"; #判断是WES 还是WGS
    my $dataDir     = "$report_dir/sv/circosData";
    my @samples     = package::utils::get_sample($hashConfig, 'case', 'control');
    my %posInfo     = get_pos($dataDir,\@samples);

    # 数据状态检测
    my $excel    = "$report_dir/document/sv_evidence.xlsx";
    if(package::utils::is_file_ok($excel))
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }  

    get_bam_evidence($hashPara, $hashConfig, $SVdbtype, \%posInfo);
}

##获取样本CNV和SV位置信息
sub get_pos{
    my $dataDir = shift @_;
    my $samples = shift @_;

    my %posInfo = ();
    foreach my $sample(@$samples)
    {
        my $sample_loss = "$dataDir/$sample.loss";##loss gain 储存cnv信息
        my $sample_gain = "$dataDir/$sample.gain";
        my $sample_sv   = "$dataDir/$sample.sv";

        if(package::utils::is_file_ok($sample_loss) or package::utils::is_file_ok($sample_gain))
        {
            open CNVFile,"cat $sample_loss $sample_gain |";
            while(my $line = <CNVFile>)
            {
                $line =~ s/[\r\n]//g;
                my ($chr,$pos1,$pos2,$tmp) = split /\t/,$line;
                $chr =~ s/hs//g;
                $posInfo{"CNV"}{"$chr:$pos1"}{"$sample"}++;
                $posInfo{"CNV"}{"$chr:$pos2"}{"$sample"}++;
            }
            close CNVFile;
        }

        next if(package::utils::is_file_ok($sample_sv) == 0);
        open SVFile,$sample_sv;
        while(my $line = <SVFile>)
        {
            $line =~ s/[\r\n]//g;
            my ($chr1,$pos1,$tmp1,$chr2,$pos2,$tmp2) = split /\t/,$line;
            $chr1 =~ s/hs//g;
            $chr2 =~ s/hs//g;
            $posInfo{"SV"}{"$chr1:$pos1"}{"$sample"}++;
            $posInfo{"SV"}{"$chr2:$pos2"}{"$sample"}++;
        }
        close SVFile;
    }
    return %posInfo;
}

##输出bam文件在相应位置的信息
sub get_bam_evidence{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $SVdbtype   = shift @_;
    my $posInfo    = shift @_;

    my $SamTools            = $hashPara   -> {'Soft'}{'SamTools'};
    my $report_dir          = $hashConfig -> {'Report'};
    my $output_dir          = $hashConfig -> {'Output'};
    ##输出每个样本所有样本汇总位点信息(All)or输出每个样本自身位点信息(Sample)
    my $sample_evidence_out = (exists $hashConfig->{"sample_evidence_out"} and $hashConfig->{'sample_evidence_out'}=~/^All$|^Sample$/)? $hashConfig->{"sample_evidence_out"} : "Sample";
    my $bam_evidence_dir    = "$report_dir/sv/evidence";
    package::utils::make_dir($bam_evidence_dir);
    my @samples = package::utils::get_sample($hashConfig, 'case', 'control');
    my $excel   = "$report_dir/document/sv_evidence.xlsx";
    
    print "Output $excel ... ";
    my $workbook = Excel::Writer::XLSX->new($excel); 
    my %format   = package::utils::formatRun($workbook);
    my $worksheet_cnv  = $workbook->add_worksheet("CNV");    
    my $worksheet_sv   = $workbook->add_worksheet("SV"); 
    my @heads = ("Sample","Chr","Pos","QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","ISIZE/TLEN","SEQ","QUAL");
    for my $i(0..$#heads)
    {
        $worksheet_cnv -> write(0,$i,$heads[$i],$format{"title"});
        $worksheet_sv  -> write(0,$i,$heads[$i],$format{"title"});
    }
    my $row_cnv = 1;
    my $row_sv  = 1;
    foreach my $sample(@samples)
    {
        my $RawBam = ();
        if ($SVdbtype eq "SVdbWGS")
        {
            $RawBam = "$output_dir/$sample/$sample\_final.bam"; 
        }
        if ($SVdbtype eq "SVdbWES")
        {
            $RawBam = "$output_dir/$sample/svbam/$sample\_sv.bam";
        }
     
        ##write cnv bam evidence
        foreach my $title(sort keys $posInfo->{"CNV"})
        {
            next if($sample_evidence_out eq "Sample" and not exists $posInfo->{"CNV"}{"$title"}{"$sample"});
            my ($chr,$pos) = split /:/,$title;
            my $start      = $pos-50;
            my $end        = $pos+50;
            my $position   = "$chr:$start\-$end";
            my $evidence_file = "$bam_evidence_dir/$sample\_cnv\_$chr\_$pos";
            system("$SamTools view $RawBam $position > $evidence_file");
            next if(package::utils::is_file_ok($evidence_file) == 0);
            open File,$evidence_file;
            while(my $line = <File>)
            {
                $line =~ s/[\r\n]//g;
                my @datas = split /\t/,$line;
                @datas = ($sample,$chr,$pos,@datas);
                foreach my $i(0..$#datas)
                {
                    $worksheet_cnv -> write($row_cnv,$i,$datas[$i],$format{"normal"});
                }
                $row_cnv++;
             }
            close File;
        }

        ##write sv bam evidence
        foreach my $title(sort keys $posInfo->{"SV"})
        {
            next if($sample_evidence_out eq "Sample" and not exists $posInfo->{"SV"}{"$title"}{"$sample"});
            my ($chr,$pos) = split /:/,$title;
            my $start      = $pos-50;
            my $end        = $pos+50;
            my $position   = "$chr:$start\-$end";
            my $evidence_file = "$bam_evidence_dir/$sample\_sv\_$chr\_$pos";
            system("$SamTools view $RawBam $position > $evidence_file");
            next if(package::utils::is_file_ok($evidence_file) == 0);
            open File,$evidence_file;
            while(my $line = <File>)
            {
                $line =~ s/[\r\n]//g;
                my @datas = split /\t/,$line;
                @datas = ($sample,$chr,$pos,@datas);
                foreach my $i(0..$#datas)
                {
                    $worksheet_sv -> write($row_sv,$i,$datas[$i],$format{"normal"});
                }
                $row_sv++;
            }
            close File;
        }
    }

    write_readme($workbook,\%format,"readme_evidence.txt");
    print "OK\n";
}

sub write_readme{
    my $workbook = shift @_;
    my $format   = shift @_;
    my $file     = shift @_; 
    
    open IN,"$file";
    my $readme = $workbook -> add_worksheet("ReadMe");
    my $row = 0;
    $readme -> set_column(0,0,15);
    $readme -> set_column(1,1,50);
    while(<IN>){
        $_ =~ s/[\r\n]//g;
        my @datas = split /\t/,$_;
        $readme -> write($row,0,$datas[0], $format->{"title"});
        $readme -> write($row,1,decode("utf-8", $datas[1]), $format->{"normal"});
        $row++;     
    }
    close IN;
}

1