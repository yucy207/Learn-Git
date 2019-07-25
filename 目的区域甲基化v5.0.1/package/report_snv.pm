package package::report_snv;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start Report SNV Result ".package::utils::get_time()." ##########\n";

    # 路径准备
    my $report_dir   = $hashConfig->{'Report'};
    my $document_dir = "$report_dir/document";
    my $snv_dir      = "$report_dir/document/07_SNP_Calling";
    package::utils::make_dir($document_dir, $snv_dir);

    my $vcf_file = "$report_dir/vcf/allSample.vcf.split.gz";
    my %hashVCF  = read_VCF($hashPara, $hashConfig, $vcf_file);        # 读取所有突变位点
    my %hashAnno = annovar_convert($hashPara, $hashConfig, \%hashVCF); # 对突变位点格式转换为annovar格式（主要是插入缺失的问题）

    add_genome_position($hashConfig, \%hashAnno);   # 根据基因组区域注释信息，注释基因组位置与方向
    add_ref_seq_base($hashConfig, \%hashAnno);      # 添加参考序列上的ref碱基
    filter(\%hashAnno);                             # 根据信息对SNV进行假阳性过滤，主要是C->T->C类型突变
    annotation($hashPara, $hashConfig, \%hashAnno); # 注释1000g和SNPID信息
    add_statistic(\%hashVCF, \%hashAnno);           # 进行统计分析，HWE、突变频率、召回率
    output($hashConfig, \%hashVCF, \%hashAnno);     # 结果输出
}

sub output{
    my $hashConfig = shift @_;
    my $hashVCF    = shift @_;
    my $hashAnno   = shift @_;
    my $output_dir = "$hashConfig->{'Report'}/document/07_SNP_Calling";
    my $snv_excel  = "$output_dir/SNV_bwa.xlsx";
    print "Output $snv_excel\n";
    my $workbook   = Excel::Writer::XLSX->new($snv_excel);
    my %format     = package::utils::sheet_format($workbook);
    my @annoNames  = ('Target','Position','Original Ref','Ref','Alt','Warning','SNP Filter','GenomeChr','GenomePosition','GenomeStrand','SNPID','1000g_Freq','Geno(0|1|2)','Alt Allele Freq','CallRate','HWE');
    my @allSamples = split /,/,$hashVCF->{'allSampleName'};
    my @titles     = (@annoNames, @allSamples);
    my $sheetGeno  = $workbook->add_worksheet("GenoType");
    my $sheetSNV   = $workbook->add_worksheet("SNV");
    my $row=0;
    for my $col(0..$#titles){
        $sheetGeno->write($row,$col,$titles[$col],$format{"title"});
        $sheetSNV->write($row,$col,$titles[$col],$format{"title"});
    }
    $row++;
    foreach my $SNVtitle(sort keys %$hashVCF){
        next if($SNVtitle eq 'allSampleName');
        next if($hashAnno->{$SNVtitle}{'Filter'} eq 'REJECT');
        my $ref = $hashAnno->{$SNVtitle}{"Ref"};
        my $alt = $hashAnno->{$SNVtitle}{"Alt"};
        my @datasGeno;
        my @datasSNV;
        my @colors;
        # 注释部分信息提取
        foreach my $annoName(@annoNames){
            my $value = (exists($hashAnno->{$SNVtitle}{$annoName})) ? $hashAnno->{$SNVtitle}{$annoName} : "";            
            push @datasGeno,$value;
            push @datasSNV,$value;
            push @colors,"normal";
        }
        # 分型部分信息提取
        foreach my $sample(@allSamples){
            my $genotype = (exists($hashVCF->{$SNVtitle}{$sample}{'GenoType'})) ? $hashVCF->{$SNVtitle}{$sample}{'GenoType'} : "";
            my $refDepth = (exists($hashVCF->{$SNVtitle}{$sample}{'RefDepth'})) ? $hashVCF->{$SNVtitle}{$sample}{'RefDepth'} : "";
            my $altDepth = (exists($hashVCF->{$SNVtitle}{$sample}{'AltDepth'})) ? $hashVCF->{$SNVtitle}{$sample}{'AltDepth'} : "";
            my ($geno,$snv) = (" "," ");
            if($genotype=~/\w/){
                $geno = "$ref/$ref"  if($genotype eq 'HOMR');
                $geno = "$ref/$alt"  if($genotype eq 'HET');
                $geno = "$alt/$alt"  if($genotype eq 'HOMA');
                my $altFreq = sprintf "%0.4f",$altDepth/($refDepth+$altDepth); 
                $snv  = "$genotype:$ref=$refDepth,$alt=$altDepth:$altFreq";        	
            }
            push @datasGeno,$geno;
            push @datasSNV,$snv;
            push @colors,"normal";
        }
        # excel输出
        for my $col(0..@datasGeno-1){
        	$sheetGeno->write($row,$col,$datasGeno[$col],$format{$colors[$col]});
        	$sheetSNV->write($row,$col,$datasSNV[$col],$format{$colors[$col]});
        }			
		$row++;
    }
    readme($workbook, \%format, "readme/readme_bwa.txt");
}

sub add_statistic{
	my $hashVCF  = shift @_;
    my $hashAnno = shift @_;
    print "Statistic HWE/alt allele freq/call rate\n";
    my @allSamples   = split /,/,$hashVCF->{'allSampleName'};
    my $allSampleNum = @allSamples;
    foreach my $title(sort keys %$hashVCF){
        next if($title eq 'allSampleName');
        my ($homrC, $hetC, $homaC) = (0,0,0);# 三种基因型数量
        foreach my $sample(@allSamples){
        	next if(not exists($hashVCF->{$title}{$sample}));
            my $genotype = $hashVCF->{$title}{$sample}{'GenoType'};
            $homrC++ if($genotype eq 'HOMR');
            $hetC++  if($genotype eq 'HET');
            $homaC++ if($genotype eq 'HOMA');
        }
        my $hwe      = package::hwe::run($hetC,$homrC,$homaC);     # hwe
        my $freq     = ($hetC+2*$homaC)/(2*($homrC+$hetC+$homaC)); # 突变频率
        my $callRate = ($homrC+$hetC+$homaC)/$allSampleNum;        # 召回率
        $hashAnno->{$title}{"HWE"}              = $hwe;
        $hashAnno->{$title}{"Alt Allele Freq"}  = $freq;
        $hashAnno->{$title}{"CallRate"}         = $callRate;
        $hashAnno->{$title}{"Geno(0|1|2)"}      = "$homrC|$hetC|$homaC";
    }
}

sub annotation{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $hashAnno   = shift @_;
    my $species    = $hashConfig->{"Species"};
    return if($species!~/Human/);
    my $library_file = "$hashConfig->{'Report'}/vcf/library";
    # 生成library文件
    open LIB,">$library_file";
    foreach my $title(sort keys %$hashAnno){
        my $chr    = $hashAnno->{$title}{'GenomeChr'};
        my $pos    = $hashAnno->{$title}{'GenomePosition'};
        my $ref    = $hashAnno->{$title}{'Original Ref'};
           $ref    = "-" if($ref !~ /\w/);
        my $alt    = $hashAnno->{$title}{'Alt'};
        my $strand = $hashAnno->{$title}{'GenomeStrand'};
        my @tmps   = split/\-/,$pos;
        my($start, $end) = @tmps[0,$#tmps];
        if($strand eq "-"){ # 负链需进行碱基互补后再注释
            $ref =~ tr/ACGTacgt/TGCAtgca/ if($ref =~ /\w/);            
            $alt =~ tr/ACGTacgt/TGCAtgca/;
        }
        print LIB "$chr\t$start\t$end\t$ref\t$alt\t$title\n";
    }
    close LIB;
    # Annovar注释   
    my $annovarDir   = $hashPara->{'Soft'}{'AnnovarDIR'};
    my $AnnovarBuild = $hashPara->{$species}{'AnnovarBuild'};
    system("$annovarDir/annotate_variation.pl -filter -dbtype avsnp150 --buildver $AnnovarBuild $library_file /home/genesky/database/annovar/$AnnovarBuild/avsnp");
    system("$annovarDir/annotate_variation.pl -filter -dbtype 1000g2015aug_all --buildver $AnnovarBuild $library_file /home/genesky/database/annovar/$AnnovarBuild/1000g");
    
    # 注释结果读取
    my $dbsnp_result = "$library_file.$AnnovarBuild\_avsnp150_dropped";
    my $g1000_result = "$library_file.$AnnovarBuild\_ALL.sites.2015_08_dropped";
    readAnno($dbsnp_result, "SNPID", $hashAnno);
    readAnno($g1000_result, "1000g_Freq", $hashAnno);
}

sub readAnno{
    my $file     = shift @_;
    my $key      = shift @_;
    my $hashAnno = shift @_;
    open IN,$file;
    while(<IN>){
        $_=~s/[\r\n]//g;
        my @datas = split/\t/,$_;
        my $anno  = $datas[1];
        my $title = $datas[$#datas];
        $hashAnno->{$title}{$key} = $anno;
    }
    close IN;
}

sub filter{
    my $hashAnno = shift @_;
    foreach my $title(sort keys %$hashAnno){
        my $originalRef = $hashAnno->{$title}{'Original Ref'};
        my $ref         = $hashAnno->{$title}{'Ref'};
        my $alt         = $hashAnno->{$title}{'Alt'};
        my $mark        = "PASS";
        my $warning     = " ";
        $mark    = 'REJECT'                     if($originalRef eq 'C' and $alt eq 'C');# 参考是C
        $warning = 'Alt allele C is methylated' if($originalRef eq 'A' and $alt eq 'C');# 参考是A
        $warning = 'Alt allele could be T or C' if($originalRef eq 'A' and $alt eq 'T');# 参考是A
        $warning = 'Alt allele C is methylated' if($originalRef eq 'G' and $alt eq 'C');# 参考是G
        $warning = 'Alt allele could be T or C' if($originalRef eq 'G' and $alt eq 'T');# 参考是G
        $warning = 'Alt allele C is methylated' if($originalRef eq 'T' and $alt eq 'C');# 参考是T
 
        $hashAnno->{$title}{'Filter'}  = $mark;
        $hashAnno->{$title}{'Warning'} = $warning;
    }
}

sub add_ref_seq_base{
	my $hashConfig = shift @_;
	my $hashAnno   = shift @_;
    my %hashSeq    = package::utils::read_fasta($hashConfig->{'TargetFasta'});
	print "Add referrence genome base\n";
	foreach my $title(keys %$hashAnno){
        my $target = $hashAnno->{$title}{'Target'};# 片段名称
        my $start  = $hashAnno->{$title}{'Start'}; # 突变在片段上起点
        my $end    = $hashAnno->{$title}{'End'};   # 突变在片段上终点
        my $ref    = $hashAnno->{$title}{'Ref'};   # 突变甲基化后的ref
        my $originalRef = ($ref eq '-') ? '' : substr($hashSeq{$target},$start-1,$end-$start+1);
        $hashAnno->{$title}{'Original Ref'} = uc($originalRef);
	}
}

sub add_genome_position{
	my $hashConfig    = shift @_;
	my $hashAnno      = shift @_;
	my %hashSeqRegion = package::utils::read_seq_region("$hashConfig->{'Report'}/seq.region"); # 片段在参考基因组上的位置信息 
	print "Add referrence genome position\n";
    foreach my $title(keys %$hashAnno){
        my $target   = $hashAnno->{$title}{'Target'};# 片段名称
        my $start    = $hashAnno->{$title}{'Start'}; # 突变在片段上起点
        my $end      = $hashAnno->{$title}{'End'};   # 突变在片段上终点
        my $position = $hashAnno->{$title}{'Position'};   # 突变在片段上位置
        my $targetChr    = $hashSeqRegion{$target}{'Chr'};  # 片段的染色体
        my $targetStart  = $hashSeqRegion{$target}{'Start'};# 片段在参考基因组上的起点
        my $targetEnd    = $hashSeqRegion{$target}{'End'};  # 片段在参考基因组上的终点
        my $targetStrand = $hashSeqRegion{$target}{'Target_Strand'};# 片段在参考基因组上的方向
        my $genomeStart = "";# 突变在参考基因组上的起点
        my $genomeEnd   = "";# 突变在参考基因组上的终点
        if($targetStart=~/\d/){# 有参考基因组时，对基因组位置重新计算
            $genomeStart = ($targetStrand eq '+') ? $targetStart+$start-1 : $targetStart-$start+1;
            $genomeEnd   = ($targetStrand eq '+') ? $targetStart+$end-1 : $targetStart-$end+1;
        }
        $hashAnno->{$title}{'GenomeChr'}      = $targetChr;   # 突变在参考基因组上的染色体
        $hashAnno->{$title}{'GenomeStart'}    = $genomeStart; # 突变在参考基因组上的起始
        $hashAnno->{$title}{'GenomeEnd'}      = $genomeEnd;   # 突变在参考基因组上的终点
        $hashAnno->{$title}{'GenomePosition'} = ($position=~/-/) ? "$genomeStart-$genomeEnd" : $genomeStart;   # 突变在参考基因组上的位置
        $hashAnno->{$title}{'GenomeStrand'}   = $targetStrand;# 突变在参考基因组上的方向
    }
}

sub annovar_convert{
    my $hashPara    = shift @_;
    my $hashConfig  = shift @_;
    my $hashVCF     = shift @_;
    my $report_dir  = $hashConfig->{'Report'};
    my $annovar_dir = $hashPara->{'Soft'}{'AnnovarDIR'};
    my $vcf_dir     = "$report_dir/vcf";
    my $vcf_input   = "$vcf_dir/convert.input";  # 存放需要转换的信息
    my $vcf_output  = "$vcf_dir/convert.output"; # 转换后的信息
    print "Convert VCF -> annovar Format\n";
    # 输出需要转换的位点
    open CONVERT, ">$vcf_input";
    print CONVERT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    foreach my $title(sort keys %$hashVCF){
    	next if($title eq 'allSampleName');
        my ($target, $pos, $ref, $alt) = split /\|/, $title;
        print CONVERT "$target\t$pos\t.\t$ref\t$alt\t.\t.\t$title\n";
    }
    close CONVERT;
    system("$annovar_dir/convert2annovar.pl -format vcf4 $vcf_input -outfile $vcf_output -includeinfo");
    # 读取注释后的结果
    my %hashAnno;
    open CONVERTOUTPUT, $vcf_output;
    while(<CONVERTOUTPUT>){
        $_=~s/[\r\n]//g;
        my ($target, $start, $end, $ref, $alt, @tmps) = split /\t/,$_;
        $end = $start + length($ref) - 1;
        $ref = "-" if($ref=~/\d/ and $ref==0); # vcf文件如果alt是-的话，转换后就是0（这个-是vcf的*字符替换的）
        $alt = "-" if($alt=~/\d/ and $alt==0);
        my $title = $tmps[$#tmps];             # VCF哈希中的title
        $hashAnno{$title}{"Target"}     = $target;
        $hashAnno{$title}{"Start"}      = $start;
        $hashAnno{$title}{"End"}        = $end;
        $hashAnno{$title}{"Ref"}        = $ref;
        $hashAnno{$title}{"Alt"}        = $alt;
        $hashAnno{$title}{"Position"}   = ($ref eq '-' or $alt eq '-') ? "$start-$end" : $start;
        $hashAnno{$title}{"SNP Filter"} = (exists $hashVCF->{$title}{'SNP Filter'})? $hashVCF->{$title}{'SNP Filter'} : "";
    }
    close CONVERTOUTPUT;
    return %hashAnno;
}

sub read_VCF{
	# 读取VCF文件
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $vcf_file   = shift @_;
    my %hashVCF;
    my @heads;   # 标题
    my @samples; # 样本
    print "Reading $vcf_file ... ";
    open VCF, "gzip -dc $vcf_file|";
    while(<VCF>){
        $_=~s/[\r\n]//g;
        next if($_=~/^##/);
        if($_=~/^#CHROM/){
            $_ =~s/^#//;
            @heads = split/\t/, $_;
            @samples = @heads[9..$#heads];
            next;
        }
        # 把数据读到临时哈希中
        my @datas=split /\t/,$_;
        my %hashTmp;
        foreach my $col(0..$#heads){
            $hashTmp{$heads[$col]}=$datas[$col];
        }
        my $target     = $hashTmp{'CHROM'};  # 染色体
        my $pos        = $hashTmp{'POS'};    # 位置
        my $filter     = $hashTmp{'FILTER'}; # 硬过滤注释
        my $refAllele  = $hashTmp{'REF'};    # 文件中参考碱基
        my $altAllele  = $hashTmp{'ALT'};    # 文件中突变碱基
        next if($altAllele eq '*');          # *表示缺失
        my $title = "$target|$pos|$refAllele|$altAllele";
        my $flag  = 0;
        foreach my $sample(@samples){
            my ($gt,$ad,$dp,$tmp) = split /:/, $hashTmp{$sample};# 样本的详细分型信息
            my @depths=split /,/,$ad;                            # 每一个等位基因的深度
            my $refDepth = $depths[0];                           # 参考碱基深度
            my $altDepth = $depths[1];                           # 突变碱基深度
            my $freq     = ( ($refDepth + $altDepth) > 0 ) ? $altDepth/($refDepth + $altDepth) : 0;
            next if( ($refDepth + $altDepth) <= 20 ); #测序深度<=20,不要
            next if($gt eq './.');                    # 分型失败,不要
            my $genoType = "";                        # 对当前位点定义分型
               $genoType = "HET"  if($gt eq "0/1" or $gt eq "1/0");
               $genoType = "HOMR" if($gt eq "0/0");
               $genoType = "HOMA" if($gt eq "1/1");
            next if($genoType ne "HOMR" and $freq < 0.2);
            $flag = 1 if($genoType eq "HET" or $genoType eq "HOMA");
            $hashVCF{$title}{$sample}{"AltDepth"}  = $altDepth;
            $hashVCF{$title}{$sample}{"RefDepth"}  = $refDepth;
            $hashVCF{$title}{$sample}{"GenoType"}  = $genoType;
        }
        $hashVCF{$title}{'SNP Filter'} = $filter if(exists $hashVCF{$title});
        delete $hashVCF{$title} if(exists $hashVCF{$title} and $flag==0);
    }
    $hashVCF{'allSampleName'}=join ",",@samples;
    close VCF;
    print "OK\n";
    return %hashVCF;
}

sub readme{
	my ($workbook,$format,$file)=@_;
    # package路径
    my $package_dir = Cwd::abs_path(package::utils::get_dirname(File::Spec->rel2abs(__FILE__)));
    return if(package::utils::is_file_ok("$package_dir/$file")==0);
    

	my $sheet=$workbook->add_worksheet("Read Me");
	$sheet->set_row(0, 65);
	$sheet->set_column('A:A', 35);
	$sheet->set_column('B:B', 115);
	my $row=0;
	open FILE,"$package_dir/$file";
	while(<FILE>){
		$_=~ s/\^/\n/g;
		my @split_line=split /\t/,$_;
		my $col=0;
		foreach(@split_line)
		{
			my $text = Encode::decode("UTF-8",$_);
			if($row==0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme1'});}
			if($row==0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme2'});}
			if($row>0 and $col==0){$sheet->write($row,$col,$text,$format->{'readme3'});}
			if($row>0 and $col==1){$sheet->write($row,$col,$text,$format->{'readme5'});}
			$col++;
		}
		$row++;
	}
	close FILE;
}

1
