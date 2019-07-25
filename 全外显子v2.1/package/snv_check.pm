package package::snv_check;
use strict;
use warnings;
use Encode;
use Excel::Writer::XLSX;
use Statistics::R;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    print "########## Start SNV check ".package::utils::get_time()." ##########\n";
    my $report_dir   = $hashConfig->{'Report'}; 
    my $document_dir = "$report_dir/document";
    my $snv_dir      = "$document_dir/2_SNV";
    my $snv_excel    = "$snv_dir/SNVcount.xlsx";
    
    #####
    # 写日志
    #####
    package::utils::write_log_simple("$report_dir/run.log", "SNV check");
    if(package::utils::is_file_ok($snv_excel) == 1 and package::utils::is_force_sample($hashConfig) == 0)
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }

    # 文件、路径准备、检测
    my @samples      = package::utils::get_sample($hashConfig, 'case', 'control');
    my $analysis_txt = "$report_dir/analysis_all.txt"; # 突变基因型文件
    my $database_snv = "$report_dir/database.snv";     # 突变信息文件
    die "[Error] Lost $analysis_txt" if(package::utils::is_file_ok($analysis_txt) == 0);
    die "[Error] Lost $database_snv" if(package::utils::is_file_ok($database_snv) == 0);

    # snv统计及绘图
    snv_count_analysis($hashConfig, $hashPara, \@samples, $analysis_txt, $database_snv);
}

# snv统计及绘图
sub snv_count_analysis{
    my $hashConfig   = shift @_;
    my $hashPara     = shift @_;
    my $samples      = shift @_;
    my $analysis_txt = shift @_;
    my $database_snv = shift @_;
    
    my $report_dir   = $hashConfig->{'Report'}; 
    my $document_dir = "$report_dir/document";
    my $snv_dir      = "$document_dir/2_SNV";
    my $snv_excel    = "$snv_dir/SNVcount.xlsx";
    package::utils::make_dir($document_dir);
    package::utils::make_dir($snv_dir);

    # fastqc文件拷贝
    prepare_QC_image($hashConfig, $samples); 

    # fastqc_clean 文件拷贝
    prepare_QC_Clean_image($hashConfig, $samples); 

    # 文件读取和数据统计汇总
    my %hashSNV     = readLibrary($database_snv);
    my %hashGeno    = readAnalysis($analysis_txt);
    my %hashSummary = summary_snv(\%hashSNV, \%hashGeno, $samples);
    
    # excel表格输出
    my $workbook = Excel::Writer::XLSX->new($snv_excel);
    my %format   = package::format::run($workbook);
    output_gene_region_count($workbook, \%format, \%hashSummary, $samples); # 基因区域sheet
    output_function_count($workbook, \%format, \%hashSummary, $samples); # 突变功能sheet
    output_gatk($workbook, \%format, \%hashSummary, $samples); # GATK 突变snv/indel sheet
    output_hom_het_novel($workbook, \%format, \%hashSummary, $samples); # HomHetNovel sheet
    output_ts_vs_tv($workbook, \%format, \%hashSummary, $samples); # Ts Tv sheet
    output_indel_length($workbook, \%format, \%hashSummary, $samples); # InDel Length sheet
    
    # 绘图结果输出
    my $Rbin = $hashPara->{'Soft'}{'R'};
    my $Rlib = $hashPara->{'Soft'}{'RLib'};
    my $R = Statistics::R->new(bin => $Rbin);
    $R->startR;
    $R->send(qq` .libPaths("$Rlib") \n`) ;
    $R->send(qq` library(ggplot2)\n`);
    $R->send(qq` library(scales)\n`);
    plot_pie(\%hashSummary, $snv_dir, 'Region_Summary', 'Region_pie', $R);
    plot_pie(\%hashSummary, $snv_dir, 'Function_Summary', 'Function_pie', $R);
    plot_pie(\%hashSummary, $snv_dir, 'TsTv_Summary', 'TsTv_pie', $R);
    plot_bar(\%hashSummary, $snv_dir, 'Indel_Length_Summary_Plot', 'Indel_distribution', $R);
    plot_violin_plot(\%hashSummary, $snv_dir, 'Alt_Freq_Summary', 'snv_frequency_violin_plot', $R, $samples);
    $R->stop();

    # 亲缘关系检测
    cluster_analysis($hashConfig, $hashPara, $samples, \%hashSNV, \%hashGeno) if( @$samples >= 2);
}


# 亲缘关系检测
sub cluster_analysis{
    my $hashConfig   = shift @_;
    my $hashPara     = shift @_;
    my $samples      = shift @_;
    my $hashSNV      = shift @_;
    my $hashGeno     = shift @_;
    
    # 软件
    my $king     = $hashPara->{'Soft'}{'king'};
    my $plink    = $hashPara->{'Soft'}{'plink'};
    my $fasttree = $hashPara->{'Soft'}{'fasttree'};
    
    # 数据准备和生成样本的fasta序列
    my %hashseq  = getSeq($hashConfig, $hashSNV, $hashGeno, $samples);
    generated_fasta($hashConfig, $samples, \%hashseq);
    
    # 数据分析：tree,king
    data_analysis($hashConfig, $fasttree, $plink, $king); 
    
    # 输出cluster.xlsx及绘图
    plot_cluster($hashConfig, $hashPara, $samples);
    report_cluster($hashConfig);
    
    # 输出不配对样本（主要针对肿瘤流程）
    unpaired_print($hashConfig);
}

# 输出cluster.xlsx
sub report_cluster{
    my $hashConfig   = shift @_;
    my $report_dir   = $hashConfig->{"Report"};
    my $relationFile ="$report_dir/FastTree/relationShip.kin0";
    my $excel    = "$report_dir/document/1_QC/cluster.xlsx";
    my $workbook = Excel::Writer::XLSX->new($excel);
    my %format   = package::format::run($workbook);
    my $sheet    = $workbook->add_worksheet("cluster");
    my %hash;
    open IN,$relationFile;
    my $title = <IN>;
       $title  =~ s/[\r\n]//;chomp $title;
    my @titles = split/\t/,$title;
    while(<IN>){
        $_ =~ s/[\r\n]//;chomp;
        my @datas = split/\t/,$_;
        my($id1,$id2) = @datas[0,2];
        foreach my $i(0..$#titles){
            $hash{"$id1|$id2"}{$titles[$i]} = (exists $datas[$i])?$datas[$i]:"";
        }
    }
    close IN;
    my $row=0;
    foreach my $col(0..$#titles){
        $sheet->write($row,$col,$titles[$col],$format{'title'});
    }
    $row++;
    foreach my $key(sort keys %hash){
        foreach my $col(0..$#titles){
            my $value = " ";
               $value = $hash{$key}{$titles[$col]} if(exists $hash{$key}{$titles[$col]});
            $sheet->write($row,$col,$value,$format{'normal'});
        }
        $row++;
    }
    ReadMe($workbook,\%format,"readme_cluster.txt");

    #亲缘关系为同一个体/双胞胎样本输出
    my $string = "Warnings : These samples' kinship > 0.354 : ";
    my $flag = 0;
    foreach my $key(sort keys %hash){
        my $value = (exists $hash{$key}{"Kinship"}) ? $hash{$key}{"Kinship"} : 0;
        if($value > 0.354){
            my($sample1, $sample2) = split/\|/,$key;
            $string.="$sample1,$sample2; ";
            $flag = 1;
        }
    }
    system "echo -e \"\\033[41;37m $string \\033[0m\"" if($flag == 1); # 警告 红底白字
}

sub ReadMe{
    my ($workbook,$format,$file) = @_;
    return if(package::utils::is_file_ok($file) == 0);
    my $sheet=$workbook->add_worksheet("Read Me");
    $sheet->set_row(0, 65);
    $sheet->set_column('A:A', 35);
    $sheet->set_column('B:B', 60);
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

# 输出不配对样本（主要针对肿瘤流程）
sub unpaired_print{
    my $hashConfig = shift @_;
    if(exists $hashConfig->{"SomaticPair"}){
        print "Check SomaticPair\n";
        my $somaticPair = $hashConfig->{"SomaticPair"};
        my @pairs = split /\t/,$somaticPair;
        my %hashpair = ();
        for(my $i=0;$i<=$#pairs;$i++){
            my @line = split /,/,$pairs[$i];
            $hashpair{$line[0]}{$line[1]}++;
        }
        my $report_dir   = $hashConfig->{"Report"};
        my $relationFile = "$report_dir/FastTree/relationShip.kin0";
        open IN,$relationFile;
        while(<IN>){
            $_ =~ s/[\r\n]//g;
            next if($_ =~ /^FID1/);
            my @data = split /\s+/,$_;
            if(exists($hashpair{$data[1]}{$data[3]}) or exists($hashpair{$data[3]}{$data[1]}))
            {
                next if($data[6] > 0.354);
                print "[Warning] $data[1]\t$data[3]\t$data[6] < 0.354\n";
            }
        }    
        close IN;
        print "Check SomaticPair Finished\n";
    }
}

# 绘图cluster.pdf
sub plot_cluster{
    my $hashConfig = shift @_;
    my $hashPara   = shift @_;
    my $samples    = shift @_;
    my $sampleNum  = @$samples;
    my $report_dir = $hashConfig->{"Report"};
    my $document   = "$report_dir/document";package::utils::make_dir($document);
    my $qc_dir     = "$document/1_QC";package::utils::make_dir($qc_dir);
    my $treefile   = "$report_dir/FastTree/FastTree.fasta.tre";
    my $relationFile    = "$report_dir/FastTree/relationShip.kin0";
    my $relationFileNew = readRelation($relationFile,$samples);
    my $pdf      = "$qc_dir/cluster.pdf";
    my $pdfWidth = ($sampleNum>10) ? $sampleNum : 10;
    print "Processing $pdf ...";
    my $Rbin = $hashPara->{'Soft'}{'R'};
    my $Rlib = $hashPara->{'Soft'}{'RLib'};
    my $R = Statistics::R->new(bin => $Rbin);
    $R->startR;
    $R->send(qq` .libPaths("$Rlib") \n`) ;
    $R->send(qq` library("ggplot2") \n`);
    $R->send(qq` library("ggtree") \n`);
    $R->send(qq` library("grid") \n`);
    $R->send(qq` sampleNum = $sampleNum \n`);
    $R->send(qq` tree = read.tree("$treefile") \n`);
    $R->send(qq` relation=read.table("$relationFileNew",sep="\t",header=T,row.name=1) \n`);
    $R->send(qq` colnames(relation)=rownames(relation) \n`);
    $R->send(qq` relation=as.matrix(relation) \n`);
    $R->send(qq` pdf("$pdf", width =$pdfWidth, height =$pdfWidth) \n`);
    $R->send(qq` ggtree(tree, layout="circular") + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = 2+sampleNum/10) + geom_tiplab2(size=4+sampleNum/10) + ggtitle("Fan Dendrogram")\n`); ##扇形
    $R->send(qq` ggtree(tree) + geom_text2(aes(label=as.numeric(label)*100,subset=!isTip), hjust=-.3,size = 2+sampleNum/10) + geom_tiplab(aes(x=branch), size=4+sampleNum/10, vjust=-0.3) + ggtitle("Cluster Dendrogram")\n`);
    $R->send(qq` panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){usr <- par("usr"); on.exit(par(usr)); par(usr = c(0, 1, 0, 1)); z=x[!is.na(y)]; txt=as.numeric( sprintf( "\%0.4f", z[length(z)] ) ); if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt);color=1;if(txt>=0.354) color=2;if(txt>=0.177 && txt<0.354) color=3;if(txt>=0.0884 && txt<0.177) color=4;if(txt>=0.0442 && txt<0.0884) color=5; text(0.5, 0.5, txt, cex = cex.cor,col=color);} `);
    $R->send(qq` panel.text<-function(x, y, labels, cex, font, ...){cex.cor <- 0.8/strwidth(labels);text(0.5, 0.5, labels, cex = cex.cor);} `);
    $R->send(qq` pairs(relation ,lower.panel = NULL ,upper.panel = panel.cor,text.panel=panel.text ,font.labels = 2, main="Sample Relationship (Based On King software)") `);    
    $R->send(qq` info=c(">0.354                = duplicate/MZ twin\n\n","[0.177, 0.354]     = 1st-degree\n\n","[0.0884, 0.177]   = 2nd-degree\n\n","[0.0442, 0.0884] = 3rd-degree\n\n") `);    
    $R->send(qq` textExpand=sampleNum\n`);
    $R->send(qq` if(sampleNum<10) textExpand=10 \n`);
    $R->send(qq` mtext(info,side=1,adj=0,cex=2*textExpand/10,line=c(-2*textExpand/10,0,2*textExpand/10,4*textExpand/10),col=c(2,3,4,5)) `);     
    $R->send(qq` dev.off() `);
    $R->stop();    
    print "OK\n";    
}

# 把样本关系文件转格式，使R能够有效读取
sub readRelation{
    my $file    = shift @_;
    my $samples = shift @_;
    my (%hashRelation,@heads,%hashsamples);
    my $row=0;
    open IN,$file;
    while(<IN>){
        $_ =~ s/[\r\n]//g;
        next if($_ !~ /\w/);
        my @datas = split /\s+/,$_;
        $row++;
        if($row == 1){
            @heads = @datas;
            next;
        }
        my %hashtmp;
        foreach my $col(0..$#heads){
            my $value = (exists($datas[$col])) ? $datas[$col] : "";
            $hashtmp{$heads[$col]} = $value;
        }
        $hashRelation{$hashtmp{"ID1"}}{$hashtmp{"ID2"}} = $hashtmp{"Kinship"};
        $hashRelation{$hashtmp{"ID2"}}{$hashtmp{"ID1"}} = $hashtmp{"Kinship"};
        $hashsamples{$hashtmp{"ID1"}}++;
        $hashsamples{$hashtmp{"ID2"}}++;
    }
    close IN;

    my $outfile = "$file.reform";
    open OUT,">$outfile";
    print OUT "Sample";
    map{print OUT "\t$_";} @$samples; 
    print OUT "\n";
    $row=0;
    foreach my $sample_row(@$samples){
        print OUT "$sample_row";
        $row++;
        my $col=0;
        foreach my $sample_col(@$samples){
            $col++;
            my $value = (exists($hashRelation{$sample_row}{$sample_col})) ? $hashRelation{$sample_row}{$sample_col} : 0;
            $value = 1 if($row == $col);
            $value = "" if($row > $col);
            print OUT "\t$value";
        }
        print OUT "\n";
    }
    close OUT;
    return $outfile;
}

# 数据分析：tree,king
sub data_analysis{
    my $hashConfig   = shift @_;
    my $fasttree     = shift @_;
    my $plink        = shift @_;
    my $king         = shift @_;
    my $report_dir   = $hashConfig->{'Report'};
    my $fastafile = "$report_dir/FastTree/fasta";
    my $treefile  = "$report_dir/FastTree/FastTree.fasta.tre";
    system("$fasttree -nt $fastafile > $treefile");  
    system("$plink --file $report_dir/FastTree/sample --make-bed --out $report_dir/FastTree/bsample --noweb");# 二进制转化
    system("$king -b $report_dir/FastTree/bsample.bed  --homo --prefix $report_dir/FastTree/relationShip");# 关系计算 
}

# 生成样本的fasta序列
sub generated_fasta{
    my $hashConfig =shift @_;
    my $samples    = shift @_;
    my $hashseq    = shift @_;
    my $report_dir = $hashConfig->{'Report'};
    my $fasttree_dir = "$report_dir/FastTree";
    package::utils::make_dir($fasttree_dir);
    open IN,">$fasttree_dir/fasta";
    foreach my $sample(@$samples){
        print IN">$sample\n";      
        my $seq=$$hashseq{$sample};
        print IN"$seq\n";
    }
    close IN;
}

sub getSeq{
    my $hashConfig   = shift @_;
    my $library      = shift @_;
    my $analysis     = shift @_;
    my $samples      = shift @_;

    my $report_dir = $hashConfig->{'Report'};
    my (%hashseq, %analysis_filter) = ();
    # 生成序列
    foreach my $title(keys %{$library})
    {
        foreach my $alt(keys %{$library->{$title}})
        {
            my $blast = $library->{$title}{$alt}{"HOMOLOGY HITS"};
            my $freq  = $library->{$title}{$alt}{"Freq_Alt (1000g)"};
            my ($chr,$pos,$ref) = split /\|/,$title;
            ##过滤
            next if($blast > 1);
            next if($freq =~ /\d+/ and $freq>0.5);
            next if($title =~ /[XY]/);
            my $ref_length = length($ref);
            my $alt_length = length($alt);
            next if($ref_length > 200 or $alt_length > 200);

            my $s0 = $analysis->{$title}{$alt}{"Mutation 0"};
            my $s1 = $analysis->{$title}{$alt}{"Mutation 1"};
            my $s2 = $analysis->{$title}{$alt}{"Mutation 2"};
            my $samplesinfo = $s0.$s1.$s2;
            my @temp = split /\,/,$samplesinfo;
            my %hashTmp;
            foreach (@temp){
                next if($_ eq "");
                my ($sample,$type,$geno,$refdp,$altdp) = split /\:/,$_;
                $hashTmp{$sample}{'geno'} = $geno;
                $hashTmp{$sample}{'type'} = $type;
            }

            # 检查当前样本群体下，该位点是否发生突变
            my $isAlt = 0;
            foreach my $sample(@$samples){
                next if(!exists($hashTmp{$sample}) or $hashTmp{$sample}{'type'} eq 'HOMR');
                $isAlt = 1;
            }
            if($isAlt == 1){
                foreach my $sample(keys %hashTmp){
                    $analysis_filter{"$title|$alt"}{$sample} = $hashTmp{$sample}{'geno'};
                }
            }
        }
    }
    
    my @snvs;
    foreach my $title(keys %analysis_filter)
    {
        my ($chr,$pos,$ref,$alt) = split /\|/,$title;
        my ($Homr,$Het,$Homa,$Loss) = reGenotype($ref,$alt); 
        my $HomrGeno = "$ref/$ref";
        my $HetGeno  = "$ref/$alt";
        my $HomaGeno = "$alt/$alt";
        push @snvs,$title;
        foreach my $sample(@$samples){
            my $geno = exists ($analysis_filter{$title}{$sample}) ? $analysis_filter{$title}{$sample} : ""; 
            my $seq = $Loss;
            $seq = $Homr if($geno eq $HomrGeno);
            $seq = $Het if($geno eq $HetGeno);
            $seq = $Homa if($geno eq $HomaGeno);
            $hashseq{$sample} .= "$seq";
        }
    }

    # 生成ped map
    my $fasttree_dir = "$report_dir/FastTree";
    package::utils::make_dir($fasttree_dir);
    my $ped = "$fasttree_dir/sample.ped";
    my $map = "$fasttree_dir/sample.map";
    open MAP,">$map";
    foreach my $title(@snvs){
        my ($chr,$position,$ref,$alt) = split /\|/,$title;
        my $start = (split /-/,$position)[0];
        print MAP "$chr\t$title\t0\t$start\n";
    }
    close MAP;
    open PED,">$ped";
    foreach my $sample(@$samples){
        print PED "$sample\t$sample\t0\t0\t0\t0";
        foreach my $title(@snvs){
            my $geno = (exists $analysis_filter{$title}{$sample}) ? $analysis_filter{$title}{$sample} : "0/0";
            $geno =~ s/\// /;
            print PED "\t$geno";
        }
        print PED "\n";
    }
    close PED;
    
    return %hashseq;
}

sub reGenotype{
    my $ref  = shift @_;
    my $alt  = shift @_;
    my $Homr = "$ref$ref";
    my $Het  = "$ref$alt";
    my $Homa = "$alt$alt";
    my $Loss = "??";
    my $max = 1;
    $max = length($ref) if(length($ref) > $max);
    $max = length($alt) if(length($alt) > $max);
    return ($Homr, $Het, $Homa, $Loss) if($max == 1);    
    $ref  = "-" x $max if($ref eq "-");
    $alt  = "-" x $max if($alt eq "-");
    $Homr = "$ref$ref";
    $Het  = "$ref$alt";
    $Homa = "$alt$alt";
    $Loss = "?" x ($max*2);
    return ($Homr, $Het, $Homa, $Loss);
}

# 对样本的突变频率绘制小提琴图，用于检查是否有污染
sub plot_violin_plot{
    my $hashSummary  = shift @_;
    my $snv_dir      = shift @_;
    my $data_name    = shift @_;
    my $output_name  = shift @_;
    my $R            = shift @_;
    my $samples      = shift @_;
    my $image_dir    = "$snv_dir/Image";
    package::utils::make_dir($image_dir);

    # 图像名称
    my $pdf_vioplot = "$image_dir/$output_name.pdf";
    my $freq_vioplot = "$image_dir/$output_name.pdf.txt";

    print "Plot $pdf_vioplot";

    # 输出频率文件
    my $sample_count = scalar(@$samples);
    open FREQ, ">$freq_vioplot";
    print FREQ "Sample\tFrequency\n";
    foreach my $sample(@$samples)
    {
        if(not exists $hashSummary->{$data_name}{$sample})
        {
            print FREQ "$sample\t0\n";
            next;
        }
        else
        {
            foreach my $freq(split /,/, $hashSummary->{$data_name}{$sample})
            {
                print FREQ "$sample\t$freq\n";
            }
        }
    }
    close FREQ;

    # 绘图
    my $pdf_width = 2 * $sample_count;
       $pdf_width = 3 * $sample_count   if($sample_count > 5 and $sample_count <= 10);
       $pdf_width = 5 * $sample_count   if($sample_count > 2 and $sample_count <= 5);
       $pdf_width = 7 * $sample_count   if($sample_count <= 2);

    $R->send(qq` mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)\n`); # 颜色模版
    $R->send(qq` mycol <- colors()[rep(mycol, 50)] \n`);

    $R->send(qq` plot_data <- read.table("$freq_vioplot", head = T, colClasses = "character") \n`); # 样本名存在0开头的纯数字格式，R会自动删除，所以按照字符读入，后面再转成数字
    $R->send(qq` plot_data\$Frequency <- as.numeric(plot_data\$Frequency) \n`);
    $R->send(qq` plot_data\$Sample <- factor(plot_data\$Sample) \n`);
    $R->send(qq` sample_name <- unique(plot_data\$Sample) \n`);
    $R->send(qq` sample_color <- mycol[1 : length(sample_name)] \n`);

    $R->send(qq` pdf("$pdf_vioplot", width = $pdf_width, height = 10)\n`);
    $R->send(qq` ggplot(data = plot_data, aes(x = Sample, y = Frequency)) +
        geom_violin(aes(fill = Sample) ) +
        scale_fill_manual(breaks = sample_name, values = sample_color) +
        guides(fill = FALSE) + 
        geom_hline(aes(yintercept = 0.75), colour="#990000", linetype="dashed") + 
        geom_hline(aes(yintercept = 0.5), colour="#990000", linetype="dashed") + 
        geom_hline(aes(yintercept = 0.25), colour="#990000", linetype="dashed") + 
        theme(text = element_text(size = 20)) +
        ggtitle("Sample SNV frequency violin plot") + theme(plot.title = element_text(hjust = 0.5)) \n`);
    $R->send(qq`dev.off()\n`);    

    system("rm $freq_vioplot");
    print "\n";
}

# 柱状图绘制
sub plot_bar{
    my $hashSummary  = shift @_;
    my $snv_dir      = shift @_;
    my $data_name    = shift @_;
    my $output_name  = shift @_;
    my $R            = shift @_;
    my $image_dir    = "$snv_dir/Image";
    package::utils::make_dir($image_dir);

    # 图像名称
    my $png_bar = "$image_dir/$output_name.png";

    print "Plot $png_bar";
    
    # 绘图数据准备
    my @values;
    my @value_names;
    foreach my $value_name(sort keys %{$hashSummary->{$data_name}})
    {
        push @values, $hashSummary->{$data_name}{$value_name};
        push @value_names, "$value_name";
    }
    
    # 判定是否有数据
    if(scalar(@values) == 0)
    {   
        print "[no plot]\n";
        return ;
    }

    my $value_string      = join ",", @values;
    my $value_name_string = join ",", @value_names;

    # 绘图
    $R->send(qq` plot_data <- data.frame(length = c($value_name_string), count = c($value_string))\n`);
    $R->send(qq` png("$png_bar", width = 800, height = 600)\n`);
    $R->send(qq` ggplot(data = plot_data, aes(x = length, y = count)) +
        geom_bar(fill = 'lawngreen', stat = 'identity') +
        theme(text = element_text(size = 20)) +
        ggtitle("InDel length distribution") +
        guides(fill = FALSE) \n`);
    $R->send(qq`dev.off()\n`);

    print "\n";
}

# 饼图绘制
sub plot_pie{
    my $hashSummary  = shift @_;
    my $snv_dir      = shift @_;
    my $data_name    = shift @_;
    my $output_name  = shift @_;
    my $R            = shift @_;
    my $image_dir    = "$snv_dir/Image";
    package::utils::make_dir($image_dir);

    # 图像名称
    my $png_pie = "$image_dir/$output_name.png";

    print "Plot $png_pie";
    # 绘图数据准备
    my @values;
    my @value_names;
    if($data_name eq 'TsTv_Summary') # TsTv_Summary
    {
        my @type_lists = ("Ts", "Tv", "Novel_Ts", "Novel_Tv");
        foreach my $type(@type_lists)
        {
            next if(not exists $hashSummary->{"TsTv_Summary"}{$type});
            push @values, $hashSummary->{"TsTv_Summary"}{$type};
            push @value_names, "\'$type\'";
        }
 
    }
    else # Region_Summary Function_Summary 
    {
        foreach my $value_name(sort keys %{$hashSummary->{$data_name}})
        {   
            next if($value_name eq 'Total');
            push @values, $hashSummary->{$data_name}{$value_name}{'ALL'};
            push @value_names, "\'$value_name\'";
        }
    }

    # 判定是否有数据
    if(scalar(@values) == 0)
    {   
        print "[no plot]\n";
        return ;
    }

    my $value_string      = join ",", @values;
    my $value_name_string = join ",", @value_names;

    # 绘图
    $R->send(qq` mycol <- c(119,132,147,454,89,404,123,529,463,461,128,139,552,28,54,100,258,558,376,43,652,165,31,610,477,256,588,99,632,81,503,104,562,76,96,495,598,645,507,657,33,179,107,62)\n`); # 颜色模版
    $R->send(qq` mycol = colors()[rep(mycol,20)] \n`);

    $R->send(qq` plot_data = data.frame(Group = c($value_name_string), Count = c($value_string))\n`); # 原始数据
    $R->send(qq` plot_data\$Percent       = round(100 * plot_data\$Count / sum(plot_data\$Count), 2) \n`);  # 百分比
    $R->send(qq` plot_data\$Legend        = paste(plot_data\$Group, '(', plot_data\$Percent, '%)', sep = '') \n`); # 标签
    $R->send(qq` plot_data                = plot_data[order(plot_data\$Count, decreasing = TRUE), ] \n`); # 由大到小排序
    $R->send(qq` plot_data\$Group         = factor(plot_data\$Group, levels = plot_data\$Group) \n`); # factor定义顺序
    $R->send(qq` plot_data\$Color_Defined = mycol[0 : nrow(plot_data)] \n`); # 颜色

    $R->send(qq` png("$png_pie", width = 800, height = 600)\n`);
    $R->send(qq` ggplot(plot_data, aes( x = "", y = Count, fill = Group )) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar( theta = "y" ) +
        labs( x = "", y = "", title = "") +
        theme(text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.title = element_blank() , legend.position = "right" ) +
        scale_fill_manual(breaks = plot_data\$Group, values = plot_data\$Color_Defined, labels = plot_data\$Legend) \n`);   # 将原来的图例标签换成现在的myLabel
    $R->send(qq`dev.off()\n`);

    print "\n";
}


sub output_inbreeding_coef{
    my $workbook     = shift @_;
    my $format       = shift @_;
    my $hashInb      = shift @_;
    my $samples      = shift @_;

    print " [Inbreeding coefficients]\n";
    my $sheet = $workbook->add_worksheet("Inbreeding coef");
       $sheet->set_row(0, 60);
       $sheet->set_column(0, 0, 20);
    
    # 表头
    my @titles = ("FID","IID","O(HOM)","E(HOM)","N(NM)","F");
    my $row = 0;
    foreach my $col(0..$#titles)
    {
        $sheet->write($row, $col, $titles[$col], $format->{"title"}); 
    }
    $row++;
    
    # 数据输出
    foreach my $sample(@$samples)
    {
        my $info  = $hashInb->{$sample};
        my @datas = split /\t/,$info;
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;  
    }
}

# InDel Length突变数量统计
sub output_indel_length{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $hashSummary = shift @_;
    my $samples     = shift @_;

    print " [InDel Length]\n";
    my $sheet = $workbook->add_worksheet("InDel Length");
       $sheet->set_row(0, 60);
       $sheet->set_column(0, 0, 20);

    my @titles = ('Insertion/Deletion Length', 'InDel Count', 'Percent');
    my @indel_length_lists = ("5~1",
        "10~6",
        "15~11",
        "20~16",
        ">20",
        "-5~-1",
        "-10~-6",
        "-15~-11",
        "-20~-16",
        "<-20",
        ); 
    # 表头
    my $row = 0;
    foreach my $col(0..$#titles){
        $sheet->write($row, $col, $titles[$col], $format->{"title"}); 
    }
    $row++;

    foreach my $indel_length(@indel_length_lists){
        my @datas;
        foreach my $title(@titles)
        {
            my $value = "";
               $value = exists $hashSummary->{"Indel_Length_Summary"}{$indel_length}{$title} ? $hashSummary->{"Indel_Length_Summary"}{$indel_length}{$title} : "NA";
               $value = $indel_length if($title eq 'Insertion/Deletion Length');

            push @datas, $value;
        }

        # 数据输出
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;  
    }

}

# Ts Tv突变数量统计
sub output_ts_vs_tv{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $hashSummary = shift @_;
    my $samples     = shift @_;

    print " [Ts Tv]";
    my $sheet = $workbook->add_worksheet("Ts Tv");
       $sheet->set_row(0, 60);
       $sheet->set_column(0, 0, 20);

    my @titles = ('SNV Type','Count');
    my @type_lists = ("Ts",
        "Tv",
        "Ts/Tv",
        "Novel_Ts",
        "Novel_Tv",
        "Novel_Ts/Tv",
        ); 

    # 表头
    my $row = 0;
    foreach my $col(0..$#titles){
        $sheet->write($row, $col, $titles[$col], $format->{"title"}); 
    }
    $row++;

    foreach my $type(@type_lists){
        my $type_value = exists $hashSummary->{'TsTv_Summary'}{$type} ? $hashSummary->{'TsTv_Summary'}{$type} : "";
        my @datas = ($type, $type_value);

        # 数据输出
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;  
    }

}

# HomHetNovel突变数量统计
sub output_hom_het_novel{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $hashSummary = shift @_;
    my $samples     = shift @_;

    print " [HomHetNovel]";
    my $sheet = $workbook->add_worksheet("HomHetNovel");
       $sheet->set_row(0, 60);
       $sheet->set_column(0, 0, 20);

    my @titles = ('Sample', 'Het_SNV', 'Hom_SNV', 'Novel_SNV', 'Het_InDel', 'Hom_InDel', 'Novel_InDel', 'Het_SNV/Hom_SNV');
    my $row = 0;
    foreach my $col(0..$#titles){
        $sheet->write($row, $col, $titles[$col], $format->{"title"}); 
    }
    $row++;

    foreach my $sample(@$samples){
        my @datas;
        foreach my $title(@titles)
        {
            my $value = "";
               $value = exists $hashSummary->{"HomHetNovel_SNVInDel_Summary"}{$title}{$sample} ? $hashSummary->{"HomHetNovel_SNVInDel_Summary"}{$title}{$sample} : 0;
               $value = $sample if($title eq 'Sample');

            push @datas, $value;
        }

        # 数据输出
        $sheet->write_string($row, 0, $datas[0], $format->{'normal'});
        foreach my $col(1..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;  
    }

    # 行业内部质控警告输出
    my $string = "Warnings : Het_SNV/Hom_SNV > 1.6 : ";
    my $flag = 0;
    foreach my $sample(@$samples){
        my $value = exists $hashSummary->{"HomHetNovel_SNVInDel_Summary"}{"Het_SNV/Hom_SNV"}{$sample} ? $hashSummary->{"HomHetNovel_SNVInDel_Summary"}{"Het_SNV/Hom_SNV"}{$sample} : "NA";
        if($value=~/\d/ and $value > 1.6){
            $string.="$sample,";
            $flag = 1;
        }        
    }
    system "echo -e '\\033[41;37m $string \\033[0m'" if($flag==1);

}

# GATK突变数量统计
sub output_gatk{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $hashSummary = shift @_;
    my $samples     = shift @_;

    print " [GATK]";
    my $sheet = $workbook->add_worksheet("GATK");
       $sheet->set_row(0, 60);
       $sheet->set_column(0, 0, 20);

    # 表头
    $sheet->merge_range(0, 0, 1, 0, "Sample", $format->{"title"});
    $sheet->merge_range(0, 1, 0, 2, "GATK", $format->{"title"});
    $sheet->write(1, 1, 'SNV', $format->{"title"});
    $sheet->write(1, 2, 'InDel', $format->{"title"});

    my $row = 2;
    foreach my $sample(("ALL", @$samples))
    {   
        my $snv_count   = exists $hashSummary->{"GATK_Summary"}{"SNV"}{$sample}   ? $hashSummary->{"GATK_Summary"}{"SNV"}{$sample}   : 0;
        my $indel_count = exists $hashSummary->{"GATK_Summary"}{"InDel"}{$sample} ? $hashSummary->{"GATK_Summary"}{"InDel"}{$sample} : 0;
        my @datas = ($sample, $snv_count, $indel_count);
       
        # 数据输出
        $sheet->write_string($row, 0, $datas[0], $format->{'normal'});
        foreach my $col(1..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;        
    }

}

# 基因功能统计
sub output_function_count{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $hashSummary = shift @_;
    my $samples     = shift @_;

    print " [function]";
    my $sheet = $workbook->add_worksheet("Function Count");
       $sheet->set_row(0, 60);
       $sheet->set_column(0, 0, 20);
    my @titles         = ("Region", "ALL", @$samples);
    my @function_lists = ("nonsynonymous SNV",
        "synonymous SNV",
        "unknown",
        "nonframeshift deletion",
        "stopgain",
        "frameshift deletion",
        "nonframeshift insertion",
        "frameshift insertion",
        "stoploss",
        "Total",
        ); 

    # 检查是否有列表以外的区域信息
    my @unknowns;
    foreach my $function(keys %{$hashSummary->{"Function_Summary"}})
    {
        next if($function ~~ @function_lists);
        push @unknowns, $function;
    }
    print "[Waring] Find Unknown function value : @unknowns\n" if(@unknowns > 0);


    # 表头
    my $row = 0;
    foreach my $col(0..$#titles)
    {
        $sheet->write_string($row, $col, $titles[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $function(@function_lists)
    {   
        # 数据准备
        my @datas = ($function);
        foreach my $sample(("ALL", @$samples))
        {
            my $value = (exists $hashSummary->{"Function_Summary"}{$function} and exists $hashSummary->{"Function_Summary"}{$function}{$sample}) ? $hashSummary->{"Function_Summary"}{$function}{$sample} : 0;
            push @datas, $value;
        }

        # 数据输出
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;
    }

}

# 基因区域统计
sub output_gene_region_count{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $hashSummary = shift @_;
    my $samples     = shift @_;

    print " [region]";
    my $sheet = $workbook->add_worksheet("Gene Region Count");
       $sheet->set_row(0, 60);
       $sheet->set_column(0, 0, 20);
    my @titles       = ("Region", "ALL", @$samples);
    my @region_lists = ("intronic",
        "exonic",
        "intergenic",
        "ncRNA_intronic",
        "UTR3","splicing",
        "ncRNA_exonic",
        "UTR5",
        "upstream",
        "downstream",
        "ncRNA_splicing",
        "upstream;downstream",
        "exonic;splicing",
        "ncRNA_exonic;splicing", 
        "UTR5;UTR3",
        "ncRNA_UTR5",
        "ncRNA_UTR3",
        "Total",
        ); 

    # 检查是否有列表以外的区域信息
    my @unknowns;
    foreach my $region(keys %{$hashSummary->{"Region_Summary"}})
    {
        next if($region ~~ @region_lists);
        push @unknowns, $region;
    }
    print "[Waring] Find Unknown region value : @unknowns\n" if(@unknowns > 0);

    # 表头
    my $row = 0;
    foreach my $col(0..$#titles)
    {
        $sheet->write_string($row, $col, $titles[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $region(@region_lists)
    {   
        # 数据准备
        my @datas = ($region);
        foreach my $sample(("ALL", @$samples))
        {
            my $value = (exists $hashSummary->{"Region_Summary"}{$region} and exists $hashSummary->{"Region_Summary"}{$region}{$sample}) ? $hashSummary->{"Region_Summary"}{$region}{$sample} : 0;
            push @datas, $value;
        }

        # 数据输出
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;
    }

}

# 突变统计汇总
sub summary_snv{
    my $hashSNV  = shift @_;
    my $hashGeno = shift @_;
    my $samples  = shift @_;

    my %hashSummary; # 统计结果保存

    my @snv_lists = keys %$hashSNV;
    my $snv_count = @snv_lists;
    my $count     = 0;

    my %hashTs = ('AG'=>'Ts','TC'=>'Ts','GA'=>'Ts','CT'=>'Ts'); # 转换列表
    foreach my $snv_title(@snv_lists)
    {
        $count++;
        my $process = sprintf "%0.2f", 100 * $count / $snv_count;
        print "\rStatistic : $process%";

        foreach my $alt(keys %{$hashSNV->{$snv_title}})
        {   
            my $ref = $hashSNV->{$snv_title}{$alt}{"Ref Allele"};
            # 读取样本分型
            my %hashGeno_this_snv;
            read_sample_geno(\%hashGeno_this_snv, $hashGeno->{$snv_title}{$alt}{"Mutation 0"}); 
            read_sample_geno(\%hashGeno_this_snv, $hashGeno->{$snv_title}{$alt}{"Mutation 1"});
            read_sample_geno(\%hashGeno_this_snv, $hashGeno->{$snv_title}{$alt}{"Mutation 2"});
 
            # 获取突变样本名称
            my @mut_samples = grep{ exists $hashGeno_this_snv{$_} and ($hashGeno_this_snv{$_}{"Type"} eq 'HET' or $hashGeno_this_snv{$_}{"Type"} eq 'HOMA') } @$samples;
 
            # 没有突变，跳过
            next if(scalar(@mut_samples) == 0);

            # 区域统计
            my $alt_region = exists $hashSNV->{$snv_title}{$alt}{"Gene Region"} ? $hashSNV->{$snv_title}{$alt}{"Gene Region"} : "";
            if($alt_region =~ /\w/)
            {
                $hashSummary{"Region_Summary"}{$alt_region}{'ALL'}++;
                $hashSummary{"Region_Summary"}{"Total"}{'ALL'}++;
                foreach my $sample(@mut_samples)
                {
                    $hashSummary{"Region_Summary"}{$alt_region}{$sample}++;
                    $hashSummary{"Region_Summary"}{"Total"}{$sample}++;
                }           
            }

            # 功能统计
            my $alt_function = exists $hashSNV->{$snv_title}{$alt}{"Function"} ? $hashSNV->{$snv_title}{$alt}{"Function"} : "";
            if($alt_function =~ /\w/)
            {
                $hashSummary{"Function_Summary"}{$alt_function}{'ALL'}++;
                $hashSummary{"Function_Summary"}{"Total"}{'ALL'}++;
                foreach my $sample(@mut_samples)
                {
                    $hashSummary{"Function_Summary"}{$alt_function}{$sample}++;
                    $hashSummary{"Function_Summary"}{"Total"}{$sample}++;
                }           
            }

            # GATK SNV/indel 突变数量汇总
            my $mutation_type = 'SNV'; 
               $mutation_type = 'InDel' if($ref eq "-" or $alt eq "-");
            $hashSummary{"GATK_Summary"}{$mutation_type}{'ALL'}++;
            foreach my $sample(@mut_samples)
            {
                $hashSummary{"GATK_Summary"}{$mutation_type}{$sample}++;
            }

            # 杂合纯合统计
            my $snp_id = $hashSNV->{$snv_title}{$alt}{"SNP ID"};
            foreach my $sample(@mut_samples)
            {
                my $hom_het_novel = 'Novel';
                   $hom_het_novel = "Het" if($snp_id =~ /rs\d+/ and $hashGeno_this_snv{$sample}{"Type"} eq 'HET');
                   $hom_het_novel = "Hom" if($snp_id =~ /rs\d+/ and $hashGeno_this_snv{$sample}{"Type"} eq 'HOMA');

                $hashSummary{"HomHetNovel_SNVInDel_Summary"}{"$hom_het_novel\_$mutation_type"}{$sample}++;
            }
            
            # Ts Tv统计
            if($ref ne "-" and $alt ne "-")
            {
                my $alt_type = exists $hashTs{"$ref$alt"} ? "Ts" : 'Tv';
                my $alt_mark = ($snp_id =~ /rs\d+/) ? "$alt_type" : "Novel_$alt_type";
                $hashSummary{'TsTv_Summary'}{$alt_mark}++;        
            }

            # 插入缺失长度统计
            if($ref eq "-" or $alt eq "-")
            {
                my $indel_length = ($ref eq '-') ? length($alt) : length($ref);
                my $in_or_del    = ($ref eq '-') ? '' : '-'; # - 表示缺失
                $hashSummary{"Indel_Length_Summary"}{$in_or_del . "5~$in_or_del" . "1"}{'InDel Count'}++   if($indel_length >= 1  and $indel_length <= 5);
                $hashSummary{"Indel_Length_Summary"}{$in_or_del . "10~$in_or_del" . "6"}{'InDel Count'}++  if($indel_length >= 6  and $indel_length <= 10);
                $hashSummary{"Indel_Length_Summary"}{$in_or_del . "15~$in_or_del" . "11"}{'InDel Count'}++ if($indel_length >= 11 and $indel_length <= 15);
                $hashSummary{"Indel_Length_Summary"}{$in_or_del . "20~$in_or_del" . "16"}{'InDel Count'}++ if($indel_length >= 16 and $indel_length <= 20);
                $hashSummary{"Indel_Length_Summary"}{">20"}{'InDel Count'}++                           if($indel_length > 20 and $in_or_del eq '');                     
                $hashSummary{"Indel_Length_Summary"}{"<-20"}{'InDel Count'}++                          if($indel_length > 20 and $in_or_del eq '-');                     
                $hashSummary{"Indel_Length_Summary_Plot"}{"$in_or_del$indel_length"}++; # 绘图数据                     
            }

            # 突变频率统计, 绘制小提琴图，检测样品是否有污染, 只统计可靠位点
            if($hashSNV->{$snv_title}{$alt}{"HOMOLOGY HITS"} <= 1 and $hashSNV->{$snv_title}{$alt}{"SNP ID"} !~ /STR|INDEL/) 
            {
                foreach my $sample(@mut_samples)
                {   
                    $hashSummary{"Alt_Freq_Summary"}{$sample} .= "$hashGeno_this_snv{$sample}{'Freq'},";
                }
            } 
        }
    }

    # HET/HOMA 比例, 根据经验，外显子项目中 该比例> 1.6 ，则存在污染
    if(exists $hashSummary{"HomHetNovel_SNVInDel_Summary"})
    {   
        foreach my $sample(@$samples)
        {
            my $Het_SNV = (exists $hashSummary{"HomHetNovel_SNVInDel_Summary"}{'Het_SNV'} and exists $hashSummary{"HomHetNovel_SNVInDel_Summary"}{'Het_SNV'}{$sample}) ? $hashSummary{"HomHetNovel_SNVInDel_Summary"}{'Het_SNV'}{$sample} : 0;
            my $Hom_SNV = (exists $hashSummary{"HomHetNovel_SNVInDel_Summary"}{'Hom_SNV'} and exists $hashSummary{"HomHetNovel_SNVInDel_Summary"}{'Hom_SNV'}{$sample}) ? $hashSummary{"HomHetNovel_SNVInDel_Summary"}{'Hom_SNV'}{$sample} : 0;
            my $het_vs_hom = ($Hom_SNV > 0) ? $Het_SNV/$Hom_SNV : 'NA';

            $hashSummary{"HomHetNovel_SNVInDel_Summary"}{'Het_SNV/Hom_SNV'}{$sample} = $het_vs_hom;
        }       
    }

    # TS TV 比例
    if(exists $hashSummary{'TsTv_Summary'})
    {
        my $ts_count       = exists $hashSummary{'TsTv_Summary'}{'Ts'}       ? $hashSummary{'TsTv_Summary'}{'Ts'}       : 0;
        my $tv_count       = exists $hashSummary{'TsTv_Summary'}{'Tv'}       ? $hashSummary{'TsTv_Summary'}{'Tv'}       : 0;
        my $ts_novel_count = exists $hashSummary{'TsTv_Summary'}{'Novel_Ts'} ? $hashSummary{'TsTv_Summary'}{'Novel_Ts'} : 0;
        my $tv_novel_count = exists $hashSummary{'TsTv_Summary'}{'Novel_Tv'} ? $hashSummary{'TsTv_Summary'}{'Novel_Tv'} : 0;

        my $ts_vs_tv       = ($tv_count > 0)       ? sprintf "%0.4f", $ts_count       / $tv_count       : "NA";
        my $ts_vs_tv_novel = ($tv_novel_count > 0) ? sprintf "%0.4f", $ts_novel_count / $tv_novel_count : "NA";

        $hashSummary{'TsTv_Summary'}{"Ts/Tv"}       = $ts_vs_tv;
        $hashSummary{'TsTv_Summary'}{"Novel_Ts/Tv"} = $ts_vs_tv_novel;
    }

    # Indel_Length_Summary 百分比计算
    if(exists $hashSummary{"Indel_Length_Summary"})
    {
        my $sum = 0;
        map{$sum += $hashSummary{"Indel_Length_Summary"}{$_}{'InDel Count'}; } keys %{$hashSummary{"Indel_Length_Summary"}};

        foreach my $indel_length(keys %{$hashSummary{"Indel_Length_Summary"}})
        {
            $hashSummary{"Indel_Length_Summary"}{$indel_length}{"Percent"} = sprintf "%0.2f", 100 * $hashSummary{"Indel_Length_Summary"}{$indel_length}{'InDel Count'} / $sum;
            $hashSummary{"Indel_Length_Summary"}{$indel_length}{"Percent"} .= "%";
        }
    }

    print "\n";

    return %hashSummary;

} 

# 获取样本基因分型信息
sub read_sample_geno{
    my $hashGeno_this_snv = shift @_;
    my $geno_info         = shift @_;

    foreach my $sample_snv_info((split ',', $geno_info))
    {   
        next if($sample_snv_info !~ /\w/);
        my ($sample, $type, $geno, $ref_depth, $alt_depth, $tmp) = split /:/, $sample_snv_info, 6;
        $hashGeno_this_snv->{$sample}{"Type"} = $type;
        $hashGeno_this_snv->{$sample}{"Geno"} = $geno;
        $hashGeno_this_snv->{$sample}{"Freq"} = sprintf "%0.4f", $alt_depth / ($ref_depth + $alt_depth);
    }   
}

sub readLibrary{
    my $file = shift @_;
    my %library = ();
    print "Reading $file...";
    open FILE,$file;
    while(my $line = <FILE>){
        $line =~ s/[\r\n]//g;
        my ($title,$q,$rs,$freq,$ref,$alt,$chr,$pos,$gene,$region,$function,$hgvs,$seq5,$seq3,$blast) = split /\t/,$line;
        $library{$title}{$alt}{"SNP Calling Quality"} = $q;
        $library{$title}{$alt}{"SNP ID"}              = $rs;
        $library{$title}{$alt}{"Freq_Alt (1000g)"}    = $freq;
        $library{$title}{$alt}{"Ref Allele"}          = $ref;
        $library{$title}{$alt}{"Alt Allele"}          = $alt;
        $library{$title}{$alt}{"Chrs"}                = $chr;
        $library{$title}{$alt}{"Position"}            = $pos;
        $library{$title}{$alt}{"Position2"}           = (split /\-/,$pos)[0];
        $library{$title}{$alt}{"Gene"}                = $gene;
        $library{$title}{$alt}{"Gene Region"}         = $region;
        $library{$title}{$alt}{"Function"}            = $function;
        $library{$title}{$alt}{"Predicted Protein Variants"} = $hgvs;
        $library{$title}{$alt}{"5' FLANKING SEQUENCE"}       = $seq5;
        $library{$title}{$alt}{"3' FLANKING SEQUENCE"}       = $seq3;
        $library{$title}{$alt}{"HOMOLOGY HITS"}              = $blast;
    }
    close FILE;
    print "OK\n";
    return %library;
}

sub readAnalysis{
    my $file = shift @_;
    my %analysis = ();
    print "Reading $file...";
    open FILE,$file;
    while(my $line = <FILE>){
        $line =~ s/[\r\n]//g;
        my ($title,$alt,$genoQ,$s0,$s1,$s2) = split /\t/,$line;
        $analysis{$title}{$alt}{"Genotyping Quality"} = $genoQ;
        $analysis{$title}{$alt}{"Mutation 0"}         = $s0;
        $analysis{$title}{$alt}{"Mutation 1"}         = $s1;
        $analysis{$title}{$alt}{"Mutation 2"}         = $s2;
    }
    close FILE;
    print "OK\n";
    return %analysis;
}

# fastqc拷贝
sub prepare_QC_image{
    my $hashConfig = shift @_;
    my $samples    = shift @_;
    my $report_dir     = $hashConfig->{'Report'};
    my $document_dir   = "$report_dir/document";
    my $qc_dir         = "$document_dir/1_QC";
    my $fastqc_cp_dir  = "$qc_dir/fastqc";  # 拷贝路径
    my $fastqc_dir     = "$report_dir/fastqc"; # 原图路径

    package::utils::make_dir($document_dir);
    package::utils::make_dir($qc_dir);
    package::utils::make_dir($fastqc_cp_dir);

    print "Copy fastqc image ... ";
    foreach my $sample(@$samples){
        my $R1QC       = "$fastqc_dir/$sample\_R1_fastqc/Images/per_base_quality.png";
        my $R1ATCG     = "$fastqc_dir/$sample\_R1_fastqc/Images/per_base_sequence_content.png";
        my $R2QC       = "$fastqc_dir/$sample\_R2_fastqc/Images/per_base_quality.png";
        my $R2ATCG     = "$fastqc_dir/$sample\_R2_fastqc/Images/per_base_sequence_content.png";
        my $Errorate   = "$fastqc_dir/$sample\_ErrorRate.png";
        my $InsertSize = "$report_dir/status/$sample.insert.pdf";
        system("cp $R1QC       $fastqc_cp_dir/$sample\_R1_per_base_quality.png") if(-e $R1QC);
        system("cp $R2QC       $fastqc_cp_dir/$sample\_R2_per_base_quality.png") if(-e $R2QC);
        system("cp $R1ATCG     $fastqc_cp_dir/$sample\_R1_per_base_sequence_content.png") if(-e $R1ATCG);
        system("cp $R2ATCG     $fastqc_cp_dir/$sample\_R2_per_base_sequence_content.png") if(-e $R2ATCG);
        system("cp $Errorate   $fastqc_cp_dir/") if(-e $Errorate);
        system("cp $InsertSize $fastqc_cp_dir/") if(-e $InsertSize);
    }
    print "OK\n";
}

# fastqc clean 拷贝
sub prepare_QC_Clean_image{
    my $hashConfig = shift @_;
    my $samples    = shift @_;
    my $report_dir     = $hashConfig->{'Report'};
    my $document_dir   = "$report_dir/document";
    my $qc_dir         = "$document_dir/1_QC";
    my $fastqc_cp_dir  = "$qc_dir/fastqc_clean";  # 拷贝路径
    my $fastqc_dir     = "$report_dir/fastqc_clean"; # 原图路径

    return if(package::utils::is_dir_ok($fastqc_dir) == 0);

    package::utils::make_dir($document_dir);
    package::utils::make_dir($qc_dir);
    package::utils::make_dir($fastqc_cp_dir);

    print "Copy fastqc clean image ... ";
    foreach my $sample(@$samples){
        my $R1QC       = "$fastqc_dir/$sample.trim.R1_fastqc/Images/per_base_quality.png";
        my $R1ATCG     = "$fastqc_dir/$sample.trim.R1_fastqc/Images/per_base_sequence_content.png";
        my $R2QC       = "$fastqc_dir/$sample.trim.R2_fastqc/Images/per_base_quality.png";
        my $R2ATCG     = "$fastqc_dir/$sample.trim.R2_fastqc/Images/per_base_sequence_content.png";
        my $Errorate   = "$fastqc_dir/$sample\_ErrorRate.png";
        system("cp $R1QC       $fastqc_cp_dir/$sample\_R1_per_base_quality.png") if(-e $R1QC);
        system("cp $R2QC       $fastqc_cp_dir/$sample\_R2_per_base_quality.png") if(-e $R2QC);
        system("cp $R1ATCG     $fastqc_cp_dir/$sample\_R1_per_base_sequence_content.png") if(-e $R1ATCG);
        system("cp $R2ATCG     $fastqc_cp_dir/$sample\_R2_per_base_sequence_content.png") if(-e $R2ATCG);
        system("cp $Errorate   $fastqc_cp_dir/") if(-e $Errorate);
    }    
    print "OK\n";
}

1