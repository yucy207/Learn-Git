package package::region_plot;
use strict;
use warnings;
$|=1;

sub run{
    my $hashPara   = shift @_;
    my $hashConfig = shift @_;
    my $SVdbtype   = exists $hashConfig->{'GenomeBed'} ? "SVdbWGS" : "SVdbWES";
    return if ($SVdbtype eq "SVdbWES");
    print "########## Start Region plot ".package::utils::get_time()." ##########\n";
    my $report     = $hashConfig->{'Report'};
    my $bedtools   = $hashPara->{'Soft'}{'bedtools'};
    my $document   = "$report/document"; package::utils::make_dir($document);
    my $sv_dir     = "$document/4_SV"; package::utils::make_dir($sv_dir);
    my $outDir     = "$sv_dir/region_plot"; package::utils::make_dir($outDir);
    my $pdf        = "$outDir/Region_plot.pdf";
    if(package::utils::is_file_ok($pdf))
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }    
    my @cases      = package::utils::get_sample($hashConfig, "case");
    my $plotDirTemp = "$report/sv/region_plot"; package::utils::make_dir($plotDirTemp);  #临时文件存储目录，后续需要删除
    
    # CNVseq结果与SV cnv结果 取交集
    my $cnvSeqDir = "$report/cnvseq";
    if (@cases > 0) {
        my $threshold = 5;
        my $pm = Parallel::ForkManager->new($threshold);
        foreach my $sample (@cases) {
            $pm -> start($sample) and next;
            getIntersect ($bedtools, $sample, $cnvSeqDir, $plotDirTemp);
            $pm -> finish;    
        }
        $pm->wait_all_children;
   }
   else 
   {
        print "[Note] This Project Don't have Case Sample! Please Check ConfigFile\n";
   }  
   
   # 生成区域绘图文件 
   my $textDir     = "$outDir/region_text";  package::utils::make_dir($textDir);
   my %hashSamples = generatePloaData($plotDirTemp, $textDir, \@cases); 
   
   # 绘图  
   plot ($pdf, $hashPara, $hashConfig, \%hashSamples, \@cases);
   `rm -r $plotDirTemp`;
}

sub plot{
    my $pdf         = shift @_;
    my $hashPara    = shift @_;
    my $hashConfig  = shift @_;
    my $hashSamples = shift @_;
    my $cases       = shift @_;
    my $samplecount = @$cases;
    my $pdf_heigth  = 3.2 * $samplecount + 1.7;
    my $height      = 500;
    my $Rbin        = $hashPara -> {'Soft'}{'R'};
    my $Rlib        = $hashPara -> {'Soft'}{'RLib'};
    my $R           = Statistics::R->new(bin => "$Rbin");
    $R -> startR;
    $R -> send(qq` .libPaths("$Rlib") \n`);
    $R -> send(qq` pdf("$pdf", width=18, height=$pdf_heigth) \n`);
    foreach my $region (sort keys %$hashSamples){ print "Plot $region ...";
        my ($chr, $start, $end) = split /[:-]/, $region;
        my $file       = $hashConfig  -> {'Report'}."/document/4_SV/region_plot/region_text/Chr$chr-$start-$end.txt";
        my $type       = $hashSamples -> {$region}{$$cases[0]}{'Type'};
        my $title      = "$chr\_$start\_$end\_$type";
        $R -> send(qq` par(mfrow=c($samplecount,1), mar=c(5,8,3,5), oma=c(0,0,4,0), ps=22) \n`);
        $R -> send(qq` data=read.table("$file", header=T, check.names = F) \n`);
        foreach my $sample (@$cases)
        {
            my $cnv = ($sample eq 'RefSample') ? "" : "_CN".$hashSamples -> {$region}{$sample}{'CNV'}; 
            $R -> send(qq` Position=data["Position"][[1]] \n`);
            $R -> send(qq` Depth=data["$sample"][[1]] \n`);
            $R -> send(qq` quan = quantile(Depth) \n`);
            $R -> send(qq` top_cutoff = quan[4] + 3 * (quan[4] - quan[2]) \n`);
            $R -> send(qq` plot(Position,Depth,type='l',ylim=c(0,top_cutoff),main="$sample$cnv",cex.main=0.75,cex.lab=0.75,cex.axis=0.75) \n`);
            $R -> send(qq` abline(v=c($start,$end),lwd=2,col='red') \n`);           
        }   
        $R -> send(qq` mtext("$title",side=3,line=0,outer=T) \n`); print "ok\n";
    }
    $R -> send(qq` dev.off() \n`);
    $R -> stop();
}

sub generatePloaData {
    my $plotDirTemp  = shift @_;
    my $textDir      = shift @_;
    my $cases        = shift @_;
    my %hashdata     = ();
    my %hashSamples  = ();
    foreach my $sample (@$cases){
        my $bed = "$plotDirTemp/$sample\_intersect.bed";
        print "Reading $bed ...";
        open IN, $bed;
        while(<IN>){
            my ($chr, $start, $end, $regionstart, $regionend, $type, $cnv, $cnvseqChr, $cnvseqStart, $cnvseqend, $test, $ref, $tmp)=split /\t/, $_, 13;
            my $region = "$chr\-$regionstart\-$regionend";          
            if(not exists $hashdata{$region}{$cnvseqStart}) {
                $hashdata{$region}{$cnvseqStart}{'RefSample'} = $ref;   
                $hashdata{$region}{$cnvseqStart}{'Start'}     = $start;
                $hashdata{$region}{$cnvseqStart}{'End'}       = $end;
            }
            $hashdata{$region}{$cnvseqStart}{$sample}         = $test;              
            $hashSamples{$region}{$sample}{'Type'}            = $type;
            $hashSamples{$region}{$sample}{'CNV'}             = $cnv;
        }
        close IN;
        print " OK \n";
    }
    push @$cases,   "RefSample";
    foreach my $region (sort keys %hashdata)
    {
         my ($chr,$start,$end) = split /[:-]/,$region;
         my $file = "$textDir/Chr$chr-$start-$end.txt";   
         open OUT, ">$file";
         my $samles = join "\t", @$cases;
         print OUT "Position\t$samles\n";
         foreach my $position (sort {$a<=>$b} keys %{$hashdata{$region}}){ 
             print OUT "$position";
             foreach my $sample (@$cases){
                 my $depth = exists $hashdata{$region}{$position}{$sample} ? $hashdata{$region}{$position}{$sample} : 0;
                 print OUT "\t$depth";   
             }
             print OUT "\n";                        
         }
         close OUT;         
    }
    return %hashSamples;

}

sub getIntersect{
    my $bedtools    = shift @_;
    my $sample      = shift @_;
    my $cnvSeqDir   = shift @_;
    my $plotDirTemp = shift @_;    
    my $cnvFile     = (glob "$cnvSeqDir/$sample.*.cnv")[0];
    my $cnvBed      = "$plotDirTemp/$sample\_cnv.bed";
    open CNV, $cnvFile;   
    open BED, ">$cnvBed";
    <CNV>;
    while(<CNV>) {
        my ($chr, $infos) = split /\t/, $_, 2;
        $chr =~ s/[\"]//g;
        print BED "$chr\t$infos";
    }
    close CNV;
    close BED;
    my $regionBed    = "$plotDirTemp/$sample\_region.bed";
    my $intersectBed = "$plotDirTemp/$sample\_intersect.bed";
    system("$bedtools intersect -loj -a $regionBed -b $cnvBed > $intersectBed ");

} 

1
