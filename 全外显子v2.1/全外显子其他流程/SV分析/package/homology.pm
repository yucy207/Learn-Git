package package::homology;
use strict;
use warnings;

sub run{
    my $hashPara     = shift @_;
    my $hashConfig   = shift @_;
    my $SVdbtype     = (exists $hashConfig->{'GenomeBed'}) ? "SVdbWGS" : "SVdbWES"; 
    my $HomologyDir  = $hashConfig->{'Report'}."/sv/homology"; package::utils::make_dir($HomologyDir);
    my $HomologyFile = "$HomologyDir/homology.anno";
    print "########## Start Homology ".package::utils::get_time()." ##########\n";
    if(package::utils::is_file_ok($HomologyFile))
    {
        print "[Note] Process None, for reProcess Delete Result\n";
        return;
    }    
    my $species      = package::utils::get_species($hashConfig);
    my $genome       = $hashPara->{$species}{'Genome'};
    my $blast_idx    = $hashPara->{$species}{'blast_idx'};
    my $svFile       = $hashConfig->{'Report'}."/sv/sample.all.final.sv.vcf";
    my %hashGenome   = package::utils::readgenome($genome);
    my %hashSV       = package::utils::readSV($svFile, $hashConfig, $SVdbtype);
    
    ########
    # 记录断点，同源性初始化
    ########
    my %hashHomology = ();
    foreach my $title(sort keys %hashSV){
        my ($chr1, $pos1, $chr2, $pos2, $svtype) = split /\|/,$title;
        $hashHomology{"$chr1|$pos1"} = 0;
        $hashHomology{"$chr2|$pos2"} = 0;
    }
    ########
    #提取断点处的序列
    ########
    my $breakPointSeqFile ="$HomologyDir/breakPointSeq.fa";# 断点处，左右扩展50bp参考基因组序列
    open BREAKPOINTSEQ,">$breakPointSeqFile";
    foreach my $point(keys %hashHomology){
        my ($chr,$position)=split /\|/,$point;
        my $start = $position-50;
        my $end   = $position+50;
        my $seq   = substr($hashGenome{$chr},$start-1,$end-$start+1);
        print BREAKPOINTSEQ ">$point\n$seq\n";
    }
    close BREAKPOINTSEQ;
    ########
    # megablast比对,并确认同源性
    ########
    print "Generate $breakPointSeqFile.blast ...";
    system $hashPara->{"Soft"}{"BlastDIR"}."/megablast -a 20 -d $blast_idx -i $breakPointSeqFile -o $breakPointSeqFile.blast -F F -D3 -p 85 -s 50 ";
    print " ok\n";
    open BLAST,"$breakPointSeqFile.blast";
    while (<BLAST>){
        next if(/^#/);
        my ($name,$chr,$identity,$length,$temp)=split(/\t/,$_,5);
        $hashHomology{$name}++ if($length>90);
    }
    close BLAST;
    ########
    # 同源性结果输出
    ########
    open OUT,">$HomologyFile";
    foreach my $title(keys %hashHomology){
        my ($chr, $pos) = split /\|/, $title;
        my $count = $hashHomology{$title};
        print OUT "$chr\t$pos\t$count\n";        
    }
    close OUT;
}

1