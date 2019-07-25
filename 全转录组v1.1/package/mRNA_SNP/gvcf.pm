package package::mRNA_SNP::gvcf;
use Parallel::ForkManager;
use strict;
use warnings;
$|=1;

sub run
{

    my $metadata = shift;
    my $base     = shift;
    my @samples  = split /,/, $metadata->{'samples'};
    my $organ    = qq{$metadata->{organ}};

    my $Java       = qq{$base->{Java}};
    my $GATK4_Loc  = qq{$base->{GATK4_Loc}};
    my $Genome     = qq{$base->{$organ}{genome_fasta}};
    my $RealignBed = qq{$base->{$organ}{RealignBed}};

    my $map      = qq{$metadata->{project}/mapping/result};
    my $result   = qq{$metadata->{project}/mRNA/snp}; 
    my $report   = qq{$metadata->{report}/03_mRNA_Analysis/09_mRNA_SNP_Analysis};

    system qq{mkdir -p $report}         if not -d qq{$report};
    system qq{mkdir -p $result/run}     if not -d qq{$result/run};
    system qq{mkdir -p $result/result}  if not -d qq{$result/result};
    system qq{mkdir -p $result/log}     if not -d qq{$result/log};

    my @res_samples = res_check(qq{$result/result}, \@samples);

    if (scalar @res_samples == 0) {
        print qq{gvcf已经分析完成!\n};
        return 0;
    }

    my $max_threads = 10;
    my $pm = Parallel::ForkManager->new($max_threads);

    foreach my $x (@res_samples) {

        my $pid = $pm->start and next;

        my $tmp_dir   = qq{/home/tmp};
        my $gvcf_dir  = qq{$result/result/$x};
        system qq{mkdir -p $gvcf_dir} if not -d $gvcf_dir;

        my $bam_final = qq{$map/$x/$x\_final.bam};
        my $gVCF      = qq{$gvcf_dir/$x.g.vcf.gz};

        # 运行
        my $cmd     = qq{$Java -jar -Xmx6g $GATK4_Loc HaplotypeCaller -I $bam_final -O $gVCF -R $Genome -ERC GVCF --min-base-quality-score 20 --sample-ploidy 2 --tmp-dir $tmp_dir -L $RealignBed --max-reads-per-alignment-start 0 --dont-use-soft-clipped-bases true};
        my $finish  = qq{touch $gvcf_dir/$x.gvcf.finish};

        open SAVE, qq{>$result/run/$x.gvcf.sh} or die "Can't open $result/run/$x.gvcf.sh\n";
        print SAVE qq{$cmd\n};
        print SAVE qq{$finish\n};
        close SAVE;

        system qq{bash $result/run/$x.gvcf.sh &> $result/log/$x.gvcf.log};
        $pm->finish;

    }

    $pm->wait_all_children;
    print qq{gvcf分析完成!\n};

}

sub res_check
{
    my $output = shift;
    my $sample = shift;

    my @temp_sample = ();
    foreach my $x (@{$sample}) {
        next if -e qq{$output/$x/$x.gvcf.finish};
        push @temp_sample, $x;
    }

    return @temp_sample;
}

1