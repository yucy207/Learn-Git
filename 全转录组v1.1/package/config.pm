package package::config;

sub read_config
{
	my $config_txt = shift;
	my %hash   = ();
	open EXO, $config_txt or die "Can;t open $config_txt\n";
	while (<EXO>) {
		chomp;
		s/^\s+//;
		s/\s+$//;
		next if /^#/;
		next if /^\s*$/;

		my @arr = split /\s+/;		
		if (scalar @arr == 2 && $arr[0] ne "time") {
			$hash{$arr[0]} = $arr[1];
			printf("%-10s", $arr[0]);
			print qq{: $arr[1]\n};
		} 
		if (/^group/) {
			push @{$hash{'group'}}, [$arr[1], $arr[2], $arr[3]]; 
			next;
		}
		if (/^overlap/) {
			push @{$hash{'overlap'}}, [$arr[1], $arr[2]];
			next;
		}
		if (/^anova/) {
			push @{$hash{'anova'}}, [$arr[1], $arr[2]];
			next;
		}
		if (/^time/) {
			push @{$hash{'time'}}, $arr[1];
			next;
		}

	}
	close EXO;
	my $sample_num = split /,/, $hash{'samples'};
	printf("%-10s", "Sample Num");
	print qq{: $sample_num\n};
	return \%hash;

}


sub base_config
{
	my $id = `id`;
	my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;

	my %hash  = ();
	$hash{'util'}         = qq{/home/genesky/pipeline/whole_transriptome_sequencing/v1.1/util};
	$hash{'rscript_bin'}  = qq{docker run --rm -v /home:/home -u $uid:$gid  xudl/r-3.4.4:v1.0};
	
	######################################## default options ##########################################
	$hash{'thread_qc'}    = 5;
	$hash{'thread_map'}   = 4;
	$hash{'max_threads'}  = 5;
	$hash{'log2fc'}       = 1;
	$hash{'log2fc_circ'}  = 0;

	########################################## softwares ##############################################
	# QC
	$hash{'fastqc_bin'}    = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/fastqc:v0.11.7};
	$hash{'fastx_bin'}     = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/fastx_toolkit:v0.0.13};
	$hash{'seqtk_bin'}     = qq{/home/genesky/software/seqtk/1.3-r106/seqtk};
	$hash{'fastp_bin'}     = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/fastp:0.15.0};
	$hash{'adapter_1'}{'nextera'} = qq{CTGTCTCTTATACACATCTCCGAGCCCACGAGAC};
	$hash{'adapter_2'}{'nextera'} = qq{CTGTCTCTTATACACATCTGACGCTGCCGACGA};
	$hash{'adapter_1'}{'truseq'}  = qq{AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC};
	$hash{'adapter_2'}{'truseq'}  = qq{AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT};
	# mapping
	$hash{'hisat2_bin'}    = qq{/home/genesky/software/hisat2/2.1.0/hisat2};
	$hash{'samtools_bin'}  = qq{/home/genesky/software/samtools/1.9/samtools};
	$hash{'sambamba_bin'}  = qq{/home/genesky/software/sambamba/0.6.7/sambamba};
	$hash{'rseqc_bin'}     = qq{docker run --rm -v /home:/home xudl/rseqc:latest};
	$hash{'rnaseqc_bin'}   = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/rna-seqc:v1.1.8};
	# 定量
	$hash{'stringtie_bin'} = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/stringtie:v1.3.4a};
	$hash{'ballgown_bin'}  = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/ballgown:v2.10.0};
	$hash{'pca_bin'}       = qq{docker run --rm -v /home:/home -u $uid:$gid transcriptome_pca:v1.0};
	$hash{'Rtsne_bin'}     = qq{docker run --rm -v /home:/home -u $uid:$gid wangly/r-3.4.4:v1.0};
	$hash{'cor_bin'}       = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/r-3.4.4:v1.2};
	# 富集分析
	$hash{'clusterprofiler_bin'}   = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/clusterprofiler:v3.6};
	$hash{'venndiagram_bin'}       = qq{docker run --rm -v /home:/home xudl/venndiagram:v1.6.20};
	# 可变剪切
	$hash{'rmats4_bin'}            = qq{docker run --rm -v /home:/home xudl/rmats4:latest};
	$hash{'rmats2sashimiplot_bin'} = qq{docker run --rm -v /home:/home xudl/rmats2sashimiplot:v2.0.3};
	# 转录本分析
	$hash{'cufflinks_bin'}         = qq{docker run --rm -v /home:/home xudl/cufflinks:v2.1.1};
	$hash{'gtf_to_fasta_bin'}      = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/tophat:v2.1.0};
	$hash{'transdedcoder_bin'}     = qq{/home/genesky/software/transdedcoder/3.0.0/TransDecoder.LongOrfs};
	$hash{'diamond_bin'}           = qq{/home/xudl/soft/diamond/diamond};
	# 融合基因
	$hash{'mono_bin'}              = qq{docker run --rm -v /home:/home -u $uid:$gid chengsy/mono:latest};
	$hash{'fusionmap_bin'}         = qq{/home/genesky/software/oshell/oshell.exe};
	# 蛋白交互(ppi)
	$hash{'stringdb_bin'}          = qq{docker run --rm -v /home:/home -u $uid:$gid chengsy/stringdb:v10};
	# wgcna
	$hash{'wgcna_bin'}             = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/wgcna:v1.64};
	# lncRNA
	$hash{'lncRNA_gtf_to_fasta_bin'} = qq{/home/xudl/soft/tophat-2.1.0.Linux_x86_64/gtf_to_fasta};
	# circRNA
	$hash{'bwa_bin'}               = qq{/home/genesky/software/bwa/0.7.17/bwa};
	$hash{'ciri_bin'}              = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/ciri:v2.0.6};
	$hash{'iresfinder_bin'}        = qq{docker run --rm -v /home:/home chengsy/iresfinder};
	$hash{'cpat_bin'}              = qq{docker run --rm -v /home:/home chengsy/cpat};
	$hash{'miranda_bin'}           = qq{docker run --rm -v /home:/home -u $uid:$gid xudl/miranda:v3.3a};
	# 时间序列
	$hash{'masigpro'}              = qq{docker run --rm -v /home:/home -u $uid:$gid dongxj/masigpro:v2};
	#SNP
	$hash{'sambamba'}              = qq{/home/genesky/software/sambamba/0.6.7/sambamba};
	$hash{'cpulimit'}              = qq{/home/genesky/software/cpulimit/0.2/cpulimit};
	$hash{'Java'}                  = qq{/home/genesky/software/java/1.8.0_181/bin/java};
	$hash{'Picard'}                = qq{/home/genesky/software/picard/2.18.29/picard.jar};
	$hash{'GATK3'}                 = qq{/home/genesky/software/gatk/3.5/GenomeAnalysisTK.jar};
	$hash{'GATK4_Loc'}             = qq{/home/genesky/software/gatk/4.0.12.0/gatk-package-4.0.12.0-local.jar};
	$hash{'bcftools'}              = qq{/home/genesky/software/bcftools/1.9/bin/bcftools};
	$hash{'AnnovarDir'}            = qq{/home/genesky/software/annovar/2018Apr16};
	$hash{'MegaBlast'}             = qq{/home/genesky/software/blast/20110130/bin/megablast};
	$hash{'Python27'}              = qq{/home/genesky/software/python/2.7.13/bin/python};
	$hash{'InterVar'}              = qq{/home/genesky/software/intervar/2.1.2};
	$hash{'snpEFF'}                = qq{/home/genesky/software/snpeff/4_3t/snpEff/snpEff.jar};
	

	####################################### database(species) ############################################
	$hash{'nr_db'}                 = qq{/home/pub/database/NT_NR/2017_6_13/nr};
	$hash{'gi_to_protein'}         = qq{/home/genesky/database/self_build_database/WTS/protein_acc_2_annotation.txt};

	# human(hg19)
	$hash{'hg19'}{'three_letter'}         = qq{hsa};
	$hash{'hg19'}{'official_name'}        = qq{Homo_sapiens};
	$hash{'hg19'}{'common_name'}          = qq{human};
	$hash{'hg19'}{'hisat2_index'}         = qq{/home/genesky/database/self_build_database/WTS/hg19/hisat2_trans_db/hg19_trans};
	$hash{'hg19'}{'genome_bed12'}         = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_refGene_bed12.txt};
	$hash{'hg19'}{'genome_chr_length'}    = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19.chr.size};
	$hash{'hg19'}{'genome_fasta'}         = qq{/home/pub/database/Human/hg19/bowtie2_db/hg19.fa};
	$hash{'hg19'}{'genome_gencode'}       = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_gencode.gtf};
	$hash{'hg19'}{'genome_mRNA_gtf'}      = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_mRNA_refGene.gtf};
	$hash{'hg19'}{'go_info_db'}           = qq{/home/genesky/database/self_build_database/WTS/hg19/go_info};
	$hash{'hg19'}{'kegg_info_db'}         = qq{/home/genesky/database/self_build_database/WTS/hg19/kegg.annotation.txt};
	$hash{'hg19'}{'pathway_position_db'}  = qq{/home/genesky/database/self_build_database/WTS/hg19/hsa.gene.position.on.pathway.txt};
	$hash{'hg19'}{'genome_gtf'}           = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_mRNA_lncRNA.gtf};
	$hash{'hg19'}{'fusionmap_db'}         = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_fusionmap_db};
	$hash{'hg19'}{'fusionmap_db_release'} = qq{hg19_lib};
	$hash{'hg19'}{'fusionmap_db_model'}   = qq{hg19_model};

	$hash{'hg19'}{'genome_lncRNA_gtf'}    = qq{/home/genesky/database/self_build_database/WTS/hg19/lncipedia_5_2_hg19.lncRNA.gtf};
	$hash{'hg19'}{'genome_lncRNA_one_gtf'}= qq{/home/genesky/database/self_build_database/WTS/hg19/lncipedia_5_2_hg19.lncRNA.one.trans.gtf};
	$hash{'hg19'}{'genome_mRNA_one_gtf'}  = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_mRNA_one_trans.gtf};

	$hash{'hg19'}{'circbase_db'}          = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_circRNA.txt};
	$hash{'hg19'}{'circ_ref_gtf'}         = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_mRNA_refGene.gtf};
	$hash{'hg19'}{'bwa_index'}            = qq{/home/xudl/database/hg19/bwa_db/hg19.fa};
	$hash{'hg19'}{'miRNA_fasta'}          = qq{/home/genesky/database/self_build_database/WTS/hg19/hsa.mature.dna.fa};

	$hash{'hg19'}{'mRNA_miRNA_db'}        = qq{/home/genesky/database/self_build_database/WTS/hg19/mRNA_miRNA_interaction.xls};
	$hash{'hg19'}{'lncRNA_miRNA_db'}      = qq{/home/genesky/database/self_build_database/WTS/hg19/lncRNA_miRNA_interaction.xls};

	$hash{'hg19'}{'InDel'}                = qq{/home/genesky/database/gatk/hg19/1000G_phase1_indels/1000G_phase1.indels.hg19.sites.vcf.gz};  # 1000g indel 用于GATK校正
    $hash{'hg19'}{'InDelGold'}            = qq{/home/genesky/database/gatk/hg19/mills_and_1000g_gold_standard_indels/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz}; # 1000g indel 用于GATK校正
    $hash{'hg19'}{'DBsnp'}                = qq{/home/genesky/database/gatk/hg19/dbsnp/dbsnp_138.hg19.vcf.gz};  # DBSNP 用于GATK校正
    $hash{'hg19'}{'RealignBed'}           = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_ref.bed};
    $hash{'hg19'}{'TargetBed'}            = qq{/home/genesky/database/self_build_database/WTS/hg19/hg19_target.bed};
    $hash{'hg19'}{'AnnovarBuild'}         = qq{hg19};
    $hash{'hg19'}{'blast_idx'}            = qq{/home/genesky/database/ucsc/hg19/genome/blast+_idx/hg19.fa};
    $hash{'hg19'}{'snpEFF_config'}        = qq{/home/genesky/database/snpeff/4_3t/snpEff.config};
    $hash{'hg19'}{'snpEFF_build'}         = qq{hg19};

	# mouse(mm10)
	$hash{'mm10'}{'three_letter'}         = qq{mmu};
	$hash{'mm10'}{'official_name'}        = qq{Mus_musculus};
	$hash{'mm10'}{'common_name'}          = qq{mouse};
	$hash{'mm10'}{'hisat2_index'}         = qq{/home/genesky/database/self_build_database/WTS/mm10/hisat2_trans_db/mm10_trans};
	$hash{'mm10'}{'genome_bed12'}         = qq{/home/genesky/database/self_build_database/WTS/mm10/mm10_refGene_bed12.txt};
	$hash{'mm10'}{'genome_chr_length'}    = qq{/home/genesky/database/self_build_database/WTS/mm10/mm10.chr.size};
	$hash{'mm10'}{'genome_fasta'}         = qq{/home/xudl/database/mm10/hisat2_db/mm10.fa};
	$hash{'mm10'}{'genome_gencode'}       = qq{/home/genesky/database/self_build_database/WTS/mm10/mm10_gencode_nochr.gtf};
	$hash{'mm10'}{'genome_mRNA_gtf'}      = qq{/home/genesky/database/self_build_database/WTS/mm10/mm10_mRNA_refGene.gtf};
	$hash{'mm10'}{'go_info_db'}           = qq{/home/genesky/database/self_build_database/WTS/mm10/go_info};
	$hash{'mm10'}{'kegg_info_db'}         = qq{/home/genesky/database/self_build_database/WTS/mm10/kegg.annotation.txt};
	$hash{'mm10'}{'pathway_position_db'}  = qq{/home/genesky/database/self_build_database/WTS/mm10/mmu.gene.position.on.pathway.txt};
	$hash{'mm10'}{'genome_gtf'}           = qq{/home/genesky/database/self_build_database/WTS/mm10/mm10_mRNA_refGene.gtf};
	$hash{'mm10'}{'fusionmap_db'}         = qq{/home/genesky/database/self_build_database/WTS/mm10/mm10_fusionmap_db};
	$hash{'mm10'}{'fusionmap_db_release'} = qq{Mouse.mm10};
	$hash{'mm10'}{'fusionmap_db_model'}   = qq{RefGene};

	$hash{'mm10'}{'genome_lncRNA_gtf'}    = qq{/home/genesky/database/self_build_database/WTS/mm10/NONCODE2016_mouse_mm10_lncRNA.gtf};
	$hash{'mm10'}{'genome_lncRNA_one_gtf'}= qq{/home/genesky/database/self_build_database/WTS/mm10/NONCODE2016_mouse_mm10_lncRNA_one_trans.gtf};
	$hash{'mm10'}{'genome_mRNA_one_gtf'}  = qq{/home/genesky/database/self_build_database/WTS/mm10/mm10_mRNA_one_trans.gtf};

	$hash{'mm10'}{'circbase_db'}          = qq{/home/genesky/database/self_build_database/WTS/mm10/mm9_circRNA.txt};
	$hash{'mm10'}{'circ_ref_gtf'}         = qq{/home/genesky/database/self_build_database/WTS/mm10/parse_circRNA_mm10.gtf};
	$hash{'mm10'}{'bwa_index'}            = qq{/home/xudl/database/mm10/hisat2_db/mm10.fa};
	$hash{'mm10'}{'miRNA_fasta'}          = qq{/home/genesky/database/self_build_database/WTS/mm10/mmu.mature.dna.fa};

	$hash{'mm10'}{'mRNA_miRNA_db'}        = qq{/home/genesky/database/self_build_database/WTS/mm10/mRNA_miRNA_interaction.xls};
	$hash{'mm10'}{'lncRNA_miRNA_db'}      = qq{/home/genesky/database/self_build_database/WTS/mm10/lncRNA_miRNA_interaction.xls};

	# rat(rnor6)
	$hash{'rnor6'}{'three_letter'}         = qq{rno};
	$hash{'rnor6'}{'official_name'}        = qq{Rattus_norvegicus};
	$hash{'rnor6'}{'common_name'}          = qq{rat};
	$hash{'rnor6'}{'hisat2_index'}         = qq{/home/genesky/database/self_build_database/WTS/rnor6/hisat2_trans_db/rnor6_trans};
	$hash{'rnor6'}{'genome_bed12'}         = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_refGene_bed12.txt};
	$hash{'rnor6'}{'genome_chr_length'}    = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6.chr.size};
	$hash{'rnor6'}{'genome_fasta'}         = qq{/home/pub/database/Rnor6/Rnor6.fa};
	$hash{'rnor6'}{'genome_gencode'}       = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_gencode.gtf};
	$hash{'rnor6'}{'genome_mRNA_gtf'}      = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_mRNA_refGene.gtf};
	$hash{'rnor6'}{'go_info_db'}           = qq{/home/genesky/database/self_build_database/WTS/rnor6/go_info};
	$hash{'rnor6'}{'kegg_info_db'}         = qq{/home/genesky/database/self_build_database/WTS/rnor6/kegg.annotation.txt};
	$hash{'rnor6'}{'pathway_position_db'}  = qq{/home/genesky/database/self_build_database/WTS/rnor6/rno.gene.position.on.pathway.txt};
	$hash{'rnor6'}{'genome_gtf'}           = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_mRNA_lncRNA.gtf};
	$hash{'rnor6'}{'fusionmap_db'}         = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_fusionmap_db};
	$hash{'rnor6'}{'fusionmap_db_release'} = qq{Rnor6};
	$hash{'rnor6'}{'fusionmap_db_model'}   = qq{Rnor6};

	$hash{'rnor6'}{'genome_lncRNA_gtf'}    = qq{/home/genesky/database/self_build_database/WTS/rnor6/NONCODEv5_rat_rn6_lncRNA.gtf};
	$hash{'rnor6'}{'genome_lncRNA_one_gtf'}= qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_lncRNA_one_trans.gtf};
	$hash{'rnor6'}{'genome_mRNA_one_gtf'}  = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_mRNA_one_trans.gtf};

	$hash{'rnor6'}{'circ_ref_gtf'}         = qq{/home/genesky/database/self_build_database/WTS/rnor6/rnor6_mRNA_refGene.gtf};
	$hash{'rnor6'}{'bwa_index'}            = qq{/home/pub/database/Rnor6/Rnor6.fa};
	$hash{'rnor6'}{'miRNA_fasta'}          = qq{/home/genesky/database/self_build_database/WTS/rnor6/rno.mature.dna.fa};

	$hash{'rnor6'}{'mRNA_miRNA_db'}        = qq{/home/genesky/database/self_build_database/WTS/rnor6/mRNA_miRNA_interaction.xls};
	$hash{'rnor6'}{'lncRNA_miRNA_db'}      = qq{/home/genesky/database/self_build_database/WTS/rnor6/lncRNA_miRNA_interaction.xls};


	return \%hash;
	
}

sub error
{
	my $message = qq{
        (()__(()
        /       \ 
       ( /    \  \
        \ o o    /
        (_()_)__/ \
       / _,==.____ \
      (   |--|      )
      /\_.|__|'-.__/\_
     / (        /     \ 
     \  \      (      /
      )  '._____)    /
   (((____.--(((____/
	};
	
	print qq{$message};
}

1;
