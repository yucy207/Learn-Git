package config;

use strict;
use warnings;

sub read_config
{
	my $config = shift;
	my %hash   = ();
	open FILE, $config or die "Can;t open $config,Please check it!\n";
	while (<FILE>) {

		s/[\r\n]//g;
		s/^\s+//g;
		s/\s+$//g;
		
		next if /^#/;
		next if /^\s*$/;
		my @arr = split /\s+/;
		if (scalar @arr == 2) {
			$hash{$arr[0]} = $arr[1];
			printf("%-15s", $arr[0]);
			print qq{: $arr[1]\n};
		}
		if ($_ =~/^group/) {
			my $num = @arr;
			for (my $i = 2; $i < $num; $i++) {
				printf("%-15s", $arr[0]);
				print qq{: $arr[1]/$arr[$i]\n};
				push @{$hash{'group'}}, qq{$arr[1]/$arr[$i]};
			}
		}
		if ($_ =~/^env/) {
			my $num = @arr;
			for (my $i = 2; $i < $num; $i++) {
				printf("%-15s", $arr[0]);
				print qq{: $arr[1]/$arr[$i]\n};
				push @{$hash{'env'}}, qq{$arr[1]/$arr[$i]};
			}
		}
		
	}
	close FILE;
	return \%hash;
}

sub base_config
{
	my %hash  = ();
	my $id = `id`;
	my ($uid, $gid) = $id =~ /uid=(\d+).+?gid=(\d+)/;
	my $docker_bin  = qq{docker run --rm -v /home:/home -u $uid:$gid};
	my $scriptdir   = qq{/home/genesky/pipeline/metagenome_16s_18s_its_sequencing/v2.3.1/script/};

	# soft
	$hash{'seq_crumbs_bin'}    = qq{$docker_bin chengsy_16s/seq_crumbs:v0.1.9};
	$hash{'trim_galore_bin'}   = qq{$docker_bin xudl/trim_galore:v0.4.5};
	$hash{'flash2_bin'}        = qq{$docker_bin chengsy_16s/flash2:v2.2.00};
	$hash{'fastx_toolkit_bin'} = qq{$docker_bin chengsy_16s/fastx_toolkit:v0.0.13};
	$hash{'mothur_bin'}        = qq{$docker_bin chengsy_16s/mothur:v.1.39.3};
	$hash{'alpha_mothur_bin'}  = qq{docker run --rm -v /home:/home chengsy_16s/mothur:v.1.39.3};
	$hash{'usearch_bin'}       = qq{$docker_bin chengsy_16s/usearch:v10};
	$hash{'qiime_bin'}         = qq{docker run --rm -v /home:/home yoshikiv/basespace-qiime-191-dev};
	$hash{'fasta_number_bin'}  = qq{$scriptdir/fasta_number.py};
	$hash{'muscle_bin'}        = qq{$docker_bin zhengy/muscle:v3.8.31};
	$hash{'FastTree_bin'}      = qq{$docker_bin zhengy/fasttree:v2.1.10};
	$hash{'kronatools_bin'}    = qq{$docker_bin zhengy/kronatools:v2.7};
	$hash{'Rscript_bin'}       = qq{$docker_bin zhengy/r-3.4.3:v1 Rscript};
	$hash{'picrust_bin'}       = qq{$docker_bin chengsy_16s/picrust};
	$hash{'lefse_bin'}         = qq{$docker_bin lefse};
	$hash{'vegan_bin'}         = qq{$docker_bin chengsy/vegan:v2.4-5};
	$hash{'rdptools'}          = qq{$docker_bin zhengy/rdptools:v-2.0.2-1};
	$hash{'blast'}             = qq{$docker_bin zhengy/ncbi-blast-2.7.1};
	$hash{'Rscript_3_5_3'}     = qq{$docker_bin shenlx_r-3.5.3:v1.4 Rscript};


	# dir
	$hash{'scriptdir'} = qq{$scriptdir};

	# database
	#primer
	$hash{'primerdb'} = qq{$scriptdir/primer.list};

	# RDP
	$hash{'RDP_16S'}{'ChimeraFA'}    = qq{/home/zhengy/bin/modules/database/SILVA/Chimera/silva_132_rep_set99_16S.fasta};
	$hash{'RDP_16S'}{'TaxonomyFA'}   = qq{/home/zhengy/bin/modules/database/RDP/RDP_11.5_Bacteria_Archaea.fasta};
	$hash{'RDP_16S'}{'Taxonomy'}     = qq{/home/zhengy/bin/modules/database/RDP/RDP_11.5_Bacteria_Archaea.tax};
	$hash{'RDP_16S'}{'Taxonomy_udb'} = qq{/home/zhengy/bin/modules/database/RDP/new_usearch_RDP_11.5_Bacteria_Archaea.fasta.udb};
	
	# SILVA
	$hash{'SILVA_16S'}{'ChimeraFA'}    = qq{/home/zhengy/bin/modules/database/SILVA/Chimera/silva_132_rep_set99_16S.fasta};
	$hash{'SILVA_16S'}{'TaxonomyFA'}   = qq{/home/zhengy/bin/modules/database/SILVA/SILVA_132_SSURef_Nr99_16S_trunc.fasta};
	$hash{'SILVA_16S'}{'Taxonomy'}     = qq{/home/zhengy/bin/modules/database/SILVA/SILVA_132_SSURef_Nr99_16S_trunc.tax};
	$hash{'SILVA_16S'}{'Taxonomy_udb'} = qq{/home/zhengy/bin/modules/database/SILVA/new_usearch_SILVA_132_SSURef_Nr99_16S_trunc.fasta.udb};
	$hash{'SILVA_18S'}{'ChimeraFA'}    = qq{/home/zhengy/bin/modules/database/SILVA/Chimera/silva_132_rep_set99_18S.fasta};
	$hash{'SILVA_18S'}{'TaxonomyFA'}   = qq{/home/zhengy/bin/modules/database/SILVA/SILVA_132_SSURef_Nr99_Eukaryota_trunc.fasta};
	$hash{'SILVA_18S'}{'Taxonomy'}     = qq{/home/zhengy/bin/modules/database/SILVA/SILVA_132_SSURef_Nr99_Eukaryota_trunc.tax};
	$hash{'SILVA_18S'}{'Taxonomy_udb'} = qq{/home/zhengy/bin/modules/database/SILVA/new_usearch_SILVA_132_SSURef_Nr99_Eukaryota_trunc.fasta.udb};
	
	# UNITE
	$hash{'UNITE_ITS1'}{'ChimeraFA'}   = qq{/home/zhengy/bin/modules/database/UNITE/uchime_reference_dataset_ITS1_28.06.2017.fasta};
	$hash{'UNITE_ITS2'}{'ChimeraFA'}   = qq{/home/zhengy/bin/modules/database/UNITE/uchime_reference_dataset_ITS2_28.06.2017.fasta};
	$hash{'UNITE_ITS'}{'TaxonomyFA'}   = qq{/home/zhengy/bin/modules/database/UNITE/UNITEv6_sh_99_new.fasta};
	$hash{'UNITE_ITS'}{'Taxonomy'}     = qq{/home/zhengy/bin/modules/database/UNITE/UNITEv6_sh_99_new.tax};
	$hash{'UNITE_ITS'}{'Taxonomy_udb'} = qq{/home/zhengy/bin/modules/database/UNITE/new_usearch_UNITEv6_sh_99_new.fasta.udb};
	
	# GREENGENE
	$hash{'GREENGENE_16S'}{'ChimeraFA'}  = qq{/home/panrf/database/SILVA/Chimera/Silva_128_rep_set99_16S.fasta};
	$hash{'GREENGENE_16S'}{'TaxonomyFA'} = qq{/home/panrf/database/Greengenes/99_otu_new.fasta};
	$hash{'GREENGENE_16S'}{'Taxonomy'}   = qq{/home/panrf/database/Greengenes/99_otu_taxonomy_new.tax};
	
	# PR2
	$hash{'PR2_18S'}{'ChimeraFA'}  = qq{/home/panrf/database/SILVA/Chimera/Silva_128_rep_set99_18S.fasta};
	$hash{'PR2_18S'}{'TaxonomyFA'} = qq{/home/panrf/database/PR2/pr2_gb203_version_4.5.fasta};
	$hash{'PR2_18S'}{'Taxonomy'}   = qq{/home/panrf/database/PR2/pr2_gb203_version_4.5.tax};

	# Greengenes for pi
	$hash{'Greengenes'}{'fasta'}    = qq{/home/panrf/database/Greengenes/gg_13_5_otus/rep_set/97_otus.fasta};
	$hash{'Greengenes'}{'Taxonomy'} = qq{/home/panrf/database/Greengenes/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt};
	
	# maarjAM
	$hash{'maarjAM'}{'fasta'}        = qq{/home/zhengy/bin/modules/database/maarjAM/maarjAM_2018.fasta};
	$hash{'maarjAM'}{'Taxonomy'}     = qq{/home/zhengy/bin/modules/database/maarjAM/maarjAM_2018.taxonomy.txt};
	$hash{'maarjAM'}{'Taxonomy_udb'} = qq{/home/zhengy/bin/modules/database/maarjAM/new_usearch_maarjAM_2018.fasta.udb};

	# FUNGENE
	# fungene database	
	$hash{'nifh'}{'database'} = qq{/home/panrf/database/FunGene/nifh/9.0/fungene_9.0_nifH_1989_unaligned_protein_seqs.4Analysis.fasta};
	$hash{'nifh'}{'blastdb'}  = qq{/home/panrf/database/FunGene/nifh/9.0/fungene_9.0_nifH_1989};
		
	$hash{'nirS'}{'database'} = qq{/home/panrf/database/FunGene/nirS/9.0/fungene_9.0_nirS_1124_unaligned_protein_seqs.4Analysis.fasta};
	$hash{'nirS'}{'blastdb'}  = qq{/home/panrf/database/FunGene/nirS/9.0/fungene_9.0_nirS_1124};
		
	$hash{'nosZ'}{'database'} = qq{/home/panrf/database/FunGene/nosZ/9.0/fungene_9.0_nosZ_1655_unaligned_protein_seqs.4Analysis.fasta};
	$hash{'nosZ'}{'blastdb'}  = qq{/home/panrf/database/FunGene/nosZ/9.0/fungene_9.0_nosZ_1655};
		
	$hash{'nirK'}{'database'} = qq{/home/panrf/database/FunGene/nirK/9.0/fungene_9.0_nirK_442_unaligned_protein_seqs.4Analysis.fasta};
	$hash{'nirK'}{'blastdb'}  = qq{/home/panrf/database/FunGene/nirK/9.0/fungene_9.0_nirK_442};
		
	$hash{'narG'}{'database'} = qq{/home/panrf/database/FunGene/narG/9.0/fungene_9.0_narG_32810_unaligned_protein_seqs.4Analysis.fasta};
	$hash{'narG'}{'blastdb'}  = qq{/home/panrf/database/FunGene/narG/9.0/fungene_9.0_narG_32810};
		
	$hash{'norB'}{'database'} = qq{/home/panrf/database/FunGene/norB/9.1/fungene_9.1_norB_6004_unaligned_protein_seqs.4Analysis.fasta};
	$hash{'norB'}{'blastdb'}  = qq{/home/panrf/database/FunGene/norB/9.1/fungene_9.1_norB_6004"};
		
	$hash{'phoD'}{'database'} = qq{/home/chengsy/data/fungene/phoD/9.4/fungene_9.4_phoD_3358_unaligned_protein_seqs.4Analysis.fasta};
	$hash{'phoD'}{'blastdb'}  = qq{/home/chengsy/data/fungene/phoD/9.4/fungene_9.4_phoD_3358};
	
	return \%hash;
}
1;
