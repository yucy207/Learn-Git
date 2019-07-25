Sheet Name	标注与说明	
染色体	Genetic Models	按照不同遗传模式进行突变位点筛选，Recessive(ChrX)：X染色体上的隐性遗传；Recessive(AutoHom)：同一来源的常染色体隐性遗传； Recessive(AutoComHet)：不同来源的常染色体隐性遗传； IncompleteDominance：不完全显性遗传；De Novo Mutation：新突变
染色体	SNV NO	Single Nucleotide Variant (SNV, 包括插入或缺失）编号
染色体	First Priority	最优先考虑的位点，First1：^1） 在HGMD中已有突变位点,或者以下^2） 保守^3） 在1000genome中频率低于0.001或者在自身control中小于0.01^4） 在ESP6500中频率小于0.01^5） SNP Calling Quality不为L 而且分型质量<20的样本比例不大于50%^6） 同源性为1^7） 在genesky Database中的突变频率小于0.005^First2：^1） 在1000genome中频率低于0.001或者在自身control中小于0.01^2） 在ESP6500中频率小于0.01^3） SNP Calling Quality不为L 而且分型质量<20的样本比例不大于50%^4） 同源性为1^5） 在genesky Database中的突变频率小于0.005^Second^1） SNP Calling Quality不为L 而且分型质量<20的样本比例不大于50%^2）  在1000genome中频率低于0.01或者在自身control中小于0.01^3） 同源性<3^Third:^其他剩余位点
染色体	SNP Calling Quality	我们采用GATK方法进行SNPcalling，如果该位点被calling到，则为H，若该位点为STR/附近50bp内存在STR/附近3bp内存在INDEL,则降一级为M，如果位点在vcf中“FILTER”一列不为“PASS”或位点两侧存在同源序列，再降一级
染色体	Genotyping Quality	突变位点整体分型质量。为项目中，该位点分型质量>20的样本比例
染色体	SNP ID	位点在ncbi的dbsnp中的编号。如果编号上标注“STR”表示位点很有可能是STR，如果编号上标注“STR_around”或“INDEL_around”表示位点周围存在STR或INDEL
染色体	Freq_Alt (1000)	非参照等位基因在千人基因组项目中的频率，参考“http://www.1000genomes.org/home”
染色体	1000g_chbs	中国北方南方样本突变频率
染色体	Ref Allele	位点在参照序列上的等位基因
染色体	Alt Allele	位点的另一个等位基因
染色体	Ref_dp	参照等位基因的测序深度
染色体	Alt_dp	另一等位基因的测序深度
染色体	Chrs	染色体号
染色体	Position	在对应染色体上的位置，以Reference GRCh37.p5 Primary Assembly为准
染色体	Gene Strand Orientation	显示基因在染色体上的排列方向，+表明基因序列同参考序列一致，-表明与参考序列互补
染色体	Gene Region	变异所在的基因区域（如exonic=coding；splicing=2bp aroudn splicing junction；ncRNA；UTR5=5'untranslated region；UTR3=3'untranslated region；intronic；upstream =1kb upstream of transcription strart site；downstream=1kb downstream of transcription end site；intergenic)
染色体	 Function	基因功能（包括：nonsynonymous SNV，frameshift insertion，nonframeshift insertion,Stop Gain,Stop Loss）
染色体	 tfbsConsSites Score	保守的转录因子结合区域 来源参考 http://ucscbrowser.genap.ca/cgi-bin/hgTables?hgsid=213993_qIcXI6SyNhqV1qBESrQyMiXkifvT&hgta_doSchemaDb=hg19&hgta_doSchemaTable=tfbsConsSites
染色体	SIFT Score	SIFT (Sorting Intolerant From Tolerant) 是一个用于预测氨基酸改变是否影响功能的程序，参考http://sift.jcvi.org/www/SIFT_chr_coords_submit.html；Nucleic Acids Research, 2003, Vol. 31, No. 13；通常认为小于0.05的SIFT数值表明该变异会严重影响蛋白功能的。
染色体	POLYPHEN  Score	Polymorphism Phenotyping 类似于SIFT，是另外一种用来预测非同义突变造成的氨基酸改变是否影响功能的算法，参考http://genetics.bwh.harvard.edu/pph/pph_help_text.html
染色体	MutationTaster  Score	位点时突变的可能性分值，分值越高越显著 参考：http://sites.google.com/site/jpopgen/dbNSFP. dbNSFP: a lightweight database of human nonsynonymous SNPs and their functional predictions.Hum Mutat. 2011 Aug;32(8):894-9. doi: 10.1002/humu.21517.
染色体	SIFT Score Pred	SIFT官方注释,D=DAMAGING,T=TOLERATED
染色体	POLYPHEN  Score Pred	POLYPHEN官方注释,B=benign,P=possibly damaging,D=probably damaging
染色体	MutationTaster  Score Pred	MutationTaster官方注释,A=disease_causing_automatic,D=disease_causing,N=polymorphism,P=polymorphism_automatic
染色体	Cadd_Raw	SNV危险性评分，with higher values indicating that a variant is more likely to be simulated (or "not observed") and therefore more likely to have deleterious effects。参考http://cadd.gs.washington.edu.Cutoff is usually set 4.
染色体	Cadd_Phred	CADD phred-scaled scores;Basically, for scaled scores, 10 means 10% percentile highest scores, 20 means 1% percentile highest scores, and 30% means 0.1% percentile highest scores. So we do find one <1% percentile variant in the dataset(Cadd_Raw=~4)
染色体	Dann	SNV危险性评分，DANN scores whole-genome variants by training a deep neural network (DNN). DNNs can capture non-linear relationships among features and are better suited than SVMs for problems with a large number of samples and features. It seems to perform much better than CADD。参考http://bioinformatics.oxfordjournals.org/content/31/5/761.Cutoff is usually set 0.93.
染色体	Eigen	SNV危险性评分，Eigen?score uses a spectral approach integrating functional genomic annotations for coding and noncoding variants. Unlike other approaches such as CADD/DANN, Eignen did not use labelled training data and should be considered as a unsupervised learning approach.参考：http://www.ncbi.nlm.nih.gov/pubmed/26727659?
染色体	dbscSNV_ADA_SCORE	基于AdaBoost 算法计算的突变对剪切位点产生影响的可能性评分。值越大，概率越大，常用阈值>0.6
染色体	dbscSNV_RF_SCORE	基于Random Forest 算法计算的突变对剪切位点产生影响的可能性评分。值越大，概率越大，常用阈值>0.6
染色体	Hrcr1	全基因组突变频率，Haplotype Reference Consortium database with 40 million variants from 32K samples in haplotype reference consortium
染色体	Kaviar	全基因组突变频率，Kaviar database with 170 million variants from 13K genomes and 64K exomes。参考http://db.systemsbiology.net/kaviar/
染色体	ESP6500	来源于ESP( NHLBI GO Exome Sequencing Project )全外显子测序的6500个样本中的频率，请参考:http://evs.gs.washington.edu/EVS/
染色体	ExAC03	65000 exome allele frequency data for ALL. version 0.3. Left normalization done
染色体	ExAC03_EAS	65000 exome allele frequency data for EAS (East Asian). version 0.3. Left normalization done
染色体	gnomAD_exome	The Genome Aggregation Database (gnomAD) is a resource developed by an international coalition of investigators.The data set spans 123,136 exome sequences and 15,496 whole-genome sequences from unrelated individuals sequenced as part of various disease-specific and population genetic studies. 
染色体	gnomAD_exome_EAS	gnomAD data for EAS (East Asian).
染色体	gnomAD_genome	The Genome Aggregation Database (gnomAD) is a resource developed by an international coalition of investigators.The data set spans 123,136 exome sequences and 15,496 whole-genome sequences from unrelated individuals sequenced as part of various disease-specific and population genetic studies. 
染色体	gnomAD_genome_EAS	gnomAD data for EAS (East Asian).
染色体	Predicted Protein Variants	根据参照基因序列预测的基因突变引起的蛋白的氨基酸的变化
染色体	Interpro_domain	protein domain for variants
染色体	OMIM	基因在OMIM数据库中的信息，格式为染色体区段|遗传模式|OMIM编号:疾病名称
染色体	HPO	The Human Phenotype Ontology (http://www.human-phenotype-ontology.org) aims to provide a standardized vocabulary of phenotypic abnormalities encountered in human disease. 
染色体	HGMD_site	位点在HGMD数据库中的疾病/表型|pubmed编号
染色体	HGMD_site_class	致病性分类。DM:Disease-causing mutations;DM?:denoting a probable/possible pathological mutation;DP:Disease-associated polymorphisms;FP:Functional polymorphisms;DFP:Disease-associated polymorphisms with supporting functional evidence;R:Retired record;
染色体	HGMD_gene	基因在HGMD数据库中的疾病/表型，格式为Phenotype1;Phenotype2;Phenotype3
染色体	InterVar	InterVar根据ACMG评分规则进行的突变临床解释信息
染色体	InterVar_evidence	InterVar根据ACMG评分规则进行的28项打分结果
染色体	Predicted Protein Variants(snpEFF)	使用snpEFF进行突变注释，注释信息为：Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
染色体	MGI	MGI(http://www.informatics.jax.org/) is the international database resource for the laboratory mouse, providing integrated genetic, genomic, and biological data to facilitate the study of human health and disease.
染色体	ClinVar	NCBI疾病突变数据库,注释内容分别为CLINSIG,CLNDBN,CLNACC,CLNDSDB,CLNDSDBID
染色体	MalaCards	MalaCards(http://www.malacards.org/) is an integrated database of human maladies and their annotations, modeled on the architecture and richness of the popular GeneCards database of human genes.
染色体	COSMIC	Catalogue Of Somatic Mutations In Cancer. It includes somatic mutations reported in literature in various types of cancers
染色体	NCI60	The NCI-60 cell lines are the most frequently studied human tumor cell lines in cancer research.
染色体	ICGC	International Cancer Genome Consortium.该位点在国际癌症基因组ICGC（http://www.icgc.org/）数据库中的信息.The annotations include a ICGC_ID column and a ICGC_Occurrence column. The ICGC_Occurrence column includes the project in which the mutation was identified, the number of donors affected, the total number of donors studied in the project and the frequency of the mutation, separated by ""|"".
染色体	VCF Filter	VCF过滤硬条件描述
染色体	GO_BP	Gene Ontology annotation BP : biological process
染色体	GO_MF	Gene Ontology annotation MF : Molecular Function
染色体	GO_CC	Gene Ontology annotation CC: Cellular Component
染色体	KEGG_Pathway	该基因在kegg的pathway数据库中的信息
染色体	5'Flanking Sequence	SNV 5’侧的50bp碱基
染色体	3'Flanking Sequence	SNV 3’侧的50bp碱基
染色体	Homology Hits	包含该位点在内的左右各50bp组成的序列在人基因组中的同源序列数目，同源性达到90%以上就算HIT到一次，HIT次数越多该位点可靠性越低。
染色体	HWE	哈迪-温伯格平衡(Hardy-Weinberg equilibrium,HWE)的计算值
染色体	GeneskyExonDB_Freq	天昊300样本外显子数据库频率
染色体	GeneskyGenomeDB_Freq	天昊128样本全基因组数据库频率
染色体	Case(0|1|2)	正常纯合(0)、杂合(1)、突变纯合(2)的Case样本数量
染色体	Control(0|1|2)	正常纯合(0)、杂合(1)、突变纯合(2)的Control样本数量
染色体	Control Percentage	Alt Allele在control样本中的比例
染色体	Alt Allele Freq	Alt Allele在所有样本中的比例
