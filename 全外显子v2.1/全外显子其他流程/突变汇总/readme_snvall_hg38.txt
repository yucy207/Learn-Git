Sheet Name	标注与说明	
SNV Information	Genetic Models	按照不同遗传模式进行突变位点筛选，Recessive(ChrX)：X染色体上的隐性遗传；Recessive(AutoHom)：同一来源的常染色体隐性遗传； Recessive(AutoComHet)：不同来源的常染色体隐性遗传； IncompleteDominance：不完全显性遗传；De Novo Mutation：新突变
SNV Information	SNV NO	Single Nucleotide Variant (SNV, 包括插入或缺失）编号
SNV Information	First Priority	最优先考虑的位点，First1：^1） 在HGMD中已有突变位点,或者以下^2） 保守^3） 在1000genome中频率低于0.001或者在自身control中小于0.01^4） 在ESP6500中频率小于0.01^5） SNP Calling Quality不为L 而且分型质量<20的样本比例不大于50%^6） 同源性为1^7） 在genesky Database中的突变频率小于0.005^First2：^1） 在1000genome中频率低于0.001或者在自身control中小于0.01^2） 在ESP6500中频率小于0.01^3） SNP Calling Quality不为L 而且分型质量<20的样本比例不大于50%^4） 同源性为1^5） 在genesky Database中的突变频率小于0.005^Second^1） SNP Calling Quality不为L 而且分型质量<20的样本比例不大于50%^2）  在1000genome中频率低于0.01或者在自身control中小于0.01^3） 同源性<3^Third:^其他剩余位点
SNV Information	SNP Calling Quality	我们采用GATK方法进行SNPcalling，如果该位点被calling到，则为H，若该位点为STR/附近50bp内存在STR/附近3bp内存在INDEL,则降一级为M，如果位点在vcf中“FILTER”一列不为“PASS”或位点两侧存在同源序列，再降一级
SNV Information	Genotyping Quality	突变位点整体分型质量。为项目中，该位点分型质量>20的样本比例
SNV Information	SNP ID	位点在ncbi的dbsnp中的编号。如果编号上标注“STR”表示位点很有可能是STR，如果编号上标注“STR_around”或“INDEL_around”表示位点周围存在STR或INDEL
SNV Information	Freq_Alt (1000)	非参照等位基因在千人基因组项目中的频率，参考“http://www.1000genomes.org/home”
SNV Information	1000g_eas	alternative allele frequency data in 1000 Genomes Project for EAS (East Asian)
SNV Information	Ref Allele	位点在参照序列上的等位基因
SNV Information	Alt Allele	位点的另一个等位基因
SNV Information	Ref_dp	参照等位基因的测序深度
SNV Information	Alt_dp	另一等位基因的测序深度
SNV Information	Chrs	染色体号
SNV Information	Position	在对应染色体上的位置，以Reference GRCh38.p12 Primary Assembly为准
SNV Information	Gene Strand Orientation	显示基因在染色体上的排列方向，+表明基因序列同参考序列一致，-表明与参考序列互补
SNV Information	Gene Region	变异所在的基因区域（如exonic=coding；splicing=2bp aroudn splicing junction；ncRNA；UTR5=5'untranslated region；UTR3=3'untranslated region；intronic；upstream =1kb upstream of transcription strart site；downstream=1kb downstream of transcription end site；intergenic)
SNV Information	Function	基因功能（包括：nonsynonymous SNV，frameshift insertion，nonframeshift insertion,Stop Gain,Stop Loss）
SNV Information	SIFT Score	SIFT (Sorting Intolerant From Tolerant) 是一个用于预测氨基酸改变是否影响功能的程序，参考http://sift.jcvi.org/www/SIFT_chr_coords_submit.html；Nucleic Acids Research, 2003, Vol. 31, No. 13；通常认为小于0.05的SIFT数值表明该变异会严重影响蛋白功能的。
SNV Information	POLYPHEN  Score	Polymorphism Phenotyping 类似于SIFT，是另外一种用来预测非同义突变造成的氨基酸改变是否影响功能的算法，参考http://genetics.bwh.harvard.edu/pph/pph_help_text.html
SNV Information	MutationTaster  Score	位点时突变的可能性分值，分值越高越显著 参考：http://sites.google.com/site/jpopgen/dbNSFP. dbNSFP: a lightweight database of human nonsynonymous SNPs and their functional predictions.Hum Mutat. 2011 Aug;32(8):894-9. doi: 10.1002/humu.21517.
SNV Information	SIFT Score Pred	SIFT官方注释,D=DAMAGING,T=TOLERATED
SNV Information	POLYPHEN  Score Pred	POLYPHEN官方注释,B=benign,P=possibly damaging,D=probably damaging
SNV Information	MutationTaster  Score Pred	MutationTaster官方注释,A=disease_causing_automatic,D=disease_causing,N=polymorphism,P=polymorphism_automatic
SNV Information	Cadd_Raw	SNV危险性评分，with higher values indicating that a variant is more likely to be simulated (or "not observed") and therefore more likely to have deleterious effects。参考http://cadd.gs.washington.edu.Cutoff is usually set 4.
SNV Information	Cadd_Phred	CADD phred-scaled scores;Basically, for scaled scores, 10 means 10% percentile highest scores, 20 means 1% percentile highest scores, and 30% means 0.1% percentile highest scores. So we do find one <1% percentile variant in the dataset(Cadd_Raw=~4)
SNV Information	Dann	SNV危险性评分，DANN scores whole-genome variants by training a deep neural network (DNN). DNNs can capture non-linear relationships among features and are better suited than SVMs for problems with a large number of samples and features. It seems to perform much better than CADD。参考http://bioinformatics.oxfordjournals.org/content/31/5/761.Cutoff is usually set 0.93.
SNV Information	Eigen	SNV危险性评分，Eigen?score uses a spectral approach integrating functional genomic annotations for coding and noncoding variants. Unlike other approaches such as CADD/DANN, Eignen did not use labelled training data and should be considered as a unsupervised learning approach.参考：http://www.ncbi.nlm.nih.gov/pubmed/26727659?
SNV Information	dbscSNV_ADA_SCORE	基于AdaBoost 算法计算的突变对剪切位点产生影响的可能性评分。值越大，概率越大，常用阈值>0.6
SNV Information	dbscSNV_RF_SCORE	基于Random Forest 算法计算的突变对剪切位点产生影响的可能性评分。值越大，概率越大，常用阈值>0.6
SNV Information	Hrcr1	全基因组突变频率，Haplotype Reference Consortium database with 40 million variants from 32K samples in haplotype reference consortium
SNV Information	Kaviar	全基因组突变频率，Kaviar database with 170 million variants from 13K genomes and 64K exomes。参考http://db.systemsbiology.net/kaviar/
SNV Information	ESP6500	来源于ESP( NHLBI GO Exome Sequencing Project )全外显子测序的6500个样本中的频率，请参考:http://evs.gs.washington.edu/EVS/
SNV Information	ExAC03	65000 exome allele frequency data for ALL. version 0.3. Left normalization done
SNV Information	ExAC03_EAS	65000 exome allele frequency data for EAS (East Asian). version 0.3. Left normalization done
SNV Information	gnomAD_exome	The Genome Aggregation Database (gnomAD) is a resource developed by an international coalition of investigators.The data set spans 123,136 exome sequences and 15,496 whole-genome sequences from unrelated individuals sequenced as part of various disease-specific and population genetic studies. 
SNV Information	gnomAD_exome_EAS	gnomAD data for EAS (East Asian).
SNV Information	gnomAD_genome	The Genome Aggregation Database (gnomAD) is a resource developed by an international coalition of investigators.The data set spans 123,136 exome sequences and 15,496 whole-genome sequences from unrelated individuals sequenced as part of various disease-specific and population genetic studies. 
SNV Information	gnomAD_genome_EAS	gnomAD data for EAS (East Asian).
SNV Information	Predicted Protein Variants	根据参照基因序列预测的基因突变引起的蛋白的氨基酸的变化
SNV Information	Interpro_domain	protein domain for variants
SNV Information	OMIM	基因在OMIM数据库中的信息，格式为染色体区段|遗传模式|OMIM编号:疾病名称
SNV Information	HPO	The Human Phenotype Ontology (http://www.human-phenotype-ontology.org) aims to provide a standardized vocabulary of phenotypic abnormalities encountered in human disease. 
SNV Information	HGMD_site	位点在HGMD数据库中的疾病/表型|pubmed编号
SNV Information	HGMD_site_class	致病性分类。DM:Disease-causing mutations;DM?:denoting a probable/possible pathological mutation;DP:Disease-associated polymorphisms;FP:Functional polymorphisms;DFP:Disease-associated polymorphisms with supporting functional evidence;R:Retired record;
SNV Information	HGMD_gene	基因在HGMD数据库中的疾病/表型，格式为Phenotype1;Phenotype2;Phenotype3
SNV Information	InterVar	InterVar根据ACMG评分规则进行的突变临床解释信息
SNV Information	InterVar_evidence	InterVar根据ACMG评分规则进行的28项打分结果
SNV Information	Predicted Protein Variants(snpEFF)	使用snpEFF进行突变注释，注释信息为：Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO
SNV Information	MGI	MGI(http://www.informatics.jax.org/) is the international database resource for the laboratory mouse, providing integrated genetic, genomic, and biological data to facilitate the study of human health and disease.
SNV Information	ClinVar	NCBI疾病突变数据库,注释内容分别为CLINSIG,CLNDBN,CLNACC,CLNDSDB,CLNDSDBID
SNV Information	MalaCards	MalaCards(http://www.malacards.org/) is an integrated database of human maladies and their annotations, modeled on the architecture and richness of the popular GeneCards database of human genes.
SNV Information	COSMIC	Catalogue Of Somatic Mutations In Cancer. It includes somatic mutations reported in literature in various types of cancers
SNV Information	NCI60	The NCI-60 cell lines are the most frequently studied human tumor cell lines in cancer research.
SNV Information	VCF Filter	VCF过滤硬条件描述
SNV Information	GO_BP	Gene Ontology annotation BP : biological process
SNV Information	GO_MF	Gene Ontology annotation MF : Molecular Function
SNV Information	GO_CC	Gene Ontology annotation CC: Cellular Component
SNV Information	KEGG_Pathway	该基因在kegg的pathway数据库中的信息
SNV Information	5'Flanking Sequence	SNV 5’侧的50bp碱基
SNV Information	3'Flanking Sequence	SNV 3’侧的50bp碱基
SNV Information	Homology Hits	包含该位点在内的左右各50bp组成的序列在人基因组中的同源序列数目，同源性达到90%以上就算HIT到一次，HIT次数越多该位点可靠性越低。
SNV Information	HWE	哈迪-温伯格平衡(Hardy-Weinberg equilibrium,HWE)的计算值
SNV Information	GeneskyExonDB_Freq	天昊513样本外显子数据库频率
SNV Information	Case(0|1|2)	正常纯合(0)、杂合(1)、突变纯合(2)的Case样本数量
SNV Information	Control(0|1|2)	正常纯合(0)、杂合(1)、突变纯合(2)的Control样本数量
SNV Information	Control Percentage	Alt Allele在control样本中的比例
SNV Information	Alt Allele Freq	Alt Allele在所有样本中的比例
High Freq	SNV NO	Single Nucleotide Variant (SNV, 包括插入或缺失）编号
High Freq	Chrs	染色体号
High Freq	Position	在对应染色体上的位置，以Reference GRCh38.p12 Primary Assembly为准
High Freq	Ref Allele	位点在参照序列上的等位基因
High Freq	Alt Allele	位点的另一个等位基因
High Freq	SNP ID	位点在ncbi的dbsnp中的编号。如果编号上标注“STR”表示位点很有可能是STR，如果编号上标注“STR_around”或“INDEL_around”表示位点周围存在STR或INDEL
High Freq	Freq_Alt (1000)	非参照等位基因在千人基因组项目中的频率，参考“http://www.1000genomes.org/home”
High Freq	Gene	基因名称
High Freq	Gene Region	变异所在的基因区域（如exonic=coding；splicing=2bp aroudn splicing junction；ncRNA；UTR5=5'untranslated region；UTR3=3'untranslated region；intronic；upstream =1kb upstream of transcription strart site；downstream=1kb downstream of transcription end site；intergenic)
High Freq	Function	基因功能（包括：nonsynonymous SNV，frameshift insertion，nonframeshift insertion,Stop Gain,Stop Loss）
High Freq	Case(0|1|2)	正常纯合(0)、杂合(1)、突变纯合(2)的Case样本数量
High Freq	Control(0|1|2)	正常纯合(0)、杂合(1)、突变纯合(2)的Control样本数量
High Freq	Alt Allele Freq	Alt Allele的频率（不区分case、control）
High Freq	HWE	哈迪-温伯格平衡(Hardy-Weinberg equilibrium,HWE)的计算值（不区分case、control）
High Freq	p1	基因上以case和control的HOMR数量、HET数量、HOMA数量为统计项的独立卡方检验值
High Freq	p2	基因上以case和control的HOMR数量、HET+HOMA数量为统计项的独立卡方检验值（如果Alt Allele Freq大于0.5，则将p2和p3交换）
High Freq	p3	基因上以case和control的HOMR+HET数量、HOMA数量为统计项的独立卡方检验值（如果Alt Allele Freq大于0.5，则将p2和p3交换）
High Freq	p4	基因上以case和control的Ref Allele数量、Alt Allele数量为统计项的独立卡方检验值
High Freq	Additive model OR(95%CI)	加性模型Logistic回归OR值
High Freq	Additive model P-value	加性模型Logistic回归P值
High Freq	Dominant model OR(95%CI)	显性模型Logistic回归OR值
High Freq	Dominant model P-value	显性模型Logistic回归P值
High Freq	Recessive model OR(95%CI)	隐性模型Logistic回归OR值
High Freq	Recessive model P-value	隐性模型Logistic回归P值
High Freq	Allele model OR(95%CI)	Allele模型Logistic回归OR值
High Freq	Allele model P-value	Allele模型Logistic回归P值
