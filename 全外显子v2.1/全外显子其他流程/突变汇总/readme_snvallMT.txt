Sheet Name	标注与说明	
SNV	SNV NO	Single Nucleotide Variant (SNV, 包括插入或缺失）编号
SNV	Gene	变异位点所在的基因。
SNV	Ref Allele	位点在参照序列上的等位基因。
SNV	Alt Allele	位点的另一个等位基因。
SNV	Chrs	染色体号。
SNV	Position	在对应染色体上的位置，以Reference GRCh37.p5 Primary Assembly为准。
SNV	Strand Orientation	显示基因在染色体上的排列方向，+表明基因序列同参考序列一致，-表明与参考序列互补。
SNV	Gene Region	变异所在的基因区域。
SNV	Function	基因功能（包括：nonsynonymous SNV，frameshift insertion，nonframeshift insertion,Stop Gain,Stop Loss）
SNV	Predicted Protein Variants	根据参照基因序列预测的基因突变引起的蛋白的氨基酸的变化
SNV	Codon Mutation	外显子密码子突变，书写格式：基因：转录本：codon.密码子位置.突变碱基处于密码子的位置(1/2/3).参考密码子->突变密码子
SNV	Mutation Pattern	在点突变基础上增加5’base和3’base构成三碱基突变模式，书写格式：参考->突变
SNV	tRNA_anticodon	tRNA反义密码子
SNV	phastCons100way	phastCons100way数据库保守性注释
SNV	phyloP100way	phyloP100way数据库保守性注释
SNV	HmtVar_Asia_healthy_freq	HmtVar数据库亚洲健康人群变异频率
SNV	HmtVar_Asia_pathologic_freq	HmtVar数据库亚洲患病人群变异频率
SNV	MT_Control_Freq	天昊200样本MT数据库中Alt平均频率
SNV	mitomap_disease	基于MITOMAP的线粒体疾病数据库注释。https://www.mitomap.org/MITOMAP
SNV	PolyPhen2	PolyPhen v.2.2.2数据库的分类预测，Polymorphism Phenotyping 是一种用来预测非同义突变造成的氨基酸改变是否影响功能的算法，参考http://genetics.bwh.harvard.edu/pph2
SNV	SIFT	SIFT v.5.0.3数据库的分类预测，SIFT (Sorting Intolerant From Tolerant) 是另外一个用于预测氨基酸改变是否影响功能的程序，参考http://sift.bii.a-star.edu.sg
SNV	FatHmm	FatHmm v.2.2(未加权版本)数据库的分类预测，FatHmm可用于预测突变的功能性影响，参考http://fathmm.biocompute.org.uk
SNV	FatHmmW	FatHmm v.2.3(加权版本)数据库的分类预测，同FatHmm
SNV	PROVEAN	 PROVEAN v.1.3数据库的分类预测，PROVEAN (Protein Variation Effect Analyzer) 可用于预测突变是否对蛋白质的生物学功能产生影响，参考http://provean.jcvi.org
SNV	MutationAssessor	MutationAssessor v.2.0数据库的分类预测，类似PROVEAN，评估突变对蛋白质功能的影响，参考http://mutationassessor.org
SNV	EFIN_SP	EFIN SP (使用SwissProt数据集进行训练)的分类预测，EFIN是一种基于蛋白质保守性，使用随机森林方法对氨基酸突变是否与疾病有关进行评估的工具，参考http://paed.hku.hk/efin
SNV	EFIN_HD	EFIN HD (使用HumDiv数据集进行训练)的分类预测，同EFIN_SP
SNV	CADD	CADD v.1.2数据库的分类预测，SNV的危险性评分，参考http://cadd.gs.washington.edu
SNV	Meta-SNP, PANTHER, PhD-SNP, SNAP	四个可预测与疾病相关突变的软件，参考http://snps.biofold.org/meta-snp
SNV	CAROL	CAROL consensus method的分类预测，突变的危害性预测，参考http://www.sanger.ac.uk/science/tools/carol
SNV	Condel	Condel consensus method的分类预测，类似CAROL，突变的危害性预测，参考http://bg.upf.edu/fannsdb
SNV	COVEC_WMV	COVEC v.0.4 Weighted Majority Rule consensus method的分类预测，是一种整合第三方软件结果对突变影响进行注释的工具，参考http://sourceforge.net/projects/covec
SNV	MtoolBox	MtoolBox consensus method的分类预测，是一种用于线粒体突变有害性预测的工具，参考https://github.com/mitoNGS/MToolBox/blob/master/MToolBox/data/patho_table.txt
SNV	PolyPhen2_transf	通过TransFIC v.1.0转换的PolyPhen2分类预测结果
SNV	SIFT_transf	通过TransFIC v.1.0转换的SIFT分类预测结果
SNV	MutationAssessor_transf	通过TransFIC v.1.0转换的MutationAssessor分类预测结果
SNV	MutationTaster	MutationTaster官方注释(disease_causing, disease_causing_automatic, polymorphism)
SNV	PCT_AltFreq>=0.1	突变频率大于等于0.1的样本比例
SNV	PCT_AltFreq>=freq_cutoff	突变频率大于等于特定频率阈值的样本比例（具体阈值参看报告）
SNV	AltFreq(0.1~0.9)	突变频率在0.1~0.9之间的样本。
SNV	AltFreq>=0.9	突变频率大于等于0.9的样本。
SNV	5'Flanking Sequence	SNV 5’侧的50bp碱基。
SNV	3'Flanking Sequence	SNV 3’侧的50bp碱基。
SNV	Homology Hits	包含该位点在内的左右各50bp组成的序列在人基因组中的同源序列数目，同源性达到90%以上就算HIT到一次，HIT次数越多该位点可靠性越低。
