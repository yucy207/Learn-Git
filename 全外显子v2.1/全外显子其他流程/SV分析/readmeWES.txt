Chrs	染色体号
Start	起始位置（参考基因组版本hg19）
End	终止位置（由于BND易位片段非常大，而且部分跨染色体，所以无法用染色体的起始位置和终止位置表示发送变异的片段，用起始位置代替）
SVlength	片段长度
SVtype	结构变异的类型：缺失DEL，重复DUP，易位BND,倒位INV
Strands Info	易位BND的方向信息
ControlDBFreq	  10个天昊全基因组对照样本中的变异频率
CaseFreq	  项目case样本中变异频率(当样本变异片段数SU>=5时，认为结果可信，否则认为未发生变异)
ControlFreq	  项目control样本中变异频率
Start Homology	起始位点的同源性（起始位点两侧各扩展50bp并提取参考基因组序列，然后基于megablast进行比对，统计比对结果中比对长度大于90的记录个数作为同源性量化数值）
End Homology	终止位点的同源性
Gene	结构变异所在的基因
Gene Region	结构变异位于基因的区域
Cover Region	位于同一条染色体的BND变异覆盖基因的区域
HGMD	基因在HGMD数据库中的疾病/表型，格式为Phenotype1;Phenotype2;Phenotype3
OMIM	基因在OMIM数据库中的信息，格式为Cytogenetic location|MIM Number|Disorders
MalaCards	MalaCards(http://www.malacards.org/) is an integrated database of human maladies and their annotations, modeled on the architecture and richness of the popular GeneCards database of human genes.
cytoBand	染色体条带区间
dgvFreq	DGV数据库中片段拷贝数的频率
iscaPathGainCum	ISCA（International Standards for Cytogenomic Arrays )数据库中拷贝数增加的样本数量
iscaPathLossCum	ISCA数据库中拷贝数减少的样本数量
iscaLikelyPathogenic	ISCA数据库中可能导致疾病
iscaPathogenic	ISCA数据库中导致疾病
iscaCuratedPathogenic	ISCA数据库中确认能导致疾病
CNVD	Copy Number Variation in Disease数据库注释，给出突变类型，相关基因以及对应的疾病
DECIPHER	DECIPHER数据库注释，给出缺失与多拷贝样本的频率以及样本总量
dbVar	NCBI's database of genomic structural variation,列出dbVar注释区段的拷贝数变异频数
样本	支持该结果变异的read数量

SV 特殊说明	文件包含两个表格SV Filter；SV Original。其中SV Original为原始SV分析结果（去掉样本中深度低于5的记录）。SV Filter为对原始结果进行过滤，过滤条件为：（1）如果倒位易位发生在同一个基因上，要求覆盖区域为exonic；如果是不同染色体或基因间倒位易位，则要求变异位点不能同时处于处于基因间区intergenic;(2)ControlFreq对照数据库变异频率<=0.1。(3)两个断点的同源性不能均大于10 (4)支持变异的reads数大于10。 以上条件为并列条件
CNV 特殊说明	文件包含两个表格CNV Filter；CNV Original。其中CNV Original为原始CNV分析结果（去掉样本中深度低于5的记录）。CNV Filter为对原始结果进行过滤，过滤条件为：（1）CNV变异区域位于exonic (2)ControlFreq对照数据库变异频率<=0.1。(3)两个断点的同源性不能均大于10 (4)支持变异的reads数大于10 。 以上条件为并列条件
