# correlation 
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
	"============== do OTU cor ===========
	Usage: 
	Rscript 16S.OTU.cor.r otufile outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
		      outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()

}

library(Hmisc)

otufile = normalizePath(args[1])
outputpath = normalizePath(args[2])
setwd(outputpath)


# ====================== prepare data =======================
otu = read.table(otufile, header = T, sep = "\t", row.names = 1, comment.char = "", quote = "", fill = T)
#取$OTUsize前100数据
Data <- otu[order(otu$OTUsize, decreasing = TRUE),][1:100,]

loc = grep("All", rownames(Data))
if(length(loc) > 0) Data = Data[-loc,]

loc = grep("Abundance|OTUsize", colnames(Data))
if(length(loc) > 0){
	OTUdata = Data[,1:(loc-1)]
	TAXdata = Data[,loc:ncol(Data)]
}

SamplesOTU = OTUdata[complete.cases(OTUdata),]

# =================================== COR ===================================
# 相关性分析 计算两两OTU之间的相关性
result = rcorr(t(SamplesOTU),type="spearman")


# 提取相关系数和P值
Result = data.frame(Source=rep(colnames(result$r)[1:(ncol(result$r)-1)], times = c((ncol(result$r)-1):1)), 
					Target=unlist(lapply(2:ncol(result$r), function(t){rownames(result$r)[t:ncol(result$r)]})), 
					Cor=result$r[lower.tri(result$r)], 
					Pvalue=result$P[lower.tri(result$P)])
cor = 0.8
pvalue = 0.05

ResultForPlot = Result[( abs(Result$Cor)> cor & Result$Pvalue < pvalue ), ]
#ResultForPlot = ResultForPlot[ !duplicated(ResultForPlot[,c(3,4)]), ]      # 去重复

ResultForPlot$sig = 0
ResultForPlot[ResultForPlot$Cor>0,5] = 1
ResultForPlot[ResultForPlot$Cor<0,5] = -1
write.table(ResultForPlot, paste("OTU.Correlation.Cor.", cor, "_P", pvalue, ".txt",sep=""), row.names=F, sep="\t", quote=F)

# ========================== Prepare Data For Plot ==========================
# nodes information for cytoscape
OTUs = unique(c(as.character(ResultForPlot[,1]), as.character(ResultForPlot[,2])))
Size = TAXdata[OTUs, 1]
taxs = TAXdata[OTUs, 3:ncol(TAXdata)]

OTU.size.for.Cytoscape = lapply(taxs[1:ncol(taxs)], function(t) cbind(OTUs, Size, as.character(t)))

for (i in 1:length(OTU.size.for.Cytoscape)){

	colnames(OTU.size.for.Cytoscape[[i]])=c("Id", "Size", names(OTU.size.for.Cytoscape)[i])
	OTU.size.for.Cytoscape[[i]][,3] = gsub("\"", "", OTU.size.for.Cytoscape[[i]][,3])      # 去除名称中的引号
	# write.table(as.data.frame(OTU.size.for.Cytoscape[[i]]), paste("OTU.size.", names(OTU.size.for.Cytoscape)[i], ".for.Cytoscape.txt",sep=""), sep="\t", row.names=F, quote=F)

	# 统计各分类包含的OTU数目 
	Table  = as.matrix(table(OTU.size.for.Cytoscape[[i]][,3]))

	# 无注释的及 Unassigned 统一命名为Unclassified，在后续cytoscape中单独着色
	index1 = grep("Unassigned", rownames(Table))
	index2 = which(rownames(Table) == "")
	unclassifiedTax = rownames(Table)[c(index1,index2)]
	OTU.size.for.Cytoscape[[i]][ OTU.size.for.Cytoscape[[i]][,3] %in% unclassifiedTax, 3] = "Unclassified"
	Table  = as.matrix(table(OTU.size.for.Cytoscape[[i]][,3]))

	if ( length(Table) < 20 ){     # 分类数目小于20，全画

		others = NULL

	}else{     # 分类数目大于20个，则相对丰度低于1%的合并为others
		a = OTU.size.for.Cytoscape[[i]]
		A = data.frame(OTU = unlist(a[,1]), Size = as.numeric(a[,2]), tax = unlist(a[,3]), stringsAsFactors = F )
		relative_abundance = apply(A, 1, function(t) as.numeric(t[2])/sum(as.numeric(A$Size)))
		others = OTU.size.for.Cytoscape[[i]][which(relative_abundance < 0.01), 1]

	}

	# 注意矩阵中元素的提取和替换
	OTU.size.for.Cytoscape[[i]][ OTU.size.for.Cytoscape[[i]][,1] %in% others, 3] = "Others"
	OTU.size.for.Cytoscape[[i]][ OTU.size.for.Cytoscape[[i]][,3] %in% "candidate_division_WPS-2", 3] = "WPS-2"
	OTU.size.for.Cytoscape[[i]][ OTU.size.for.Cytoscape[[i]][,3] %in% "candidate_division_WPS-1", 3] = "WPS-1"
	OTU.size.for.Cytoscape[[i]][ OTU.size.for.Cytoscape[[i]][,3] %in% "Cyanobacteria/Chloroplast", 3] = "Cyanobacteria"

	# 输出OTU、size、taxon
	write.table(OTU.size.for.Cytoscape[[i]], paste("OTU.size.", names(OTU.size.for.Cytoscape)[i], ".for.Cytoscape.classif.txt", sep=""), sep="\t", row.names=F, quote=F)

}

