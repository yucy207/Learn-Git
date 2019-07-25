# correlation between two dissimilarity matrices using Mantel Test
# 检验样本群落分布与环境因子整体之间的相关性

rm(list=ls())

args = commandArgs(trailingOnly=T)
if(length(args) != 3){

	cat(
	" ============== correlation between env and otu using Mantel Test ============
	Usage: 
	Rscript 16S.EnvironmentFactors.MantelTest.r envfile otufile outputpath
		parameters ->
			envfile: [file -> always env.xls with row env.xls and colum environmental factors];
			otufile: [file -> always otu.tax.0.03.xls];
		     outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

library(vegan)

envfile = normalizePath(args[1])
otufile = normalizePath(args[2])
outputpath = normalizePath(args[3])

envdata = read.table(envfile, header=T, sep="\t", row.names=1)
otudata = read.table(otufile, header=T, sep="\t", row.names = 1, check.name = F, comment.char = "", quote = "")
otudata = otudata[, rownames(envdata)]

fga.dist = vegdist(t(otudata), method = "bray")
env.dist = vegdist(scale(envdata), method = "euclidean")  # scale 对数据进行标准化!     
result  = mantel(fga.dist, env.dist, permutations = 9999)   # method="spearman" 默认方法是 pearson

ResultForOutput = data.frame(MatrixA = "Community", MatrixB = "Environmental factors", r = result$statistic, P_value = result$signif, Permutations = 9999)
#write.csv(ResultForOutput, paste(outputpath, "MantelTestResult.csv", sep = "/"), row.names=F, quote=F)
write.table(ResultForOutput, paste(outputpath, "MantelTestResult.xls", sep = "/"), sep = "\t", row.names = F, quote = F)
