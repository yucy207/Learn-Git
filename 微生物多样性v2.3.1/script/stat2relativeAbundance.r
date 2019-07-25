rm(list = ls())

args = commandArgs(trailingOnly=T) 
if(length(args) != 2){

	cat("Usage: 
	Rscript stat2relativeAbundance.r otustatfile outputpath
		parameters -> 
			otustatfile: [file -> always otu.tax.*.stat.xls derived from 16S.OTU.Modify.pl];
		         outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
}
 
otustat = normalizePath(args[1])
outputpath = normalizePath(args[2])

filename = unlist(strsplit(otustat, "/|\\.xls"))[length(unlist(strsplit(otustat, "/|\\.xls")))]

# otustat = normalizePath("/home/yangkj/work/metagenomics/17B0204A/newqc/tmp/testLefse/otu.tax.0.03.stat.xls")
# outputpath = normalizePath("/home/yangkj/work/metagenomics/17B0204A/newqc/tmp/testLefse")

# =============================== Data prepare ===============================
otu_stat = read.table(otustat, header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "", quote = "") 

# 去除Abundance列
loc = grep("Abundance|OTUsize", colnames(otu_stat))
if(length(loc) > 0) {
	otu_stat1 = as.matrix(otu_stat[,1:(loc-1)])
	rownames(otu_stat1) = rownames(otu_stat)
	colnames(otu_stat1) = colnames(otu_stat)[1:(loc-1)]
	otu_stat = otu_stat1
}


# 总的丰度
index = grep("\\|", rownames(otu_stat))
if(length(index) == (nrow(otu_stat)-1)){
	totalAbun = as.matrix(otu_stat[-index,])
}else{
	totalAbun = apply(as.matrix(otu_stat[-index,]), 2, sum)
}

# 相对丰度
ReAbun = sapply(1:ncol(otu_stat), function(x) 100*otu_stat[,x]/totalAbun[x])
rownames(ReAbun) = rownames(otu_stat)
colnames(ReAbun) = colnames(otu_stat)

ReAbunOut = cbind(rownames(ReAbun), ReAbun)
colnames(ReAbunOut)[1] = "Taxlevel"
write.table(ReAbunOut, paste(outputpath, paste(filename, ".RA.xls", sep = ""), sep = "/"), row.names = F, sep = "\t", quote = F)
