rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== plot species accumulation curve ===========
	Usage: 
	Rscript 16S.AlphaDiversity.Specaccum.r otufile groupfile outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
			  groupfile:[file -> always sample.groups];
		      outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(vegan)

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

otudata = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
SampleInfo = read.table(groupfile, header = F, sep = "\t",quote="")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
groups = unique(as.character(SampleInfo[,2]))

otudata = as.matrix(otudata[,rownames(SampleInfo)])
colnames(otudata) = as.character(SampleInfo[,1])

loc = grep("All", rownames(otudata))
if(length(loc) > 0) otudata = as.matrix(otudata[1:(loc-1),])

####所有样本一起绘图
if( ncol(otudata) < 3 ){

	stop("Sorry: You must have at least 3 samples to plot Specaccum!\n")

}else{
	otudata = otudata[rowSums(otudata) != 0, ]
}

sp1 <- specaccum(t(otudata))
sp2 <- specaccum(t(otudata), "random")

pdf(paste(outputpath, "all.sample.specaccum.pdf", sep = "/"))
par(mar = c(3, 3, 3, 1), mgp = c(1.7, 0.4, 0))
plot(sp1, ci.type = "bar", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue", xlab = "Number of sample", ylab = "OTUs", main = "species accumulation curve", cex.main = 1.1)
boxplot(sp2, col = "lightyellow", add = TRUE, pch = "+")
dev.off()

####分组绘图
for (i in 1:length(groups)) {

	sample = as.character(SampleInfo[which(as.character(SampleInfo[,2]) == groups[i]), 1])

	if(length(sample) >= 3){
		otudata1 = otudata[,sample]
		otudata1 = otudata1[rowSums(otudata1) != 0, ]

		sp1 <- specaccum(t(otudata1))
		sp2 <- specaccum(t(otudata1), "random")

		pdf(paste(outputpath, paste(groups[i], ".specaccum.pdf", sep = ""), sep = "/"))
		par(mar = c(3, 3, 3, 1), mgp = c(1.7, 0.4, 0))
		plot(sp1, ci.type = "bar", col = "blue", lwd = 2, ci.lty = 0, ci.col = "lightblue", xlab = "Number of sample", ylab = "OTUs", main = "species accumulation curve", cex.main = 1.1)
		boxplot(sp2, col = "lightyellow", add = TRUE, pch = "+")
		dev.off()
	}
	
}