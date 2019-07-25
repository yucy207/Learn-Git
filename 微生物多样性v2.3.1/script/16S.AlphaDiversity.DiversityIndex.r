rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== Alpha Diversity Index Boxplot ===========
	Usage: 
	Rscript 16S.AlphaDiversity.DiversityIndex.r otufile groupfile outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()

}

library(vegan)
library(entropart)

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

setwd(outputpath)

# ====================== prepare data =======================
SampleInfo = read.table(groupfile, header = F, sep = "\t",quote="")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]

otudata = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote = "", fill = T)
otu_data = as.matrix(otudata[,rownames(SampleInfo)])
colnames(otu_data) = as.character(SampleInfo[,1])
rownames(otu_data) = rownames(otudata)

loc = grep("All", rownames(otu_data))
if (length(loc) > 0) {
	otu_data = as.matrix(otu_data[-loc,])
	colnames(otu_data) = as.character(SampleInfo[,1])
}

# alpha diversity index
######################################################################
#alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson")
#source("/home/zhengy/bin/modules/script/estimate_richness_new.r")
#er = estimate_richness_new(t(otu_data), measures = alpha_meas)

OTU = t(otu_data)
outlist = vector("list")
measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson")
estimRmeas = c("Chao1", "Observed", "ACE")
if (any(estimRmeas %in% measures)) {
	outlist <- c(outlist, list(t(data.frame(estimateR(OTU))))) ### estimateR library(vegan)
}
if ("Shannon" %in% measures) {
	outlist <- c(outlist, list(shannon = diversity(OTU, index = "shannon"))) ### diversity library(vegan)
}
if ("Simpson" %in% measures) {
	#outlist <- c(outlist, list(simpson = diversity(OTU, index = "simpson")))
	outlist <- c(outlist, list(simpson = 1-diversity(OTU, index = "simpson")))
}
if ("InvSimpson" %in% measures) {
	outlist <- c(outlist, list(invsimpson = diversity(OTU, index = "invsimpson")))
}
out = do.call("cbind", outlist)
renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")
names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", "simpson", "invsimpson", "fisher")
namechange = intersect(colnames(out), names(renamevec))
colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
er <- as.data.frame(out)

### Coverage library(entropart)
options(warn = (-1))     # suppress warnings
Coverage = apply(otu_data, 2, function(x) Coverage(x, CheckArguments = FALSE))     # probability for a species of the community to be observed in the actual sample

# combine result
Diversity = cbind(rownames(er), er[,-c(3,5)], Coverage)
colnames(Diversity)[1] = "Samples"
write.table(Diversity, "alpha.diversity.index.xls", row.names = F, quote = F,sep = "\t")
######################################################################

mdf = reshape2::melt(Diversity, measure.vars = colnames(Diversity)[-1])
mdf = cbind(rep(as.character(SampleInfo$Group), ncol(Diversity)-1), mdf)
colnames(mdf)[1] = "Group"
 
cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
mycol= colors()[rep(cols,20)][1:length(levels(mdf$Group))]

box_plot <- function(value,group,name,col){
	data = data.frame( value = value , group = group)
	b = boxplot( value ~ group, data, col = col, lwd = 2,xaxt = "n")
	#text( x = c(1:length(levels(group))), y = par("usr")[3], srt = 45, cex = 1, labels = levels(group), xpd = TRUE,pos = 1,offset = 1)
	text( x = c(1:length(levels(group))), y = (par("usr")[3]-(par("usr")[4] - par("usr")[3])*0.05), adj = c(1,1),srt = 45, cex = 1, labels = levels(group), xpd = TRUE)
	title(name)

}

Groups = unique(SampleInfo$Group)

########all groups#############
if(length(Groups) >= 3 & length(SampleInfo$Sample) > length(Groups) ){

	# Kruskal-Wallis test for multiple groups
	dir1 = paste(outputpath, "All_groups", sep = "/")
	system(paste("mkdir -p ", dir1, sep = ""))
	setwd(dir1)

	p.value = sapply(2:ncol(Diversity), function(d) kruskal.test(Diversity[,d], SampleInfo$Group)$p.value)	
	fdr = p.adjust(p.value, method = "BH")
	output = cbind(colnames(Diversity)[2:ncol(Diversity)], p.value, fdr)
	colnames(output)[1] = "Index"

	write.table(as.matrix("Difference comparison using Kruskal-Wallis rank sum test"), "Kruskal-Wallis_Sig_test.xls", row.names = F, col.names = F, quote =F)
	write.table(output, "Kruskal-Wallis_Sig_test.xls", append = T, row.names = F, quote = F, sep = "\t", na = "")

	# plot box
	pdf("AlphaDiversity.pdf", width = 16, height = 12)
	par(mar = c(8,5,6,2),mfrow = c(2, 3), cex = 1, cex.axis = 1, cex.main = 1.8, las = 1,lwd = 2)

	for(i in 1:nrow(output)){

		value = mdf$value[mdf$variable == output[i,1]]
		group = mdf$Group[mdf$variable == output[i,1]]
		p = as.numeric(output[i,2])
		if(p >= 0.05){name = output[i,1]}
		if(p >= 0.01 & p < 0.05){name = paste(output[i,1],"*",collapse =" ")}
		if(p < 0.01 & p >= 0.001){name = paste(output[i,1],"**",collapse =" ")}
		if(p < 0.001 & p >= 0.0001){name = paste(output[i,1],"***",collapse =" ")}
		if(p < 0.0001){name = paste(output[i,1],"****",collapse =" ")}

		box_plot(value,group,name,mycol)
	
	}
	dev.off()

}


########pairwise groups#############

if(length(Groups) >= 2 & length(SampleInfo$Sample) > length(Groups)){
	# compare significant
	dir2 = paste(outputpath, "Pairwise_groups", sep = "/")
	system(paste("mkdir -p ", dir2, sep = ""))
	setwd(dir2)

	pairwise = function (compare.levels, level.names) {

		ix <- setNames(seq_along(level.names), level.names)
		pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
			function(k) {
				i <- ivec[k]
				j <- jvec[k]
				if (i > j) 
					compare.levels(i, j)
				else NA
			}))
		pp
	}
	
	

	write.table(as.matrix("Difference comparison using Wilcoxon rank sum test; P value"), "Wilcoxon_pairwise_Sig_test_PVAL.xls", row.names = F, col.names = F, quote = F)
	write.table(as.matrix("Difference comparison using Wilcoxon rank sum test; FDR, adjustment method: BH"), "Wilcoxon_pairwise_Sig_test_FDR.xls", row.names = F, col.names = F, quote = F)

	for (i in 2:ncol(Diversity)){

		x = Diversity[,i]
		names = colnames(Diversity)[i]
		compare.levels <- function(i, j) {
			xi <- x[as.integer(as.factor(SampleInfo$Group)) == i]
			xj <- x[as.integer(as.factor(SampleInfo$Group)) == j]
			wilcox.test(xi, xj, paired = FALSE)$p.value
		}
		PVAL = pairwise(compare.levels, levels(as.factor(Groups)))
		
		PVAL = cbind(rownames(PVAL), PVAL)
		colnames(PVAL)[1] = "Groups"
		
		# Wilcoxon rank sum test for pairwise groups
		FDR <- pairwise.table(compare.levels, levels(as.factor(Groups)), p.adjust.method = "BH")
		FDR = cbind(rownames(FDR), FDR)
		colnames(FDR)[1] = "Groups"
		
		write.table(as.matrix(paste("\nsignificance of ", names, " difference between pairwise groups", sep = "")), "Wilcoxon_pairwise_Sig_test_PVAL.xls", row.names = F, col.names = F, append = T, quote =F)
		write.table(PVAL, "Wilcoxon_pairwise_Sig_test_PVAL.xls", append = T, row.names = F, quote = F, sep = "\t", na = "")
		write.table(as.matrix(paste("\nsignificance of ", names, " difference between pairwise groups", sep = "")), "Wilcoxon_pairwise_Sig_test_FDR.xls", row.names = F, col.names = F, append = T, quote =F)
		write.table(FDR, "Wilcoxon_pairwise_Sig_test_FDR.xls", append = T, row.names = F, quote = F, sep = "\t", na = "")
		
	}

	#plot
	group2col = data.frame(as.character(levels(as.factor(Groups))),mycol,stringsAsFactors=FALSE)
	colnames(group2col) = c("gro","col")

	for(m in 1:(length(Groups)-1)){

		for (n in (m+1):length(Groups)){

			#####创建两两差异分析的文件夹###########
			dir3 <- paste(dir2, paste(Groups[m], "_vs_", Groups[n], sep = ""), sep = "/")
			system(paste("mkdir -p ", dir3, sep = ""))
			setwd(dir3)
			#####创建新的分组信息###################           
			SampleInfo1 <- rbind(SampleInfo[SampleInfo$Group == Groups[m],], SampleInfo[SampleInfo$Group == Groups[n],])
			mdf1 = mdf[which(mdf$Group %in% SampleInfo1$Group),]
			mycol1 = group2col$col[group2col$gro %in% as.character(unique(SampleInfo1$Group))]

			pdf("AlphaDiversity.pdf", width = 16, height = 12)
			par(mar = c(8,5,6,2),mfrow = c(2, 3), cex = 1, cex.axis = 1, cex.main = 1.8, las = 1,lwd = 2)
			for (i in 2:ncol(Diversity)){

				x = Diversity[,i]
				names = colnames(Diversity)[i]
				xm <- x[as.integer(as.factor(SampleInfo$Group)) == m]
				xn <- x[as.integer(as.factor(SampleInfo$Group)) == n]
				p = wilcox.test(xm, xn, paired = FALSE)$p.value
				if (p == "NaN") {p = 1}
				if(p >= 0.05){name1 = names}
				if(p >= 0.01 & p < 0.05){name1 = paste(names,"*",collapse =" ")}
				if(p < 0.01 & p >= 0.001){name1 = paste(names,"**",collapse =" ")}
				if(p < 0.001 & p >= 0.0001){name1 = paste(names,"***",collapse =" ")}
				if(p < 0.0001){name1 = paste(names,"****",collapse =" ")}

				value1 = mdf1$value[mdf1$variable == names]
				#group = mdf1$Group[mdf1$variable == names]
				group1 = as.factor(as.character(SampleInfo1$Group))

				box_plot(value1,group1,name1,mycol1)
			}
			dev.off()
		}
	}
}




















