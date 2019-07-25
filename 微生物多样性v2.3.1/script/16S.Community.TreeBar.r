#!/usr/local/bin/R
rm(list = ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== TreeBar ===========
	Usage: 
	Rscript 16S.Community.TreeBar.r taxonfile groupfile outputpath
		parameters ->
			taxonfile: [file -> always Community_structure/*.taxon.Abundance.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		    outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(vegan)

taxonfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

setwd(outputpath)

taxondata <- read.table(taxonfile, header = T, sep = "\t", row.names = 1, check.names = F, quote = "", fill = T,comment.char = "")
sampleinfo <- read.table(groupfile, header = F, sep = "\t")
colnames(sampleinfo) <- c("sample", "Group")

taxondata <- taxondata[,as.character(sampleinfo$sample)]
loc <- grep("All",rownames(taxondata))
if(length(loc)>0) taxondata <- taxondata[-loc,]

loc = grep("OTUsize|Abundance", colnames(taxondata))
if(length(loc) > 0) taxondata = taxondata[,1:(loc-1)]

# prepare data
rownames(taxondata) = gsub("\"","",rownames(taxondata))     # 去除名字中的双引号,全局替换
colnames(taxondata) = gsub("\"","",colnames(taxondata))
data1 = as.matrix(taxondata)
rownames(data1) = as.character(strsplit(rownames(data1),split="{.*}",perl=T))

loc = grep("All|p__unclassified|c__unclassified|o__unclassified|f__unclassified|g__unclassified|s__unclassified", rownames(data1))      # 去除 All和unclassified行

if(length(loc) > 0) {
	data1 = data1[-loc, ]
	data1 = matrix(as.numeric(data1),nrow=nrow(data1))
	rownames(data1) = rownames(taxondata)[-loc] 
	colnames(data1) = colnames(taxondata) 
}

if( !is.null(nrow(data1)) ){     # 只有一类的,不画

	# get relative.abundance
	data_relative_abundance = sapply(1:ncol(data1),function(x) data1[,x]/sum(data1[,x]))
	
	# 种类少于30个的全画，多于30个的，将相对丰度低的进行合并为Others进行展示！
	if(nrow(data_relative_abundance) < 30){

		data_for_plot = data_relative_abundance
		colnames(data_for_plot) = colnames(taxondata)
		data_for_plot = data_for_plot[order(data_for_plot[,1],decreasing = T),]

	}else{
 
		# 非在所有样本中相对丰度都大于1%的就归为others
		index = which(apply(data_relative_abundance > 0.01, 1, any)==T)
		A = as.matrix(data_relative_abundance[index,])
		A = cbind(A,as.matrix(apply(A,1,sum)))
		A = as.matrix(A[order(A[,ncol(A)],decreasing = T),])
		A = A[,-ncol(A)]
		data_for_plot = rbind(A,apply(as.matrix(data_relative_abundance[-index,]),2,sum))
		rownames(data_for_plot)[nrow(data_for_plot)] = "Others"

	}

	data_for_plot = as.matrix(data_for_plot)
	colnames(data_for_plot) = colnames(taxondata)
	
	# index = grep("No_Rank", rownames(data_for_plot))
	# data_for_plot = rbind(data_for_plot[-index,], data_for_plot[index,])
	# rownames(data_for_plot)[nrow(data_for_plot)] = "Unclassified"

	# set color
	mycol = c(119,132,147,454,89,404,123,529,463,552,28,54,84,100,558,43,652,31,610,477,256,588,99,81,503,104,562,76,96,495)
	mycol = colors()[rep(mycol,20)]
	color = mycol[1:nrow(data_for_plot)] 
	
	ordergroup = sort(levels(factor(sampleinfo$Group)))
	ordercolor = rep(1, length(sampleinfo$Group))
	for(i in 1:length(ordergroup)){
        ordercolor[sampleinfo$Group%in%ordergroup[i]] = color[i]
}

	# set width and height
	mr = max(nchar(rownames(data_for_plot)))
	if (mr < 20) mr = 20
	mc = max(nchar(colnames(data_for_plot)))

	nr = 0.10*ncol(data_for_plot)    # 只根据样本数设定图形高度

	if(nr < 15){

		height = 15
		if(ncol(data_for_plot) < 8){
			cutpoint = c(7.5-6*nr,12*nr,7.5-6*nr)       # 样本量少于8个，则压缩绘图区域，扩大两侧区域
		}else{
			cutpoint = c(0.1,height,0.1)          # 若样本量较大，则可以把两侧空白区域压缩
		}

	}else{

		height = nr
		cutpoint = c(0.1,height,0.1)          # 若样本量较大，则可以把两侧空白区域压缩
	}

	nrow = height*2   	# 标签行数
	ncol = ceiling(nrow(data_for_plot)/nrow)
	width = 8 + 0.15*mc + 0.13*mr*ncol       # 根据 bar labels字符长度和标签列数限宽度

	name = unlist(strsplit(taxonfile,"/"))
	pdf(paste(outputpath, paste(strsplit(name[length(name)], "\\.\\w+\\.xls"), "Community.TreeBar.pdf", sep = "."), sep="/"), width = width, height = height)
	layout(rbind(c(0,0,3), c(1,2,3), c(0,0,3)), heights = cutpoint, widths = c(2, 6+0.15*mc, 0.13*mr*ncol), respect = FALSE)    # 分面,
	
	# Tree
	par(mar = c(5,1,2,0.7*mc), mgp = c(3,0.8,0.6), xpd=TRUE)
	
	fga.dist = vegdist(t(data1), method = "bray")		# BC distance
	hc = hclust(fga.dist, method = "average")		# UPGMA
	
	hcd = as.dendrogram(hc)
	
	data_for_plot = data_for_plot[,as.character(hc$labels[hc$order])]
	plot(hcd, horiz = T, axes = TRUE, yaxs = "i", ylim = c(0, ncol(data_for_plot) + 1) , cex = 13, col = ordercolor[hc$order])
	
	# Barplot
	par(mar = c(5,0,2,3), yaxt = "n", bty = "n")
	plot(sample(1:100, ncol(data_for_plot), replace=T),  1:ncol(data_for_plot), xlim = c(0, 100), ylim = c(0, ncol(data_for_plot) + 1), yaxs = "i", type = "n", xlab = "Relative abundance (%)", ylab = "")
	for (i in 1:ncol(data_for_plot)) {
		start = 0
		for (j in 1:nrow(data_for_plot)) {
			step = data_for_plot[j, i] * 100
			rect(start, i-0.4, start + step, i+0.4, border = color[j], col = color[j])
			start = start + step
		}	
	}
	
	# legend
	par(mar = c(1,1.5,3,1), xpd = TRUE)
	plot.new()
	legend("center", legend = rownames(data_for_plot), fill = color, ncol = ncol, bty = "n", cex = 1.8)
	dev.off()
	
}

