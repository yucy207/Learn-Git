# Differential Analysis Using Wilcoxon rank-sum test and ANOVA
rm(list=ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 5){

	cat(
	"Usage: 
	Rscript FUNGuild.ANOVA.DEA.r guildfile groupfile Threshold Type outputpath
		parameters ->
			  guildfile: [file -> always guilds.xls];
			  groupfile: [file -> always sample.groups without header and only two column];
			       Type: [character -> P or FDR, always P];
			  Threshold: [number -> Threshold for P to choose sig taxa, always 0.05];
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()
	
}

library(plyr)
library(pheatmap)

guildfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
Type = toupper(args[3])
Threshold = as.numeric(args[4])
outputpath = normalizePath(args[5])


# ÎïÖÖÊý¾Ý
guilddata = read.table(guildfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote = "", fill = T)

# Ñù±¾ÐÅÏ¢
SampleInfo = read.table(groupfile, header = F, sep="\t")
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("Sample", "Group")

# sort by groups
SampleInfo = SampleInfo[order(SampleInfo$Group),]
SampleInfo$Group = as.factor(SampleInfo$Group)
Groups = levels(SampleInfo$Group)

# prepare data
count = guilddata[,rownames(SampleInfo)]
count = count[rowSums(count) != 0 , ]

# Ïà¶Ô·á¶È
DataForDEA = apply(count, 2, function(d) d/sum(as.numeric(d)))

if(length(Groups) > 2){	

	dir1 = paste(outputpath, "all_group_different_analysis",sep = "/")
	system(paste("mkdir -p ", dir1, sep = ""))
	setwd(dir1)
	
	# compare	
	res = sapply(1:nrow(DataForDEA), function(g) {
		mod = lm(DataForDEA[g,] ~ SampleInfo$Group)
		am = anova(mod)
		F1 = am[4][1, 1]   # F value
		P1 = am[5][1, 1]   # P value
		FP = cbind(F1, P1)
	})

	res = t(res)
	rownames(res) = rownames(DataForDEA)
	colnames(res) = c("Fvalue", "Pvalue")	

	fdr = p.adjust(res[,2], method = "BH")

	means = sapply(Groups, function(g) apply(as.data.frame(DataForDEA[,as.character(SampleInfo[as.character(SampleInfo$Group) %in% g, 1])]), 1, mean))
	colnames(means) = paste("Mean.In.", Groups, sep = "")
	
	result = data.frame(Guild = rownames(DataForDEA), means, res, FDR = fdr)
	result = result[order(result$Pvalue),]
	write.table(result,"ANOVA_test_result.xls", sep = "\t", row.names = F, quote = F)
	
	# export significant result
	if (Type == "P"){

		SigRes = subset(result, Pvalue < Threshold)

	}else{

		SigRes = subset(result, FDR < Threshold)
	
	}

	##########plot##############
	# prepare plot data
	if(nrow(SigRes) >= 1){

		for (i in 1:nrow(SigRes)){

			tax = rownames(SigRes)[i]
			p_value = (SigRes[which(SigRes$Guild %in% tax), ])$Pvalue
			if(p_value >= 0.01){p = "*"}
			if(p_value < 0.01 & p_value >= 0.001){p = "**"}
			if(p_value < 0.001 & p_value >= 0.0001){p = "***"}
			if(p_value < 0.0001){p = "****"} 
			rownames(SigRes)[i] = paste(rownames(SigRes)[i],p,sep = "")
		}

		index = which(colnames(SigRes) %in% "Fvalue")
		DataForplot = t(SigRes[,3:index-1] * 100)

		cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
		bar_cols = colors()[rep(cols,20)][1:length(Groups)]
		height_barplot = 0.15 * nrow(SigRes) * length(Groups) + 0.15 * (nrow(SigRes) - 1) + 1 + 0.3
		height =  height_barplot + 2

		#layout
		width = 12 
		pdf("ANOVA_diff_0.05.pdf",width = width ,height = height)

		layout(matrix(c(1,2), nrow = 2, byrow = T),heights = c(2,height_barplot))
		par(mai = c(0.1,2,0.2,2), xaxs = "i", yaxs = "i")
		plot.new() ###legend绘图列
		legend("top", col = bar_cols, lty = 1, lwd = 10, xpd = T, cex = 1, ncol = 6, legend = Groups, bty = "n")
		
		#plot(1:10, type = "n", yaxt = "n", xlab = "", xaxt = "n", ylab = "", bty = "n")
		#n = 0
		#for(i in 1:length(Groups)){
		#	rect(xleft = n *1.6 + 1 , ybottom = 4, xright = n *1.6 + 1.5, ytop = 6, col = bar_cols[i],lwd = 1.5)
		#	text(x = (n *1.6 + 1.5)*1.02, y = 5, adj = 0, labels = Groups[i], cex = 1)
		#	n  = n + 1
		#}
	
		x = 0.1* max(nchar(colnames(DataForplot)))
		if (x > 4){
			xx = x
		}else{
			xx = 4
		}		
		par(las = 1, mai = c(1, xx, 0.3, 0.5),  xaxs = "i", yaxs = "i", mgp = c(3, 2, 1), lend = 2,lwd = 1.5)

		#barplot
		wd = 0.8
		if(max(DataForplot) >= 1){

			x_max <- ceiling(max(DataForplot))
			a = barplot(DataForplot, beside = T, horiz = T, width = wd ,col = bar_cols, axes = FALSE,xlim = c(0,x_max))
			cols = 1:ncol(a)
			#添加每个柱子的灰度背景
			rect(xleft = 0, ybottom = a[1, cols %% 2 == 1] - wd / 2, xright = par("usr")[2], ytop = a[length(Groups), cols %% 2 == 1] + wd / 2  , col = "gray91", border = NA)
			barplot(DataForplot, beside = T, horiz = T, width = wd, col = bar_cols , add = TRUE, axes = FALSE,xlim = c(0,x_max))
			axis(side = 1, at = c(0, x_max), labels = c(0, x_max))

		}else{

			x_max<-format(max(DataForplot),scientific=TRUE,digit=2)
			a = barplot(DataForplot, beside = T, horiz = T, width = wd ,col = bar_cols, axes = FALSE,xlim = c(0,max(DataForplot)))
			cols = 1:ncol(a)
			#添加每个柱子的灰度背景
			rect(xleft = 0, ybottom = a[1, cols %% 2 == 1] - wd / 2, xright = par("usr")[2], ytop = a[length(Groups), cols %% 2 == 1] + wd / 2  , col = "gray91", border = NA)
			barplot(DataForplot, beside = T, horiz = T, width = wd, col = bar_cols , add = TRUE, axes = FALSE,xlim = c(0,max(DataForplot)))
			axis(side = 1, at = c(0, max(DataForplot)), labels = c(0, x_max))
		}
		title(xlab = "Mean proportion (%)")

		dev.off()
	
	}

}else{

	stop("sample groups less than 2, please check it!\n")

}



