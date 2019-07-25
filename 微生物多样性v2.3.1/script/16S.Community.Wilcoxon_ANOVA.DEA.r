# Differential Analysis Using Wilcoxon rank-sum test and ANOVA
rm(list=ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 5){

	cat(
	"Usage: 
	Rscript Meta.Wilcoxon_ANOVA.DEA.r taxonfile groupfile Threshold Type outputpath
		parameters ->
			  taxonfile: [file -> always *.taxon.Abundance.xls with row taxon and column samples];
			  groupfile: [file -> always sample.groups without header and only two column];
			       Type: [character -> P or FDR, always P];
			  Threshold: [number -> Threshold for P to choose sig taxa, always 0.05];
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()
	
}

library(plyr)
library(pheatmap)

taxonfile = normalizePath(args[1])
groupfile = normalizePath(args[2])
Type = toupper(args[3])
Threshold = as.numeric(args[4])
outputpath = normalizePath(args[5])

setwd(outputpath)

# 读取物种丰度表
taxondata = read.table(taxonfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote = "", fill = T)

# 读取分组文件
SampleInfo = read.table(groupfile, header = F, sep="\t")
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("Sample", "Group")

# sort by groups
#Groups = levels(as.factor(as.character(SampleInfo$Group)))
Groups = unique(as.character(SampleInfo$Group))

taxonname = unlist(strsplit(unlist(strsplit(taxonfile,"/"))[length(unlist(strsplit(taxonfile,"/")))], "\\."))[1]

#新建目录
dir = paste(outputpath, taxonname, sep = "/")
system(paste("mkdir -p ", dir, sep = ""))
setwd(dir)

# prepare data
count = taxondata[,rownames(SampleInfo)]
rownames(count) = gsub("\"", "", rownames(count)) 
colnames(count) = gsub("\"", "", colnames(count))
rownames(count) = as.character(strsplit(rownames(count), split = "{.*}", perl = T))

count = count[rowSums(count) != 0 , ]
loc = which(rownames(count) == "All")       
if(length(loc) > 0) count = count[-loc, ]

# 获得相对丰度文件
re.count = apply(count, 2, function(d) d/sum(as.numeric(d)))

# ANOVA分析
if(length(Groups) > 2){	
	
	res = sapply(1:nrow(re.count), function(g) {
		mod = lm(re.count[g,] ~ SampleInfo$Group)
		am = anova(mod)
		F1 = am[4][1, 1]   # F value
		P1 = am[5][1, 1]   # P value
		FP = cbind(F1, P1)
	})

	res = t(res)
	rownames(res) = rownames(re.count)
	colnames(res) = c("Fvalue", "Pvalue")	
	
	fdr = p.adjust(res[,2], method = "BH")
	means = sapply(Groups, function(g) apply(as.data.frame(re.count[,as.character(SampleInfo[as.character(SampleInfo$Group) %in% g, 1])]), 1, mean))
	colnames(means) = paste("Mean.In.", Groups, sep = "")
	
	result = data.frame(Taxon = rownames(re.count), means, res, FDR = fdr)
	loc = which(result$Taxon == "No_Rank")       
	if(length(loc) > 0)  result = result[-loc, ]
	loc1 = which(result$Taxon == "Unassigned")       
	if(length(loc1) > 0)  result = result[-loc1, ]
	result = result[order(result$Pvalue),]

	# 输出所有结果
	write.table(result,paste(taxonname, ".ANOVA_test_all_result.xls", sep = ""), row.names = F, quote = F, sep = "\t")

}else{

	stop("sample groups less than 2, please check it!\n")

}


if (Type == "P"){

	SigRes = subset(result, Pvalue < Threshold)

}else{

	SigRes = subset(result, FDR < Threshold)
	
}


cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
mycol= colors()[rep(cols,20)][1:length(Groups)]
index = sapply(1:nrow(SampleInfo),function(x) which(unique(SampleInfo$Group) %in% SampleInfo[x,2]))
# print(index)
color = mycol[index]

# 输出p值小于阈值（0.05）的结果并绘图
if(nrow(SigRes) > 0){

	write.table(SigRes,paste(taxonname, ".ANOVA_test_diff_P0.05_result.xls", sep = ""), row.names = F, quote = F, sep = "\t")
	##########plot##############
	# prepare plot data
	if(nrow(SigRes) == 1){

		DataForBoxplot = as.data.frame(t(re.count[rownames(SigRes),]*100))
		rownames(DataForBoxplot) = rownames(SigRes)

	}else{

		DataForBoxplot = re.count[rownames(SigRes),]*100

	}
	DataForBoxplot = DataForBoxplot[complete.cases(DataForBoxplot),]

	DataForBoxplot1 = data.frame(Taxon = rep(rownames(DataForBoxplot), ncol(DataForBoxplot)), 
							     Group = rep(SampleInfo[,2], each = nrow(DataForBoxplot)), 
							     Relative.Abundance = as.vector(as.matrix(DataForBoxplot)))


	Taxon = rownames(DataForBoxplot)

	### box_plot
	if(length(Taxon) == 1){

		pdf(paste(taxonname, ".taxon.DA.P0.05.Boxplot.pdf", sep = ""),height = 6,width = 6)

	}else if(length(Taxon) == 2){

		pdf(paste(taxonname, ".taxon.DA.P0.05.Boxplot.pdf", sep = ""),height = 6,width = 12)
		par(mar = c(8,6,6,2), mfrow = c(1, 2),cex.axis = 1,cex.main = 1.2,cex.lab= 1.2)

	}else{

		n = ceiling(length(Taxon) / 3)
		width = 15
		height = 5*n
		pdf(paste(taxonname, ".taxon.DA.P0.05.Boxplot.pdf", sep = ""),height = height,width = width)
		par(oma = c(0,0,0,0),mar = c(8,6,6,2), mfrow = c(n, 3),cex.axis = 1.5,cex.main = 2,cex.lab= 2)

	}

	for (i in 1:length(Taxon)){

		tax = Taxon[i]
		new = DataForBoxplot1[which(DataForBoxplot1$Taxon %in% tax), ] 
		if(any(SigRes$Taxon %in% tax)){
			
			p_value = (SigRes[which(SigRes$Taxon %in% tax), ])$Pvalue
			if(p_value >= 0.01){p = "*"}
			if(p_value < 0.01 & p_value >= 0.001){p = "**"}
			if(p_value < 0.001 & p_value >= 0.0001){p = "***"}
			if(p_value < 0.0001){p = "****"} 

		}else{
			p = ""
		}

		group = factor(new[,2], levels = unique(as.character(new[,2])))
		boxplot(new$Relative.Abundance ~ group, new, 
				lwd  = 2,
				xaxt = "n", 
				main = paste(tax,p),
				ylab = "Relative Abundance (%)",
				col = mycol)	

		text( x = c(1:length(Groups)), y = (par("usr")[3]-(par("usr")[4] - par("usr")[3])*0.02), adj = c(1,1),srt = 45, cex = 1.2, labels = levels(group), xpd = TRUE)
	}

	dev.off()

	### barplor

	# m = 10
	# if(length(Taxon) == 1){

		# pdf(paste(taxonname, ".taxon.DA.P0.05.Barplot.pdf", sep = ""),height = 6,width = m)

	# }else if(length(Taxon) == 2){

		# pdf(paste(taxonname, ".taxon.DA.P0.05.Barplot.pdf", sep = ""),height = 6,width = 2*m)
		# par(mar = c(8,6,6,2), mfrow = c(1, 2),cex.axis = 1,cex.main = 1.2,cex.lab= 1.2)

	# }else{

		# n = ceiling(length(Taxon) / 3)
		# width = 3*m
		# height = 5*n
		# pdf(paste(taxonname, ".taxon.DA.P0.05.Barplot.pdf", sep = ""),height = height,width = width)
		# par(oma = c(0,0,0,0),mar = c(8,6,6,2), mfrow = c(n, 3),cex.axis = 1.5,cex.main = 2,cex.lab= 2)

	# }

	# for (i in 1:length(Taxon)){

		# tax = Taxon[i]
		# new = DataForBoxplot1[which(DataForBoxplot1$Taxon %in% tax), ] 
		# if(any(SigRes$Taxon %in% tax)){
			
			# p_value = (SigRes[which(SigRes$Taxon %in% tax), ])$Pvalue
			# if(p_value >= 0.01){p = "*"}
			# if(p_value < 0.01 & p_value >= 0.001){p = "**"}
			# if(p_value < 0.001 & p_value >= 0.0001){p = "***"}
			# if(p_value < 0.0001){p = "****"} 

		# }else{
			# p = ""
		# }
		
		# #DataForBoxplot = as.matrix(DataForBoxplot)
		# barplot(height =  as.matrix(DataForBoxplot[which(rownames(DataForBoxplot) %in% tax), ]),
		        # border = NA,
		        # las  = 2,
		        # col  = color,
		        # main = paste(tax,p), 
		        # ylab = "Relative.Abundance (%)")

	# }

	# dev.off()
	
	
	####################################### heatmap #######################################
	####只有一个不画
	if (nrow(DataForBoxplot)>=2){
		if(!require(vegan)) {
		install.packages("vegan",repos='http://cran.us.r-project.org')
		}
		if(!require(RColorBrewer)) {
			install.packages("RColorBrewer",repos='https://cran.r-project.org')
		}
		if(!require(rnaseqWrapper)) {
			install.packages("rnaseqWrapper",repos='https://cran.r-project.org')
		}
		library(vegan)
		library(RColorBrewer)
		library(rnaseqWrapper)

		taxon2 = DataForBoxplot
		rownames(taxon2) = rownames(DataForBoxplot)
		colnames(taxon2) = colnames(DataForBoxplot)
		taxonNum = nrow(taxon2)
		sampleNum = ncol(taxon2)

		DataPlot2 = taxon2
		if( taxonNum < 120 ){
			cexRow = 1.7
		}else{
			cexRow = 1.7*(1/(ceiling(taxonNum/120)))
		}
		
		rownames(DataPlot2) = as.character(strsplit(rownames(DataPlot2), split="{.*}", perl=T) )

		if( max(nchar(rownames(DataPlot2))) > 23 ){
			malen = max(nchar(rownames(DataPlot2)))
		}else{
			malen = max(nchar(rownames(DataPlot2))) + 2
		}

		malen2 = max(nchar(colnames(DataPlot2)))

		width = 0.17*sampleNum + 0.15*malen
		height = 0.4*taxonNum + 0.15*malen2+3
			
		if(  width < 9 )   width = 9
		#if( height < 12 )  height = 12
		if( height > 26 )  height = 26

		# distance and clust function for heatmap.2
		BCdis = function(X){
			vegdist(X, method = "bray")
		}
		hBC = function(Y){
			hclust(Y, method = "average")
		}

		pdf(paste(taxonname, ".taxon.DA.Heatmap.pdf", sep = ""), width = width, height = height)
		layout(rbind(c(0,3,0), c(2,1,4), c(0,5,0)), heights = c(1.5,height-3.3,1.8), widths = c(1.5,width-(malen/11.5+1.5),malen/11.5), respect=FALSE)

		Qua2 = quantile(DataPlot2)
			
		if(length(table(DataPlot2)) < 400){
			breaks2 = unique(c(seq(0,Qua2[2],length=length(DataPlot2[which(DataPlot2<=Qua2[2])])),
				seq(Qua2[2],Qua2[3],length=length(unique(DataPlot2[which(DataPlot2>Qua2[2] & DataPlot2<=Qua2[3])]))+1),
				seq(Qua2[3],Qua2[4],length=length(unique(DataPlot2[which(DataPlot2>Qua2[3] & DataPlot2<=Qua2[4])]))+1),
				seq(Qua2[4],Qua2[5],length=length(unique(DataPlot2[which(DataPlot2>Qua2[4] & DataPlot2<=Qua2[5])]))+1)))
		
		}else{
			breaks2 = unique(c(seq(0,Qua2[2],length=90),
				seq(Qua2[2],Qua2[3],length=91),
				seq(Qua2[3],Qua2[4],length=91),
				seq(Qua2[4],Qua2[5],length=91)))
		}	

		color2=colorRampPalette(c("darkblue","darkgreen","yellow","darkred"))(length(breaks2)-1)

		heatmap.mark(DataPlot2, scale="none", col=color2, trace="none", cexCol=1.7, cexRow=cexRow, distfun=BCdis, hclustfun=hBC,
			breaks=breaks2, key=F, plotNew=FALSE, margins=c(malen2*7/10,malen/11.5))    # ۡ؝ظҪսҟާքߠk

		par(mar=c(4.5,0,3.5,malen/15))
		barplot(rep(0.05,length(breaks2)-1),space=0,border=color2,col=color2,xlab="Relative Abundance of community(%)",cex.lab=2,xpd=F,axes=F)

		axis(1,at=c(0,ceiling(length(breaks2)/3*1),ceiling(length(breaks2)/3*2),length(breaks2)-1),
				lab=c("0", as.character(round(breaks2[length(breaks2)/3*1],2)),
				as.character(round(breaks2[length(breaks2)/3*2],2)),
				as.character(round(breaks2[length(breaks2)],2))),
				cex.axis = 1.6)
		dev.off()
	}
	
}






	




