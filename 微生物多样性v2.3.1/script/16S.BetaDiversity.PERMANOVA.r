rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat("Usage: 
	Rscript 16S.BetaDiversity.PERMANOVA.r otufile groupfile outputpath
		parameters ->
			    otufile: [file -> always otu.tax.0.03.xls];
			  groupfile: [file -> always sample.groups without header and only two column];
		         outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

library(vegan)
library(ggplot2)

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

# ====================== prepare data =======================
SampleInfo = read.table(groupfile, header = F, sep="\t", colClasses=c('factor', 'factor'))
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Groups = as.character(SampleInfo[,2])

OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
otu_data = OTUFile[,rownames(SampleInfo)]

loc = grep("All", rownames(otu_data)) 
if (length(loc) > 0) otu_data = otu_data[-loc, ]

cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
mycol= colors()[rep(cols,20)][1:length(levels(SampleInfo[,2]))]
group2col = data.frame(as.character(levels(SampleInfo[,2])),mycol,stringsAsFactors=FALSE)
colnames(group2col) = c("gro","col")

PERMANOVA_analysis <- function(otu,groups,col){

	DistBC = vegdist(t(otu), method = "bray")
	pa = adonis(DistBC~groups, permutations = 9999)
	pa$aov.tab[is.na(pa$aov.tab)] = " "
	data = cbind(rownames(pa$aov.tab),pa$aov.tab)
	colnames(data)[1] = ""
	#write.csv(pa$aov.tab, "ADONIS.csv", quote= FALSE)
	write.table(data, "ADONIS.xls", sep = "\t", row.names = F, quote = F)

	distIngroups = lapply(unique(groups), function(g) {A = as.matrix(DistBC)[which(groups == g), which(groups == g)]; dis = unique(A[A!=0])})
	names(distIngroups) = unique(groups)

	for (i in 1:length(distIngroups)){

		if (i == 1){
			res = as.matrix(distIngroups[[i]])
			rownames(res) = rep(names(distIngroups)[i], nrow(res))
		}else{
	
			tmp = as.matrix(distIngroups[[i]])
			rownames(tmp) = rep(names(distIngroups)[i], nrow(tmp))
			res = rbind(res, tmp)
		}
	}

	DataForPlot = data.frame(Group = rownames(res), dis = res[,1])

	# ====================== Plot =======================
	nm = max(nchar(groups))
	height = 8
	width = height + 0.3*nm

	# boxplot
	pdf("Dist.Boxplot.pdf", width = width, height = height)
	par(mar=c(0.3*nm, 0.2*nm, 1, 1), mgp=c(2, 0.7, 0))
	p = ggplot(DataForPlot, aes(x=Group, y=dis, fill=Group)) + geom_boxplot(position=position_dodge(0.8)) +  # 箱式图, position,组间距
		scale_fill_manual(values = col) +
		geom_jitter(position = position_dodge(0.8)) +     # Add points  http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
		theme_bw() +       # 白色背景
		theme(axis.text.x = element_text(angle=45, hjust=1)) +    # x轴标签倾斜
		theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) +     # 坐标轴 title 大小   
		ylab("Distance") + scale_y_continuous(limits = c(min(DataForPlot$dis), max(DataForPlot$dis)*1.08)) +
		annotate("text", y = max(DataForPlot$dis)*1.06, x = levels(DataForPlot$Group)[length(levels(DataForPlot$Group))], label = paste("P =", pa$aov.tab$Pr[1]), hjust = 0.8, size = 5) 
	print(p)
	dev.off()

}


######all groups#######
if (length(unique(Groups)) >= 3){

	dir1 = paste(outputpath, "All_groups", sep = "/")
	system(paste("mkdir -p ", dir1, sep = ""))
	setwd(dir1)
	PERMANOVA_analysis(otu_data,Groups,mycol)

	######pair groups#######
	group = unique(SampleInfo$Group)

	for(m in 1:(length(group)-1)){

		for (n in (m+1):length(group)){

			#####创建两两差异分析的文件夹###########
			dir2 <- paste(outputpath, paste(group[m], "_vs_", group[n], sep = ""), sep = "/")
			system(paste("mkdir -p ", dir2, sep = ""))
			setwd(dir2)
			#####创建新的分组信息和otu表###################           
			SampleInfo1 <- rbind(SampleInfo[SampleInfo$Group == group[m],], SampleInfo[SampleInfo$Group == group[n],])
			otu_data_1 = OTUFile[,rownames(SampleInfo1)]
			Groups_1 = as.character(SampleInfo1[,2])
			mycol_1 = group2col$col[group2col$gro %in% unique(Groups_1)]
			PERMANOVA_analysis(otu_data_1,Groups_1,mycol_1)

		}

	}

}else{

	dir1 = paste(outputpath, "All_groups", sep = "/")
	system(paste("mkdir -p ", dir1, sep = ""))
	setwd(dir1)
	PERMANOVA_analysis(otu_data,Groups,mycol)

}







