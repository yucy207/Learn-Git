rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	"============== Use a variety of distance methods to do PCoA ===========
	Usage: 
	Rscript 16S.BetaDiversity.PCoA.r otufile trefile groupfile outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
			  trefile: [file -> always subsample_otu.repseq.fasta.tre come from 16S.Resample.pl];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(vegan)
library(phyloseq)   
library(ggrepel)
library(ape)
library(pheatmap)
library(ade4)
library(scatterplot3d)

otufile = normalizePath(args[1])
trefile = normalizePath(args[2])
groupfile = normalizePath(args[3])
outputpath = normalizePath(args[4])

setwd(outputpath)
# source("/home/zhengy/bin/modules/script/PCoA.r")
# ====================== prepare data =======================
SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Info = unique(as.character(SampleInfo[[2]]))
Groups = as.character(SampleInfo[,2])

# set color by group
cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
mycol= colors()[rep(cols,20)][1:length(unique(as.character(SampleInfo[[2]])))]
group2col = data.frame(unique(as.character(SampleInfo[[2]])),mycol,stringsAsFactors=FALSE)
colnames(group2col) = c("gro","col")

PERMANOVA_analysis <- function(newSample,Group,mycol){

	OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
	otu_data = OTUFile[,rownames(newSample)]
	otu_data = otu_data[rowSums(otu_data) != 0,]
	loc = grep("All", rownames(otu_data))
	if(length(loc) > 0) otu_data = otu_data[-loc,]

	OTU = otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)

	SAMPLE = sample_data(newSample)
	sample_names(SAMPLE) = as.character(rownames(newSample))

	physeq = phyloseq(OTU, SAMPLE)

	# Add phy_tree data
	#tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
	tree = read.tree(trefile)
	physeq = merge_phyloseq(physeq, tree)
	
	index = sapply(1:nrow(newSample),function(x) which(unique(newSample$Group) %in% newSample[x,2]))
	color = mycol[index]

	# ====================== Plot =======================
	nm = max(nchar(Group))
	height = 8
	width = height + 0.1*nm

	distance = c("bray", "wunifrac", "unifrac", "jaccard")

	for (d in distance){
  
		if (d == "jaccard"){
			Dist = distance(physeq, method = "jaccard", binary = TRUE)
		}else{
			Dist = distance(physeq, method = d)
		}

		pa = adonis(Dist~Group, permutations = 9999)
		pa$aov.tab[is.na(pa$aov.tab)] = " "
		data = cbind(rownames(pa$aov.tab),pa$aov.tab)
		colnames(data)[1] = ""	

		filename = paste(d, ".ADONIS.xls", sep = "")
		write.table(data, filename, sep = "\t", row.names = F, quote = F)



		distIngroups = lapply(unique(Group), function(g) {A = as.matrix(Dist)[which(Group == g), which(Group == g)]; dis = unique(A[A!=0])})
		names(distIngroups) = unique(Group)

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
		print(head(DataForPlot))

		# ====================== Plot =======================
		nm = max(nchar(Group))
		height = 8
		width = height + 0.3*nm

		# boxplot

		pdf(paste(d,".Dist.Boxplot.pdf", sep = ""), width = width, height = height)
		par(mar=c(0.3*nm, 0.2*nm, 1, 1), mgp=c(2, 0.7, 0))
		p = ggplot(DataForPlot, aes(x=Group, y=dis, fill=Group)) + geom_boxplot(position=position_dodge(0.8)) +  # 箱式图, position,组间距
			scale_fill_manual(values = mycol) +
			geom_jitter(position = position_dodge(0.8)) +     # Add points  http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
			theme_bw() +       # 白色背景
			theme(axis.text.x = element_text(angle=45, hjust=1)) +    # x轴标签倾斜
			theme(axis.title = element_text(size = 15), axis.text = element_text(size = 12)) +     # 坐标轴 title 大小   
			ylab("Distance") + scale_y_continuous(limits = c(min(DataForPlot$dis), max(DataForPlot$dis)*1.08)) +
			annotate("text", y = max(DataForPlot$dis)*1.06, x = levels(DataForPlot$Group)[length(levels(DataForPlot$Group))], label = paste("P =", pa$aov.tab$Pr[1]), hjust = 0.8, size = 5) 
		print(p)
		dev.off()
	
	}
}




# if (length(unique(Groups)) >= 3){

# 	dir1 = paste(outputpath, "All_groups", sep = "/")
# 	system(paste("mkdir -p ", dir1, sep = ""))
# 	setwd(dir1)
# 	PERMANOVA_analysis(otu_data,Groups,mycol)

# 	######pair groups#######
# 	group = unique(SampleInfo$Group)

# 	for(m in 1:(length(group)-1)){

# 		for (n in (m+1):length(group)){

# 			#####创建两两差异分析的文件夹###########
# 			dir2 <- paste(outputpath, paste(group[m], "_vs_", group[n], sep = ""), sep = "/")
# 			system(paste("mkdir -p ", dir2, sep = ""))
# 			setwd(dir2)
# 			#####创建新的分组信息和otu表###################           
# 			SampleInfo1 <- rbind(SampleInfo[SampleInfo$Group == group[m],], SampleInfo[SampleInfo$Group == group[n],])
# 			otu_data_1 = OTUFile[,rownames(SampleInfo1)]
# 			Groups_1 = as.character(SampleInfo1[,2])
# 			mycol_1 = group2col$col[group2col$gro %in% unique(Groups_1)]
# 			PERMANOVA_analysis(otu_data_1,Groups_1,mycol_1)

# 		}

# 	}

# }else{

# 	dir1 = paste(outputpath, "All_groups", sep = "/")
# 	system(paste("mkdir -p ", dir1, sep = ""))
# 	setwd(dir1)
# 	PERMANOVA_analysis(otu_data,Groups,mycol)

# }



if(length(Info) > 5){

	Group = as.character(SampleInfo[,2])
	PERMANOVA_analysis(SampleInfo,Group,mycol)

}else{

	for(m in 2:(length(Info))){

		g = combn(Info, m)	

		for (n in 1:ncol(g)){

			Sample = g[,n]
			newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
			newGroup = as.character(newSample[,2])
			newmycol = group2col$col[group2col$gro %in% unique(newGroup)]
			
			if(length(newSample[[1]]) ==2) next
		
			#####创建差异分析的文件夹###########
			dir = paste(outputpath, paste(unique(newGroup), collapse = "_vs_"),sep = "/")
			system(paste("mkdir -p ", dir, sep = ""))
			setwd(dir)
			PERMANOVA_analysis(newSample,newGroup,newmycol)
		}
	}
}









