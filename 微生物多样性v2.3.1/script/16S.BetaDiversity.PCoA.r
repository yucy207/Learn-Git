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
source("/home/zhengy/bin/modules/script/PCoA.r")
# ====================== prepare data =======================
SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Info = unique(as.character(SampleInfo[[2]]))

# set color by group
cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
mycol= colors()[rep(cols,20)][1:length(unique(as.character(SampleInfo[[2]])))]
group2col = data.frame(unique(as.character(SampleInfo[[2]])),mycol,stringsAsFactors=FALSE)
colnames(group2col) = c("gro","col")

PCoA_plot <- function(newSample,Group,mycol){

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
		
		write.csv(as.matrix(Dist), paste(d, ".Dist.csv", sep = ""), quote=F)     # export dissimilarity matrix
		
		# heatmap
		pdf(paste(d, ".dist.Heatmap.pdf", sep = ""), width = width, height = height, onefile = FALSE) 
		pheatmap(as.matrix(Dist), show_rownames = T)
		dev.off()
		
		ord = ordinate(physeq, method = "PCoA", distance = Dist)  #(pcoa(Dist)) from  library(ape)
		#ord = PCoA(Dist)  # from source("/home/zhengy/bin/modules/script/PCoA.r")
		
		if (ncol(ord$vectors) >= 2){
			# PCoA
			position = ord$vectors[,1:2]
			write.csv(position, paste(d, ".PCoA.site.csv",sep = ""), quote=F)

			x = position[,1]
			xt = (max(x)-min(x))/4
			y = position[,2]
			yt = (max(y)-min(y))/4

			pdf(paste(d, ".PCoA.pdf", sep = ""), width = width, height = height)
			par(oma=c(5,5,2,4))
			s.class(dfxy = position ,cgrid = 0, fac = factor(Group, level = unique(Group)), col = mycol, xax = 1, yax = 2 , cellipse=1, ylim = c(min(y)-yt, max(y)+yt), xlim = c(min(x)-xt, max(x)+xt))
			axis(side = 1 ,line =5 )
			axis(side = 2 ,line =4 )
			title(main = "PCoA", xlab = paste("Axis.1  [",round(ord$values[1,2]*100,2),"%]",sep = ""), ylab = paste("Axis.2  [",round(ord$values[2,2]*100,2),"%]",sep = ""), outer = TRUE)
			dev.off()
		}
		
		if (ncol(ord$vectors) >= 3){
			# 3D.PCoA
			position = ord$vectors[,1:3]
			write.csv(position, paste(d, ".3D.PCoA.site.csv",sep = ""), quote=F)
			
			Axis.1 = position[,1]
			Axis.2 = position[,2]
			Axis.3 = position[,3]

			pdf(paste(d, ".3D.PCoA.pdf", sep = ""), width = 8, height = 9)
			layout(matrix(c(1,2)), heights = c(11, 3), respect = FALSE)    # 分面
			par(mar = c(0, 1, 1, 1))
			scatterplot3d(x = Axis.1, y = Axis.2, z = Axis.3, main = "PCoA", xlab=expression(Axis.1), ylab="", zlab=expression(Axis.3), pch = 19, color = color, cex.symbols = 1.2)   # cex.symbols means pont size
			dims <- par("usr")
			x <- dims[1]+ 0.9*diff(dims[1:2])
			y <- dims[3]+ 0.08*diff(dims[3:4])
			text(x,y,expression(Axis.2),srt=45)
			
			plot.new()
			legend("center", col = unique(color), pch = 19, xpd = T, cex = 1, ncol = 5, legend = unique(newSample$Group), bty = "n")
			dev.off()
		}
	}
}

if(length(Info) > 5){

	Group = as.character(SampleInfo[,2])
	PCoA_plot(SampleInfo,Group,mycol)

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
			PCoA_plot(newSample,newGroup,newmycol)
		}
	}
}









