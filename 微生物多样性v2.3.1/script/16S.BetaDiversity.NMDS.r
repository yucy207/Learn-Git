rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	"============== Use a variety of distance methods to do NMDS ===========
	Usage: 
	Rscript 16S.BetaDiversity.NMDS.r otufile trefile groupfile outputpath
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
library(vegan)
library(ape)
library(ade4)
library(scatterplot3d)

otufile = normalizePath(args[1])
trefile = normalizePath(args[2])
groupfile = normalizePath(args[3])
outputpath = normalizePath(args[4])

setwd(outputpath)

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


NMDS_plot <- function(newSample,Group,mycol){

	OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
	otu_data = OTUFile[,rownames(newSample)]
	otu_data = otu_data[rowSums(otu_data) != 0,]
	loc = grep("All", rownames(otu_data))
	if(length(loc) > 0) otu_data = otu_data[-loc,]

	OTU = otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)
	SAMPLE = sample_data(newSample)
	sample_names(SAMPLE) = as.character(rownames(newSample))
	physeq = phyloseq(OTU, SAMPLE)

	# prune OTUs that are not present in at least one sample
	physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

	# Add phy_tree data
	#tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
	tree = read.tree(trefile)
	physeq = merge_phyloseq(physeq, tree)

	index = sapply(1:nrow(newSample),function(x) which(unique(newSample$Group) %in% newSample[x,2]))
	# print(index)
	color = mycol[index]
	# print(index)
	# print(color)
	# cat("==============================================")
	# ====================== Plot =======================
	nm = max(nchar(Group))
	height = 8
	width = height + 0.1*nm

	# NMDS
	distance = c("bray", "wunifrac", "unifrac", "jaccard")

	for (d in distance){

		if (d == "jaccard"){
			Dist = distance(physeq, method = "jaccard", binary = TRUE)
		}else{
			Dist = distance(physeq, method = d)
		}
		
		# NMDS
		otu.mds = metaMDS(Dist)
		# print(newSample)
		# print(d)
		nmdScore = scores(otu.mds, choices = c(1,2))
		stress = otu.mds$stress
		write.csv(nmdScore, paste(d, ".NMDS.site.csv",sep = ""), quote=F)        

		x = otu.mds$points[,1]
		xt = (max(x)-min(x))/4
		y = otu.mds$points[,2]
		yt = (max(y)-min(y))/4
		sites = cbind(x, y)

		pdf(paste(d, ".NMDS.pdf", sep = ""), width = width, height = height)
		par(oma=c(5,5,2,4))
		# s.class(dfxy = sites ,cgrid = 0,fac = factor(Group, level = unique(Group)), col = mycol, xax = 1, yax = 2 , cellipse=1 )
		s.class(dfxy = sites ,cgrid = 0,fac = factor(Group, level = unique(Group)), col = mycol, xax = 1, yax = 2 , cellipse=1 , ylim = c(min(y)-yt, max(y)+yt), xlim = c(min(x)-xt, max(x)+xt))
		axis(side = 1 ,line =5 )
		axis(side = 2 ,line =4 )
		title(main = "NMDS", xlab = paste("MDS1\n","STRESS: ",round(stress,2),sep=""), ylab = "MDS2", outer = TRUE)
		dev.off()
		
		# 3D.NMDS
		otu.mds=metaMDS(Dist, k = 3)
		nmdScore = scores(otu.mds, choices = c(1,2,3))
		write.csv(nmdScore, paste(d, ".3D.NMDS.site.csv",sep = ""), quote=F)        

		MDS1 = otu.mds$points[,1]
		MDS2 = otu.mds$points[,2]
		MDS3 = otu.mds$points[,3]

		pdf(paste(d, ".3D.NMDS.pdf", sep = ""), width = 8, height = 9)
		layout(matrix(c(1,2)), heights = c(11, 3), respect = FALSE)    # 分面
		par(mar = c(0, 1, 1, 1))
		scatterplot3d(x = MDS1, y = MDS2, z = MDS3, main = "NMDS", xlab = expression(MDS1), ylab = "", zlab = expression(MDS3), pch = 19, color = color, cex.symbols = 1.2)   # cex.symbols means pont size
		
		dims <- par("usr")
		x <- dims[1]+ 0.9*diff(dims[1:2])
		y <- dims[3]+ 0.1*diff(dims[3:4])
		text(x,y,expression(MDS2),srt=45)
		
		plot.new()
		legend("center", col = unique(color), pch = 19, xpd = T,cex = 1, ncol = 5, legend = unique(newSample$Group), bty = "n")
		dev.off()

	}
	
}

if(length(Info) > 5){

	Group = as.character(SampleInfo[,2])
	NMDS_plot(SampleInfo,Group,mycol)

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
			NMDS_plot(newSample,newGroup,newmycol)
			# print(newGroup)
			# print(newmycol)

		}
	}
}










