rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== PCA ===========
	Usage: 
	Rscript 16S.BetaDiversity.PCA.r otufile groupfile outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(ade4)
library(scatterplot3d)

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

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

PCA_plot <- function(newSample,Group,mycol){

		OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
		otu_data = OTUFile[,rownames(newSample)]
		otu_data = otu_data[rowSums(otu_data) != 0,]
		loc = grep("All", rownames(otu_data))
		if(length(loc) > 0) otu_data = otu_data[-loc,]

		# PCA
		OTUtable = t(otu_data)
		num = apply(OTUtable,1,sum)
		#if (length(unique(num)) == 1){ 
		if (any(apply(OTUtable,2,function(g) length(unique(g))) == 1)){	

			otu.pca = prcomp(OTUtable) 

		}else{

			otu.pca = prcomp(OTUtable, scale = TRUE) 

		}

		# ====================== Plot =======================
		nm = max(nchar(Group))
		height = 8
		width = height + 0.1*nm

		# set color by group
		index = sapply(1:nrow(newSample),function(x) which(unique(newSample$Group)%in% newSample[x,2]))
		color = mycol[index]
		
		if ( ncol(otu.pca$x[,] ) >= 2 ) {

			x = otu.pca$x[,1]
			y = otu.pca$x[,2]
			sites = cbind(x, y)
			colnames(sites) = c("PC1", "PC2")
			write.csv(sites, "pc.sites.csv", quote=F)

			PC1.proportion = summary(otu.pca)$importance[2,1]*100      
			PC2.proportion = summary(otu.pca)$importance[2,2]*100  
			write.csv(summary(otu.pca)$importance, "pc.cont.csv", quote=F)

			pdf("PCA.pdf",width = width, height = height)
			par(oma=c(5,5,2,4))
			s.class(dfxy = sites,cgrid = 0, fac = factor(Group, level = unique(Group)), col = mycol, xax = 1, yax = 2 , cellipse=1 )
			axis(side = 1 ,line =5 )
			axis(side = 2 ,line =4 )
			title(main = "PCA", xlab = paste("PC1 [",PC1.proportion,"%]",sep = ""), ylab = paste("PC2 [",PC2.proportion,"%]", sep = ""), outer = TRUE)
			dev.off()
		}
		
		#3D PCA
		if ( ncol(otu.pca$x[,] ) >= 3 ) {

			PC1 = otu.pca$x[,1]
			PC2 = otu.pca$x[,2]
			PC3 = otu.pca$x[,3]
			sites = cbind(PC1, PC2, PC3)
			write.csv(sites, "3D.pc.sites.csv", quote=F)

			pdf("3D.PCA.pdf", width = 10, height = 10)
			layout(matrix(c(1,2)), heights = c(11, 3), respect = FALSE)    # 分面
			par(mar = c(0, 1, 1, 1))
			scatterplot3d(x = PC1, y = PC2, z = PC3, main = "PCA", xlab = expression(PC1), ylab = "", zlab = expression(PC3), pch = 19, color = color, cex.symbols = 1.2)   # cex.symbols means pont size
			
			dims <- par("usr")
			x <- dims[1]+ 0.9*diff(dims[1:2])
			y <- dims[3]+ 0.1*diff(dims[3:4])
			text(x,y,expression(PC2),srt=45)
			
			plot.new()
			legend("center", col = unique(color), pch = 19, xpd = T, cex = 1, ncol = 5, legend = unique(newSample$Group), bty = "n")
			dev.off()
		}
}


if(length(Info) > 5){

	Group = as.character(SampleInfo[,2])
	PCA_plot(SampleInfo,Group,mycol)

}else{

	for(m in 2:(length(Info))){

		g = combn(Info, m)	

		for (n in 1:ncol(g)){

			Sample = g[,n]
			newSample = SampleInfo[which(SampleInfo$Group %in% Sample), ]
			newGroup = as.character(newSample[,2])
			newmycol = group2col$col[group2col$gro %in% unique(newGroup)]

			#####创建差异分析的文件夹###########
			dir = paste(outputpath, paste(unique(newGroup), collapse = "_vs_"),sep = "/")
			system(paste("mkdir -p ", dir, sep = ""))
			setwd(dir)
			PCA_plot(newSample,newGroup,newmycol)

		}
	}
}
