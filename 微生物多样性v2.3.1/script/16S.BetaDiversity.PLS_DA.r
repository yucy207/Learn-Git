#  http://blog.sina.com.cn/s/blog_942438cf0102wth0.html
# http://mixomics.org/graphics/sample-plot/plotindiv/

rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== PLS_DA ===========
	Usage: 
	Rscript 16S.PLS_DA.r otufile groupfile outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(mixOmics)
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
Groups = SampleInfo[,2]
Info = unique(as.character(SampleInfo[[2]]))

# set color by group
cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
mycol= colors()[rep(cols,20)][1:length(unique(as.character(SampleInfo[[2]])))]
group2col = data.frame(unique(as.character(SampleInfo[[2]])),mycol,stringsAsFactors=FALSE)
colnames(group2col) = c("gro","col")


PLS_DA_plot <- function(newSample,Group,mycol) {

	OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
	otu_data = OTUFile[,rownames(newSample)]
	otu_data = otu_data[rowSums(otu_data) != 0,]
	loc = grep("All", rownames(otu_data))
	if(length(loc) > 0) otu_data = otu_data[-loc,]


	plsda = plsda(t(otu_data), Group, ncomp = 2)
	x = plsda$variates$X[,1]
	y = plsda$variates$X[,2]
	sites = cbind(x, y)
	colnames(sites) = c("X-variate 1", "X-variate 2")
	write.csv(sites, "plsda.sites.csv", quote=F)	

	PC1.proportion = round(plsda$explained_variance$X[1]*100, 2)      
	PC2.proportion = round(plsda$explained_variance$X[2]*100, 2)

	nm = max(nchar(colnames(otu_data)))
	height = 8
	width = height + 0.3*nm
	pdf("PLS-DA.pdf",width = width, height = height)
	par(oma=c(5,5,2,4))
	s.class(dfxy = sites,cgrid = 0, fac = factor(Group, level = unique(Group)), col = mycol, xax = 1, yax = 2 , cellipse=1 )
	axis(side = 1 ,line =5 )
	axis(side = 2 ,line =4 )
	title(main = "PLS_DA", xlab = paste("X-variate 1: ",PC1.proportion,"% expl.var",sep = ""), ylab = paste("X-variate 2: ",PC2.proportion,"% expl.var", sep = ""), outer = TRUE)
	dev.off()

	# nm = max(nchar(colnames(otu_data)))
	# height = 8
	# width = height + 0.3*nm
	# pdf("PLS-DA.pdf", width = width, height = height)
	# # print(plsda)
	# print(plsda)
	# plotIndiv(plsda, ellipse = T, plot.ellipse = T, legend = T, ind.names = F, title = "PLS-DA", group = Group, pch = 19, col = mycol)
	# dev.off()

}


if(length(Info) > 5){

	Group = as.character(SampleInfo[,2])
	PLS_DA_plot(SampleInfo,Group,mycol)

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
			PLS_DA_plot(newSample,newGroup,newmycol)

		}
	}
}


