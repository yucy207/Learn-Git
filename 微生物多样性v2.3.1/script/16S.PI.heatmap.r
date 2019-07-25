rm(list = ls())
args = commandArgs(trailingOnly=T)

if ( length(args)!= 4 ){

	cat(
	"Usage: 
	Rscript 16S.PI.Heatmap.r genefile groupfile N outputpath
		parameters ->
			  genefile: [file -> always * for_stamp.txt with row gene and column samples];
			  groupfile: [];
				  N: [number -> N most abundant gene to plot heatmap]
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()
	
}


library(ComplexHeatmap)
require(circlize)


genefile = normalizePath(args[1])
groupfile = normalizePath(args[2])
N = as.numeric(args[3])
outputpath = normalizePath(args[4])

inputdata = read.table(genefile, header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T) 

SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]

gene = inputdata[rowSums(inputdata)!=0, rownames(SampleInfo)]


geneNum = nrow(gene)
sampleNum = ncol(gene)

if( geneNum < N ){
	DataPlot = gene	   # 小于N,全画
}else{
	gene = cbind(gene, apply(gene, 1, sum))

	gene = gene[order(gene[,ncol(gene)], decreasing = T), ]
	DataPlot = gene[1:N, -ncol(gene)]           #  相对丰度最高的N个
}

geneNum1 = nrow(DataPlot)


DataPlot1 = sapply(1:ncol(DataPlot), function(t) DataPlot[,t]*100/sum(DataPlot[,t]))
rownames(DataPlot1) = rownames(DataPlot)
colnames(DataPlot1) = colnames(DataPlot)


if( max(nchar(rownames(DataPlot1))) > 23 ){
	malen = max(nchar(rownames(DataPlot1)))
}else{
	malen = max(nchar(rownames(DataPlot1))) + 2
}
malen2 = max(nchar(colnames(DataPlot1)))

width =  1  + 0.18*sampleNum + 0.11*malen
height = 1  + 0.21*geneNum1 + 0.02*malen2
height1 = 0.8  + 0.21*geneNum1 + 0.02*malen2


print(width)
print(height)

width = ifelse(width < 9, 9, width)

if(width >= 20){
	wd1 = 0.6*width
	legend_width = 3
	unitx = -0.42	

}else if(width > 15 && width < 20){
	wd1 = 0.6*width
	legend_width = 3
	unitx = -0.44

}else if(width <= 15 && width > 9){
	wd1 = 0.35*width
	legend_width = 3
	unitx = -0.4

}else{

	wd1 = 0.25*width
	legend_width = 3
	unitx = -0.4
}


name = unlist(strsplit(genefile,"/"))
names = strsplit(name[length(name)], "\\_\\w+\\.txt")

pdf(paste(outputpath, paste(names, "heatmap.pdf", sep = "."), sep="/"), width = width, height = height)
Heatmap(DataPlot1, 
	col = colorRamp2(seq(floor(min(DataPlot1)), ceiling(max(DataPlot1)), length=4), c("white", "#3CB1DA","#0055A4", "#00296B")),
	rect_gp = gpar(col = "gray70"), 
	row_names_gp = gpar(fontsize = 10),
	column_names_gp = gpar(fontsize = 10),
	show_heatmap_legend = FALSE,
	width = unit(wd1,"inches"))

col_fun = colorRamp2(seq(floor(min(DataPlot1)), ceiling(max(DataPlot1)), length=4), c("white", "#3CB1DA","#0055A4", "#00296B"))
col_lgd = Legend(at = round(seq(floor(min(DataPlot1)), ceiling(max(DataPlot1)), length=3),2), col_fun = col_fun,
				title = "proportion (%)", title_position="topcenter", direction = "horizontal",
				border = "black", legend_width = unit(legend_width, "cm"), 
				title_gp = gpar(fontsize = 10))

pushViewport(viewport(x = unit(unitx, "npc"), y = unit(0.4, "npc"), just = c("left", "bottom")))
grid.draw(col_lgd)
popViewport()
dev.off()



pdf(paste(outputpath, paste("uncluster", names, "heatmap.pdf", sep = "."), sep="/"), width = width, height = height1)
Heatmap(DataPlot1, 
	col = colorRamp2(seq(floor(min(DataPlot1)), ceiling(max(DataPlot1)), length=4), c("white", "#3CB1DA","#0055A4", "#00296B")),
	rect_gp = gpar(col = "gray70"), 
	row_names_gp = gpar(fontsize = 10),
	column_names_gp = gpar(fontsize = 10),
	cluster_columns = FALSE,
	show_heatmap_legend = FALSE,
	width = unit(wd1,"inches"))

col_fun = colorRamp2(seq(floor(min(DataPlot1)), ceiling(max(DataPlot1)), length=4), c("white", "#3CB1DA","#0055A4", "#00296B"))
col_lgd = Legend(at = round(seq(floor(min(DataPlot1)), ceiling(max(DataPlot1)), length=3),2), col_fun = col_fun,
				title = "proportion (%)", title_position="topcenter", direction = "horizontal",
				border = "black", legend_width = unit(legend_width, "cm"), 
				title_gp = gpar(fontsize = 10))

pushViewport(viewport(x = unit(unitx, "npc"), y = unit(0.4, "npc"), just = c("left", "bottom")))
grid.draw(col_lgd)
popViewport()

dev.off()
