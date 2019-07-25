# 16s.heatmap Using level.*sort.xls files
rm(list = ls())
args = commandArgs(trailingOnly=T)

if ( length(args)!= 3 ){

	cat(
	"Usage: 
	Rscript 16S.Community.Heatmap.r taxonfile N outputpath
		parameters ->
			  taxonfile: [file -> always *.taxon.Abundance.xls with row taxa and column samples];
				  N: [number -> N most abundant taxa to plot heatmap]
			 outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()
	
}

# 判断包是否存在并加载
if(!require(vegan)) {
	install.packages("vegan",repos='http://cran.us.r-project.org')
}
if(!require(RColorBrewer)) {
	install.packages("RColorBrewer",repos='https://cran.r-project.org')
}
if(!require(rnaseqWrapper)) {
	install.packages("rnaseqWrapper",repos='https://cran.r-project.org')
}

taxonfile = normalizePath(args[1])
N = as.numeric(args[2])
outputpath = normalizePath(args[3])

inputdata = read.table(taxonfile, header = TRUE, row.names = 1, check.name = F, quote = "", comment.char = "", sep = "\t", fill = T)    # 注意 quote="" 

# 如果存在分类注释信息,则去除
loc = grep("OTUsize|Abundance", colnames(inputdata))
if(length(loc) > 0) taxon = as.matrix(inputdata[,1:(loc-1)])

if( ncol(taxon) == 1){

	stop("Sorry: You must have at least 2 rows and 2 columns to plot heatmap!\n")

}else{

	taxon = taxon[rowSums(taxon) != 0, ]

}

library(vegan)
library(RColorBrewer)
library(rnaseqWrapper)

rownames(taxon) = gsub("\"","",rownames(taxon))     # 去除名字中的双引号
colnames(taxon) = gsub("\"","",colnames(taxon))
rownames(taxon) = as.character(strsplit(rownames(taxon), split = "{.*}", perl = T))

index = which(rownames(taxon) == "All")       # 去除 All 行
if(length(index) > 0){

	taxon = as.matrix(taxon[-index, ])
	colnames(taxon) = colnames(inputdata)[1:(loc-1)]

}

if(nrow(taxon) > 1){
	
	taxon1 = sapply(1:ncol(taxon), function(t) taxon[,t]/sum(taxon[,t]))
	rownames(taxon1) = rownames(taxon)
	colnames(taxon1) = colnames(taxon)
	taxonNum = nrow(taxon1)
	sampleNum = ncol(taxon1)

	if( taxonNum < N ){
		DataPlot = taxon1	   # 小于N,全画
	}else{
		taxon1 = cbind(taxon1, apply(taxon1, 1, sum))
		taxon1 = taxon1[order(taxon1[,ncol(taxon1)], decreasing = T), ]
		DataPlot = taxon1[1:N, -ncol(taxon1)]           #  相对丰度最高的N个
	}

	rownames(DataPlot) = as.character( strsplit(rownames(DataPlot), split="{.*}", perl=T) )
	
	if( max(nchar(rownames(DataPlot))) > 23 ){
		malen = max(nchar(rownames(DataPlot)))
	}else{
		malen = max(nchar(rownames(DataPlot))) + 2
	}

	malen2 = max(nchar(colnames(DataPlot)))

	width = 0.17*sampleNum + 0.15*malen
	height = 0.2*taxonNum + 0.15*malen2
		
	if(  width < 9 )   width = 9
	if( height < 12 )  height = 12
	if( height > 26 )  height = 26

	# distance and clust function for heatmap.2
	BCdis = function(X){
		vegdist(X, method = "bray")
	}
	hBC = function(Y){
		hclust(Y, method = "average")
	}

	name = unlist(strsplit(taxonfile,"/"))
	names = strsplit(name[length(name)], "\\.\\w+\\.xls")

	pdf(paste(outputpath, paste(names, "Community.Heatmap.pdf", sep = "."), sep="/"), width = width, height = height)
	layout(rbind(c(0,3,0), c(2,1,4), c(0,5,0)), heights = c(1.5,height-3.3,1.8), widths = c(1.5,width-(malen/11.5+1.5),malen/11.5), respect=FALSE)
	# breaks
	Qua = quantile(DataPlot)
		
	if(length(table(DataPlot)) < 400){
		breaks = unique(c(seq(0,Qua[2],length=length(DataPlot[which(DataPlot<=Qua[2])])),
			seq(Qua[2],Qua[3],length=length(unique(DataPlot[which(DataPlot>Qua[2] & DataPlot<=Qua[3])]))+1),
			seq(Qua[3],Qua[4],length=length(unique(DataPlot[which(DataPlot>Qua[3] & DataPlot<=Qua[4])]))+1),
			seq(Qua[4],Qua[5],length=length(unique(DataPlot[which(DataPlot>Qua[4] & DataPlot<=Qua[5])]))+1)))
	
	}else{
		breaks=unique(c(seq(0,Qua[2],length=90),
			seq(Qua[2],Qua[3],length=91),
			seq(Qua[3],Qua[4],length=91),
			seq(Qua[4],Qua[5],length=91)))
	}	
	
	color=colorRampPalette(c("darkblue","darkgreen","yellow","darkred"))(length(breaks)-1)

	# heatmap.mark
	heatmap.mark(DataPlot, scale="none", col=color, trace="none", cexCol=1.7, cexRow=1.7, distfun=BCdis, hclustfun=hBC,
		breaks=breaks, key=F, plotNew=FALSE, margins=c(malen2*7/10,malen/11.5))    # 横纵坐标到边界的距离
	
	# color key
	par(mar=c(4.5,0,3.5,malen/15))
	barplot(rep(0.05,length(breaks)-1),space=0,border=color,col=color,xlab="Relative abundance of community (%)",cex.lab=2,xpd=F,axes=F)
	axis(1,at=c(0,ceiling(length(breaks)/3*1),ceiling(length(breaks)/3*2),length(breaks)-1),
			lab=c("0", as.character(round(breaks[length(breaks)/3*1]*100,2)),
			as.character(round(breaks[length(breaks)/3*2]*100,2)),
			as.character(round(breaks[length(breaks)]*100,2))),
			cex.axis = 1.6)
	dev.off()



	####-------------------------
	pdf(paste(outputpath, paste("uncluster",names, "Community.Heatmap.pdf", sep = "."), sep="/"), width = width, height = height)
	color=colorRampPalette(c("darkblue","darkgreen","yellow","darkred"))(length(breaks)-1)

	layout(rbind(c(0,3,0), c(2,1,4), c(0,5,0)), heights = c(1.5,height-3.3,1.8), widths = c(1.5,width-(malen/11.5+1.5),malen/11.5), respect=FALSE)
	# heatmap.mark
	heatmap.mark(DataPlot, scale="none", col=color, trace="none", cexCol=1.7, cexRow=1.7, distfun=BCdis, hclustfun=hBC, Colv = FALSE,
		breaks=breaks, key=F, plotNew=FALSE, margins=c(malen2*7/10,malen/11.5))    # 横纵坐标到边界的距离
	
	# color key
	par(mar=c(4.5,0,3.5,malen/15))
	barplot(rep(0.05,length(breaks)-1),space=0,border=color,col=color,xlab="Relative abundance of community (%)",cex.lab=2,xpd=F,axes=F)
	axis(1,at=c(0,ceiling(length(breaks)/3*1),ceiling(length(breaks)/3*2),length(breaks)-1),
			lab=c("0", as.character(round(breaks[length(breaks)/3*1]*100,2)),
			as.character(round(breaks[length(breaks)/3*2]*100,2)),
			as.character(round(breaks[length(breaks)]*100,2))),
			cex.axis = 1.6)
	dev.off()	


}
