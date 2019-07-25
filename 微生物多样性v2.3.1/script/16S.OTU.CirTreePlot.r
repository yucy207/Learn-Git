# plot OTU phylogenetic Tree 
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
	"============== plot OTU phylogenetic Tree ===========
	Usage: 
	Rscript 16S.OTU.CirTreePlot.r otufile outputpath
		parameters ->
			otufile: [file -> always otu.tax.0.03.xls];
		     outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(ape)

otufile = normalizePath(args[1])
outputpath = normalizePath(args[2])

OTUFile = read.table(otufile, header = TRUE, row.names = 1, sep = "\t", check.name = F, comment.char = "", quote= "")

loc = grep("OTUsize", colnames(OTUFile))

if(loc > 0){

	otudata = as.matrix(OTUFile[,1:(loc-1)])
	rownames(otudata) = rownames(OTUFile)
	colnames(otudata) = colnames(OTUFile)[1:(loc-1)]

}

if(nrow(otudata) < 300){
  
	data_for_plot = otudata
  
}else{
  
	data_relative_abundance = sapply(1:ncol(otudata), function(x) otudata[,x]/sum(otudata[,x]))
  
	index = which(apply(data_relative_abundance > 0.01, 1, any) == T)
	data_for_plot = as.matrix(otudata[index,])
	colnames(data_for_plot) = colnames(otudata)
  
}

nch = max(nchar(rownames(data_for_plot)))

if ( nch < 10 ){
	width = 10; height = 10
}else if ( nch > 20 ){
	width = 20; height = 20
}else{
	width = nch; height = nch
}

filename = unlist(strsplit(otufile, "/|\\.xls"))[length(unlist(strsplit(otufile, "/|\\.xls")))]

pdf(paste(outputpath, paste(filename, "_CirTree.pdf", sep = ""), sep = "/"), width = width, height = height)

dd = dist(scale(data_for_plot), method = "euclidean")
hc = hclust(dd, method = "ward.D2")
#nclus= 3
#color=c('red','blue','green')
#color_list=rep(color,nclus/length(color))
#clus=cutree(hc,nclus)
plot(as.phylo(hc), type = "fan", label.offset = 0.1, cex = 5/nch)   #  tip.color=color_list[clus], 
dev.off()
