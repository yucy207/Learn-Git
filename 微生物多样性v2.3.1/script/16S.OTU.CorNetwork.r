# co-expression network plot using igraph
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

    cat(
    "============== do OTU cor ===========
    Usage: 
    Rscript 16S.OTU.CorNetwork.r inputfile1 inputfile2 outputpath
        parameters ->
              inputfile1: [file -> always OTU.Correlation.*.xls];
              inputfile2: [file -> always OTU.size.*.for.Cytoscape.classif.xls];
              outputpath: [path -> path for output]; \n")
    options("show.error.messages" = F) 
    stop()

}
inputfile1 = normalizePath(args[1])
#inputfile2 = normalizePath(args[2])
outputpath = normalizePath(args[3])
setwd(outputpath)

library(igraph)

inputpath = normalizePath(args[2])
file = list.files(inputpath,pattern = "*.Cytoscape.classif.txt$")
for (i in 1:length(file)){
	inputfile2 = normalizePath(file[i])


	filename = unlist(strsplit(inputfile2, "/|\\.xls"))[length(unlist(strsplit(inputfile2, "/|\\.xls")))]

	cordata = read.table(inputfile1, header=T, sep="")
	classidata = read.table(inputfile2, header=T, sep="")

	net = graph_from_data_frame(d=cordata, vertices=classidata, directed = F)

	# ================================= Plot =================================
	# plot network using igraph package

	mycol = c(119,132,147,454,89,404,123,529,463,552,28,54,84,100,558,43,652,31,610,477,256,588,99,81,503,104,562,76,96,495)
	mycol = colors()[rep(mycol,20)]
	## Node color
	color = mycol[ 1 : length(unique(classidata[,3])) ] 

	colors = rep(0, length(classidata[,3]))
	for (i in 1:length(color)){
		
		colors[grep(unique(classidata[,3])[i], as.character(classidata[,3]))] = color[i]

	}

	## Edge_color
	E(net)[cordata[[5]] == 1]$color = "red"
	E(net)[cordata[[5]] == -1]$color = "green"


	width = 10
	mr = max(nchar(unique(as.character(classidata[,3]))))     # 分类标签最长字符数	
	ncolu = ceiling(width*8/mr)                              # 根据图形宽度，限制标签列数
	h1 = ceiling(length(unique(classidata[,3]))/ncolu)        # 标签行数 
	height = width+0.5*h1                                         # 根据标签行数和宽度限定高度

	# d = degree(net, mode = "all")
	# d = d[order(d,decreasing = TRUE)][50]
	# bad.vs = V(net)[degree(net) < d]          # remove isolated nodes，即去掉和所有otu均无相关性的otu 
	# net = delete.vertices(net, bad.vs)

	deg = degree(net, mode = "all")

	pdf(paste(filename, ".pdf", sep = ""), width = width, height = height)
	# pdf(paste(filename, ".pdf", sep = ""))
	layout(matrix(c(1,2), 2, 1), heights=c(width, 0.5*h1))
	par(mar = c(0.1,0.1,0.1,0.1))
	plot(net,    vertex.color = colors,        # Node color
				  vertex.size = deg,           # Node size
		   vertex.frame.color = NA,            # Node border color
				 vertex.label = NA,            # Node label
				   edge.color = E(net)$color,  # edge color
				   edge.label = NA,            # edge label
				  edge.curved = 0.001,
					 # layout = layout.circle   
					 # layout = layout.graphopt  
					 # layout = layout.kamada.kawai 
					 # layout = layout.graphopt
					   layout = layout_nicely
		)  
	plot.new()
	par(mar = c(0.1,0.2,0.2,0.1), xpd = TRUE)
	legend("center", legend = unique(classidata[,3]), fill = color, ncol = ncolu, bty = "n", border = NA, cex = 0.8)
	dev.off()
}