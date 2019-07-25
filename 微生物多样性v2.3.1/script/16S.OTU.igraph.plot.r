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

library(igraph)
library(RColorBrewer)

edges_input = normalizePath(args[1])
nodes_input = normalizePath(args[2])
outputpath  = normalizePath(args[3])
setwd(outputpath)
	

filename = unlist(strsplit(nodes_input, "/|\\.xls"))[length(unlist(strsplit(nodes_input, "/|\\.xls")))]

links <- read.table(edges_input, header = T, as.is = T)
nodes <- read.table(nodes_input, header = T, as.is = T)
net   <- graph_from_data_frame(d = links,vertices = nodes[[1]], directed = F)

V(net)$label <- nodes[[3]]
#V(net)$size  <- nodes[[2]]

#Size
min_source_size <- min(nodes[[2]])
max_source_size <- max(nodes[[2]])

min_target_size <- 3
max_target_size <- 10
V(net)$size  <- (nodes[[2]] - min_source_size )  / (max_source_size - min_source_size) * (max_target_size - min_target_size) + min_target_size

#Edge_color
E(net)[links[[5]] == 1]$color <- "red"
E(net)[links[[5]] == -1]$color <- "green"

#Node_color
cols   = brewer.pal(length(unique(nodes[[3]])),"Set1")
name = unique(nodes[[3]])
names(cols) = name

V(net)$color <- unlist(lapply(nodes[[3]], function(t){cols[[t]]}))

#plot
pdf(paste(filename, ".pdf", sep = ""))
plot(net, vertex.frame.color = NA,
				vertex.label = NA,
				 vertex.size = V(net)$size, 
				vertex.color = V(net)$color,
				  edge.color = E(net)$color, 
				edge.curved  = 0.1,
				  edge.label = NA, 
					  layout = layout_nicely
	)
legend(x=-1.5,y=-0.75,name,pch=21,col="#777777",pt.bg=unlist(cols),cex=.8,bty="n",ncol=1)
dev.off()
