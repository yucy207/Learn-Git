rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat("============== plot tree ===========
	Usage: 
	Rscript 16S.OTU.PhyloTree.r trefile outputpath
		parameters ->
		    otufile:[file -> always otu.tax.0.03.xls];
			trefile:[file -> always TopAbundantOTUs.seq.fasta.tre];
		    outputpath:[path -> path for output]; 
	\n")		
	options("show.error.messages" = F) 
	stop()
}

otufile = normalizePath(args[1])
trefile = normalizePath(args[2])
outputpath = normalizePath(args[3])

name =  unlist(strsplit(unlist(strsplit(otufile,"/"))[length(unlist(strsplit(otufile,"/")))],"\\.xls"))

library(ggtree)
tree = read.tree(trefile)
pdf(paste(outputpath, paste(name, "_phylotree.pdf", sep = ""), sep = "/"), width = 10, height = 10)
p = ggtree(tree, layout = "circular", branch.length = "none") + geom_tiplab2(size = 3)
print(p)
dev.off()