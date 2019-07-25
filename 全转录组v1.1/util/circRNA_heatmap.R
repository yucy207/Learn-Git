#!/usr/bin/env Rscript

library("docopt")

"Usage: circRNA_heatmap.R  [options]  INPUT RESULT
Options:
   -w --width=width    the width of viewport  [default: 10]
   -h --height=width   the height of viewport [default: 10]
   -sk --skip=0   skip number                 [default: 0]
Arguments:
   INPUT   the name of input files
   RESULT  the output filename" -> doc


opts   <- docopt(doc)

input  <- opts$INPUT
output <- opts$RESULT

library("gplots")

data  <-  read.table(input, header=T, com='', sep="\t", row.names=1, check.names=F, skip = opts$sk)
index <- ncol(data)-9
data  <- data[data$type != "Not DEG", 1:index]
#head(data)

data <- data.matrix(data)
#data = log2(data+1)

myheatcol   = colorpanel(75, 'green','black','red')
#myheatcol   = colorpanel(75, 'purple','black','yellow')
#gene_dist   = dist(data, method='euclidean')
#hc_genes    = hclust(gene_dist, method='complete')
#data        = t(scale(t(data), scale=F))

#sample_dist = dist(t(data), method='euclidean')
#hc_samples  = hclust(sample_dist, method='complete')
#gene_dist   = dist(data, method='euclidean')
#hc_genes    = hclust(gene_dist, method='complete')


pdf(output, width=as.numeric(opts$w), height = as.numeric(opts$h) )
#png(output, width=600, height=800)
heatmap.2(data, 
		dendrogram   = "both", 
		#Rowv         = as.dendrogram(hc_genes), 
		#Colv         = as.dendrogram(hc_samples), 
		#Colv         = FALSE,
	      col          = myheatcol, 
	scale="row", 
      margins      = c(2,2),
      density.info = "none",
      trace="none",
      key=TRUE,
      keysize=0.6,
      cexRow=0.01,
      cexCol=0.8,
     #srtCol=180
      )

dev.off()




