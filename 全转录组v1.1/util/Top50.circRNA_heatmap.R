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
data1  <- data[data$type != "Not DEG", 1:index]


data1 <- data.matrix(data1)
#data = log2(data+1)

myheatcol   = colorpanel(75, 'green','black','red')
data <- data[data$type != 'Not DEG', ]
# head (data)
if( nrow(data) < 50 ){
               data <- data.matrix(data[ , 1:index])          # 小于50,全画
	}else{
               data <- data[order(data$pval), ]
               data <- data.matrix(data[1:50, 1:index])           #  相对丰度最高的50个
	}

#gene_dist   = dist(data, method='euclidean')
#hc_genes    = hclust(gene_dist, method='complete')
#data        = t(scale(t(data), scale=F))

#sample_dist = dist(t(data), method='euclidean')
#hc_samples  = hclust(sample_dist, method='complete')

#heatmap_file <- paste(output, "/Top50.heatmap.pdf", sep="")
pdf(output , width=as.numeric(opts$w), height = as.numeric(opts$h) )
heatmap.2(data, 
		dendrogram   = "both", 
		#Rowv         = as.dendrogram(hc_genes), 
		#Colv         = as.dendrogram(hc_samples), 
		# Colv         = FALSE,
	      col          = myheatcol, 
scale="row", 
      margins      = c(5,11),
      density.info = "none",
      trace="none",
      key=TRUE,
      keysize=0.6,
      # cexRow=0.01,
      cexRow=0.8,
      cexCol=0.8,
     #srtCol=180
      )

dev.off()



