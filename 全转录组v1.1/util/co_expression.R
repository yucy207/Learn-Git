#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(docopt))

"Usage: co_expression.R [options] MRNA LNCRNA ID1 ID2 OUTPUT 

Arguments:
  MRNA   the mRNA count table
  LNCRNA the LNCRNA count table
  ID1    the first type name
  ID2    the second type name
  OUTPUT the output file name" -> doc

opts   <- docopt(doc)
mRNA   <- opts$MRNA
lnRNA  <- opts$LNCRNA
output <- opts$OUTPUT
first  <- opts$ID1
second <- opts$ID2

oldw <- getOption("warn")
#options(warn = -1)




x <- read.table(mRNA,  header=F, sep="\t", row.names=1)
y <- read.table(lnRNA, header=F, sep="\t", row.names=1)


#if (ncol(x) < 3 || ncol(y) < 3) stop("the sample must bigger than 3!\n")



cal_cor <- function(y, a, mRNA) {

	result <- t(sapply(1:nrow(y), function(t) { 
		lncRNA <- rownames(y[t,])
		temp   <- unlist(y[t,])
		#res    <- cor(temp, a)
		suppressWarnings(res    <- cor.test(temp, a, method = "spearm"))
		c(mRNA, lncRNA, res$estimate, res$p.value)
		#c(mRNA, lncRNA, res)
	}))
	result
}

result <- lapply(1:nrow(x), function(t){

	mRNA <- rownames(x[t,])
	temp <- unlist(x[t,])
	cal_cor(y, temp, mRNA)

})

result <- data.frame(do.call(rbind, result))

#colnames(result) <- c("mRNA", "lncRNA", "cor")
colnames(result) <- c(first, second, "cor", "pvalue")



write.table(result, output, row.names=F, sep="\t", quote=F)

options(warn = oldw)
