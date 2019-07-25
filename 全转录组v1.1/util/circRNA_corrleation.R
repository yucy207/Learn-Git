#!/usr/bin/env Rscript

library(docopt)
"Usage: correlation.R [options] INPUT OUTPUT 

Arguments:
  INPUT   the input file name 
  OUTPUT  the output file name"->doc


opts   <- docopt(doc)
input  <- opts$INPUT
output <- opts$OUTPUT



x <- read.table(input, header=T, sep="\t", row.names=1, stringsAsFactors=F)
index_1 <- ncol(x) - 5
index_2 <- ncol(x) - 6

result <- data.frame(x = log10(x[[index_1]]), y=log10(x[[index_2]]))

pdf(output)
library(ggplot2)
ggplot(result, aes(x=x,y=y)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color="red") + 
  scale_x_continuous(limits=c(-2.5, 5)) + 
  scale_y_continuous(limits=c(-2.5, 5)) + 
  labs(x="log10(count)", y="log10(count)")
dev.off()

