#!/usr/bin/Rscript

library(docopt)
"Usage:ggplot_hist.R [options] INPUT OUTPUT

Options:
  -b--bins=50  hist bin number[default: 50]

Arguments:
  INPUT  the input file name
  OUTPUT the output pdf file name"->doc

opts   <- docopt(doc)
input  <- opts$INPUT
output <- opts$OUTPUT
width  <- as.integer(opts$b)

x <- read.table(input, header=F)
library(ggplot2) 

pdf(output)
ggplot(x, aes(V1)) + 
  geom_histogram(bins = width, fill = "brown") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "circRNA length", y = "count", title = "circRNA length distribution")
dev.off()
