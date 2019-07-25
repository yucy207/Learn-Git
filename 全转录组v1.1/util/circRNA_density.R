#!/usr/bin/env Rscript

library("docopt")

"Usage: circRNA_density.R [options] INPUT OUTPUT
Options:
  -w --width=width    the width of viewport  [default: 7]
  -h --height=height   the height of viewport [default: 7]
Arguments:
   INPUT   the input file
   OUTPUT  the output file" -> doc

opt <- docopt(doc)

width    <- as.integer(opt$w)
height   <- as.integer(opt$h)
input    <- opt$INPUT
output   <- opt$OUTPUT

library(ggplot2)

x  <- read.table(input, header=T, sep="\t", row.names=1, check.names=F, stringsAsFactors = F, quote = "", comment.char = "")

x1 <- data.frame(x[, -1])
colnames(x1) <- colnames(x)[-1]
rownames(x1) <- rownames(x)

x <- x1
x  <- log10(x + 0.00001)
#head(x)
data <- stack(x)
#head(data)
colnames(data) <- c("count", "sample")

pdf(output, width = width, height = height)
ggplot(data, aes(x =  count, fill = sample, colour = sample)) +
  scale_x_continuous(limits=c(-4,4)) + 
  geom_density( alpha = 0.65) + 
  labs(y = "Density", x = "log10(count)")
dev.off()
