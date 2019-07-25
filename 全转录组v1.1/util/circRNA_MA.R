#!/usr/bin/env Rscript

library(docopt)
"Usage: circRNA_MA.R [options] INPUT OUTPUT

Arguments:
  INPUT  the input file name 
  OUTPUT the output file name" -> doc

opts   <- docopt(doc)
input  <- opts$INPUT
output <- opts$OUTPUT 

x <- read.table(input, header=T, sep="\t", row.names=1, stringsAsFactors=F)

a1_inedx <- ncol(x) - 5
a2_index <- ncol(x) - 6

a3_index <- ncol(x) - 3
result <- data.frame(M = x[[a3_index]], A = (x[[a1_inedx]] + x[[a2_index]])/2, significant=x[[ncol(x)]])
result$significant <- factor(result$significant, levels = c("Up", "Down", "Not DEG"))
pdf(output)
library(ggplot2)
ggplot(result, aes(x = A, y = M, colour = significant)) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 15)) +   
  scale_y_continuous(limits = c(-20, 20)) + 
  theme(legend.title = element_blank()) + 
  labs(x = "A", y="log2(FC)", tilte="MA plot")
dev.off()

