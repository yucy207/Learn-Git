

#!/usr/bin/env Rscript

library(docopt)

"Usage:error_rate.R [options] INPUT OUTPUT SAMPLE

Options:
  -q--quality  the fastq quallity score[default: 33]

Arguments:
  INPUT   the input file name
  OUTPUT  the output file name
  SAMPLE  the sample name" -> doc

opts    <- docopt(doc)
input   <- opts$INPUT 
quality <- as.integer(opts$q)
output  <- opts$OUTPUT
sample  <- opts$SAMPLE


x <- read.table(input, header = F)[,c(1,8)]
colnames(x) <- c("position", "error")
x$position <- 1:nrow(x)
x$error <- 10 ^ ((x$error - quality) / -10) * 100
center <- nrow(x) / 2

sample_title <- paste("Error rate distribution along reads(", sample, ")", sep="")
library(ggplot2)


color <- "#32CD32"
fill_color <- "#00FF7F"




g <-ggplot(x, aes(x = position, y = error)) + 
  geom_bar(fill=fill_color,color=color, stat="identity") +
  geom_vline(xintercept = center, colour="skyblue", linetype=2) +
  labs(x = "Position along reads", y = "%Error rate", title = sample_title)  + 
  theme(plot.title = element_text(hjust = 0.5))


pdf(output)
g
dev.off()

