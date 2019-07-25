#!/usr/bin/Rscript

library(docopt)

"Usage:chr_depth.R [options] INPUT OUTPUT

Arguments:
  INPUT  the input file name
  OUTPUT the output file name" -> doc

opts   <-  docopt(doc)
input  <-  opts$INPUT
output <-  opts$OUTPUT


library(ggplot2)
library(gtools)
library(grid)

x    <- read.table(input, sep="\t", header=F)
#if (substr(x$V3[1],1,3) != "chr") {
#	x$V3 <- paste("chr", x$V3,sep="")
#}
#x$V3 <- factor(x$V3, levels = mixedsort(unique(x$V3)))
x$V2 <- log2(x$V2 + 1)

x1 <- x[x$V4 == "plus",]
x2 <- x[x$V4 == "minus",]
x2$V2 <- -x2$V2
 
x2$V3 <- factor(x2$V3, levels = mixedsort(unique(x2$V3)))
x1$V3 <- factor(x1$V3, levels = mixedsort(unique(x1$V3)))
max_val <- max(x$V1)
str(x2)
head(x2)


p <- ggplot(x2, aes(x = V1, y = V2, fill = V4, colour = V4)) + 
  geom_bar(stat="identity") + 
  geom_bar(data = x1, aes(x = V1, y= V2, fill = V4),stat="identity") + 
  facet_grid(V3~.) +
  labs(x = "chromosome position (Mb)", y = "Median of read density(log2)", title = "Reads Density in chromosome") +
  theme(
  		#axis.line.x     = element_line(size = 10, colour = "black"),
  		#axis.line.y     = element_line(size = 10, colour = "black"),
  		strip.text.y = element_text(size = 10, angle = 0),
  		strip.background = element_blank(),
  		plot.background  = element_rect(fill=NA),
  		plot.title = element_text(hjust = 0.5),
		panel.background = element_rect(fill=NA),
  		panel.grid      = element_blank(),
  		panel.border    = element_rect(fill=NA),
  		legend.position = "nonoe",
  		axis.ticks.x    = element_line(size = 1),
  		axis.ticks.y    = element_line(size = 1),
  		axis.ticks.length = unit(2, "mm"),
  		axis.text.x     = element_text(size = 8),
  		axis.text.y    = element_text(size = 8)) +
  scale_x_continuous(limits = c(0, max_val), breaks = seq(0, max_val,50000), labels = paste(seq(0, max_val, 50000)/ 10000, "M", sep="")) + 
  scale_y_continuous(limits = c(-10, 10), breaks= seq(-10, 10, 10), labels =c(-10, 0, 10))


pdf(paste(output, ".pdf", sep = ""))
p
dev.off()

png(paste(output, ".png", sep = ""))
p
dev.off()
