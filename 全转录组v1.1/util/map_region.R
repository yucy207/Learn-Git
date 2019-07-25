#!/usr/bin/Rscript


library(docopt)
library(ggplot2)

"Usage: map_region.R [optinons] INPUT OUTDIR SAMPLE

Argument:
  INPUT     the input file name
  OUTDIR    the out directory
  SAMPLE    the sample name
" -> doc

opts    <- docopt(doc)
input   <- opts$INPUT
out_dir <- opts$OUTDIR
sample  <- opts$SAMPLE




x  <- read.table(input, header=T, sep="\t", row.names=1 ,stringsAsFactors=F, check.names=F)

#rownames(x)
x <- t(data.matrix(x[rownames(x) == sample, 2:4 ]))

#x
pdf_name <- paste(out_dir, "/", sample, ".map_region.pdf",sep="")

x <- data.frame(x)
data <- data.frame(region = rownames(x), freq  = x[[1]])

data$freq <- as.numeric(as.vector(data$freq))
data$region <- factor(data$region, level = unique(data$region))
data <- data[order(data$freq), ]



label <- sprintf("%.2f%%", data$freq * 100)
data$ymax <- cumsum(data$freq)
data$ymin <- c(0, cumsum(data$freq)[-nrow(data)])
freq <- data$freq
position <- unlist(lapply(1:length(freq), function(t){ cumsum(freq)[t]-(freq[t]/2) }))

text     <- data.frame(location=position, text=label)

pdf(pdf_name)
ggplot(data=data) + 
  		geom_rect(aes(xmin=1, xmax=3, ymin=ymin, ymax=ymax, fill=region)) + 
  		geom_text(data= text,aes(x=2.8,y=location,label=text),size=3) + 
  		coord_polar(theta="y") +
  		theme_bw() + 
  		theme(
		axis.ticks.x     = element_line(size=2),
		axis.ticks.y      = element_blank(),
        	axis.text         = element_blank(),
        	panel.border      = element_rect(colour=NA),  
        	panel.background  = element_rect(colour=NA),
        	axis.title        = element_blank(),
        	panel.grid        = element_blank(),
		plot.title        = element_text(hjust = 0.5),
        	legend.key        = element_rect(colour=NA))

dev.off()

