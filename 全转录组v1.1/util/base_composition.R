#!/usr/bin/env Rscript

"Usage: base_composition.R INPUT OUTPUT NAME
Arguments:
   INPUT    the input file name  
   OUTPUT   the output filename
   NAME     the title name" -> doc

library("docopt")
opts     <- docopt(doc)
input    <- opts$INPUT
output   <- opts$OUTPUT
name     <- opts$NAME   

library("ggplot2")

x <- read.table(input, header=F, sep="\t")


data <- x[,3:7]  / x[[2]] * 100
colnames(data) <- c("A", "C", "G", "T", "N")


x2       <- stack(data)
head(x2)
x2$group <- rep(1:nrow(data), ncol(data))

x2$ind   <- factor(x2$ind,   level=unique(x2$ind))
colnames(x2) <- c("values", "type", "pos")

height_f <- 11
if ( nrow(x) < 100 )
   height_f <- 7



g <- ggplot(x2,aes(x=pos, y=values,group=type)) + 
  geom_line(aes(colour=type)) +
  labs( x    = "Position along reads",
            y    = "Percent of bases",
       title = paste("Base content along read(", name, ")" , sep=""))  +
  geom_vline(xintercept=nrow(x) / 2, linetype=2, colour="skyblue") +   
  scale_x_continuous(expand=c(0,0),limits=c(0, nrow(x)), breaks=seq(0, nrow(x), 50), labels=seq(0, nrow(x), 50)) +
  theme(plot.title = element_text(hjust = 0.5))
  
pdf(output, height = height_f, width=10)
g
dev.off()
