#!/usr/bin/env Rscript

library(docopt)
"Usage: circRNA_valcano.R [options] INPUT OUTPUT 

Arguments:
  INPUT  the input file name 
  OUTPUT the output file name"->doc


opts   <- docopt(doc)
input  <- opts$INPUT
output <- opts$OUTPUT 




x <- read.table(input, header=T, sep="\t", row.names=1, stringsAsFactors=F)

type_index <- ncol(x)
x_index    <- ncol(x) - 3
y_index    <- ncol(x) - 2

data <- x[,c(x_index,y_index,type_index)]

data[[2]] <- -log10(data[[2]])
colnames(data) <- c("x", "y", "type")

data$type <- factor(data$type, levels = c("Up", "Down","Not DEG" ))

library(ggplot2)


pdf(output)
ggplot(data, aes(x = data[[1]], y = data[[2]], color = type)) + 
  geom_point() +
  theme_bw() + 
  #scale_colour_manual(values = c("green","red")) + 
  labs(x = bquote(paste(log[2],"(fold change)",sep="")), y = bquote(paste(-log[10],"(p value)",sep="")), title = "genes: CASE/CONTROL") +
  scale_x_continuous(limits=c(-10,10)) +
  theme(plot.title = element_text(hjust = 0.5))  
 #geom_hline(yintercept = 1.3, color="black", linetype=3)

dev.off()
