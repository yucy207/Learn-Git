#!/usr/bin/Rscript

library(docopt)
"Usage:sample_rpkm.R [options] INPUT OUTPUT

Arguments:
  INPUT the input file name
  OUTPUT the output file name" -> doc


opts   <- docopt(doc)
input  <- opts$INPUT
output <- opts$OUTPUT  


x <- read.table(input, header=T,  comment.char="",  check.name=F, fill = T )
str(data)
x1    <- x[x[[26]] > 0  & x[[26]] < 0.3, ]
num1  <- nrow(x1)
rank1 <- unlist(lapply(7:26, function(t){ temp <- x1[ x1[[t]] <= 1.15 * x1[[26]]  & x1[[t]] >= 0.85 * x1[[26]] , ];nrow(temp)})) / num1

x1    <- x[x[[26]] >= 0.3  & x[[26]] < 0.6, ]
num2  <- nrow(x1)
rank2 <- unlist(lapply(7:26, function(t){ temp <- x1[ x1[[t]] <= 1.15 * x1[[26]]  & x1[[t]] >= 0.85 * x1[[26]] , ];nrow(temp)})) / num2

x1    <- x[x[[26]] >= 0.6  & x[[26]] < 3.5, ]
num3  <- nrow(x1)
rank3 <- unlist(lapply(7:26, function(t){ temp <- x1[ x1[[t]] <= 1.15 * x1[[26]]  & x1[[t]] >= 0.85 * x1[[26]] , ];nrow(temp)})) / num3


x1    <- x[x[[26]] >= 3.5  & x[[26]] < 15, ]
num4  <- nrow(x1)
rank4 <- unlist(lapply(7:26, function(t){ temp <- x1[ x1[[t]] <= 1.15 * x1[[26]]  & x1[[t]] >= 0.85 * x1[[26]] , ];nrow(temp)})) / num4

x1    <- x[x[[26]] >= 15  & x[[26]] < 50, ]
num5  <- nrow(x1)
rank5 <- unlist(lapply(7:26, function(t){ temp <- x1[ x1[[t]] <= 1.15 * x1[[26]]  & x1[[t]] >= 0.85 * x1[[26]] , ];nrow(temp)})) / num5

x1    <- x[x[[26]] >= 50, ]
num6  <- nrow(x1)
rank6 <- unlist(lapply(7:26, function(t){ temp <- x1[ x1[[t]] <= 1.15 * x1[[26]]  & x1[[t]] >= 0.85 * x1[[26]] , ];nrow(temp)})) / num6

res <- data.frame(A = rank1, B = rank2, C = rank3, D = rank4, E = rank5, F = rank6)

colnames(res) <- c(
	paste("(0-0.3) num", num1, sep=" "),
	paste("(0.3-0.6) num", num2, sep=" "),
	paste("(0.6-3.5) num", num3, sep=" "),
	paste("(3.5-15) num", num4, sep=" "),
	paste("(15-50) num", num5, sep=" "),
	paste(">=50 num", num6, sep=" ")
)

data <- stack(res)
data$pos <- rep(seq(5, 100, 5), 6)

colnames(data) <- c("val", "type", "pos")
data$type <- factor(data$type, levels = unique(data$type))
library(ggplot2)
pdf(output)
ggplot(data, aes(x = pos, y = val, colour = type, fill = type)) +
  geom_line() +
  geom_point(size=2) + 
  theme_bw() +
  theme(legend.title = element_blank(), legend.position= c(0.85, 0.15),legend.key=element_rect(colour=NA)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0,100,20), labels = seq(0, 100, 20)) +
  labs(x = "mapped reads(%)", y = "genes rpkm deviation within 15% of final value") 
dev.off()
