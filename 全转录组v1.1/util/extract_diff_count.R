#!/usr/bin/Rscript

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(docopt))

"Usage:extract_diff_count.R [options] CONTROL_SAMPLE CASE_SAMPLE INPUT DE  OUTPUT

Arguments:
  CONTROL_SAMPLE  the control sample name,eg.  
  CASE_SAMPLE     the case sample naem,eg.
  INPUT           the input file name
  DE              the sig_deseq_known.xls file
  OUTPUT          the output file name"-> doc

opts    <- docopt(doc)
control <- opts$CONTROL_SAMPLE
case    <- opts$CASE_SAMPLE
input   <- opts$INPUT
df      <- opts$DE
output  <- opts$OUTPUT

x <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)

y <- read.table(df, header=T, sep="\t", row.names=1, comment.char="", check.names=F)
y <- y[y$type %in% c("Up", "Down"),]
control_count           <- sapply(unlist(strsplit(control, split=",")), function(t){x[[t]] })

colnames(control_count) <- unlist(strsplit(control, split=","))


case_count             <- sapply(unlist(strsplit(case, split=",")), function(t){ x[[t]] })
colnames(case_count)   <- unlist(strsplit(case, split=","))

count                  <- cbind(data.frame(control_count, check.names=F), data.frame(case_count, check.names=F))

rownames(count)        <- rownames(x)
#count <- count[apply(count, 1, sum) > 0,]
count <- count[rownames(count) %in% rownames(y), ] 
write.table(count, output, col.names=F, quote=F, sep="\t")
