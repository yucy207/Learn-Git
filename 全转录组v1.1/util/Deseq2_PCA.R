#!/usr/bin/env Rscript

library(docopt)
"Usage: Deseq2_PCA.R [options] INPUT  OUTPUT


Arguments:
  INPUT          the input file name
  OUTPUT         the output directory name" -> doc


opts          <- docopt(doc)
input         <- opts$INPUT
output        <- opts$OUTPUT



library(DESeq2)

countData <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)

colData   <- data.frame(row.names = colnames(countData), 
	                    condition = rep("case", times = ncol(countData))        
	                    )

dds <- DESeqDataSetFromMatrix(countData =  countData, colData = colData, design = ~ 1)

ntd <- assay(normTransform(dds))

vsd <- assay(varianceStabilizingTransformation(dds, blind=FALSE))

rld <- assay(rlog(dds, blind=FALSE))


write.table(ntd, paste(output, "count.norm.xls", sep = "/"), quote = F, sep = "\t", col.names = NA)
write.table(vsd, paste(output, "count.vst.xls", sep = "/"), quote = F, sep = "\t", col.names = NA)
write.table(rld, paste(output, "count.rlog.xls", sep = "/"), quote = F, sep = "\t", col.names = NA)


