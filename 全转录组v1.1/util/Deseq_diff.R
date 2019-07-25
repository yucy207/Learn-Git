
#!/usr/bin/env Rscript

library(docopt)
"Usage: Deseq_diff.R [options] LOG2FC CASE CONTROL CASE_GROUP CONTROL_GROUP INPUT OUTPUT

Arguments:
  LOG2FC           the log2 fold change
  CASE             the case group sample
  CONTROL          the control group sample
  CASE_GROUP       the case group name
  CONTROL_GROUP    the control group name
  INPUT            the input count table file
  OUTPUT           the output file name" -> doc 

opts          <- docopt(doc)
log2fc        <- as.numeric(opts$LOG2FC)
case          <- opts$CASE
control       <- opts$CONTROL
case_group    <- opts$CASE_GROUP
control_group <- opts$CONTROL_GROUP
input         <- opts$INPUT
output        <- opts$OUTPUT



countTable <- read.table(input, header=T, sep="\t", check.names=F, row.names=1)
gene       <- countTable[[1]]
countTable <- countTable[, -1]

library(DESeq)
control_count           <- sapply(unlist(strsplit(control, split=",")), function(t){countTable[[t]] })
colnames(control_count) <- unlist(strsplit(control, split=","))

case_count              <- sapply(unlist(strsplit(case, split=",")), function(t){ countTable[[t]] })
colnames(case_count)    <- unlist(strsplit(case, split=","))

count              <- cbind(data.frame(control_count, check.names=F), data.frame(case_count, check.names=F))
rownames(count)    <- rownames(countTable)
countTable         <- count[apply(count, 1, sum) > 0 , ]

condition          <- factor( c(rep(control_group, ncol(control_count)), rep(case_group, ncol(case_count)) ), levels =c(control_group, case_group))



str(countTable)

condition

cds <- newCountDataSet( countTable, condition )


cds <-  estimateSizeFactors( cds )
sizeFactors(cds)

norm <- as.data.frame(counts( cds, normalized=TRUE ))
norm$gene <- gene[apply(count, 1, sum) > 0]

method <- "pooled"

if (ncol(case_count) > 1 && ncol(control_count) > 1 ) {
 method <- "pooled"
 sharingMode <- "maximum"
} else {
 method <- "blind"
 sharingMode <- "fit-only"
}
fitType <- "local"

cds <- estimateDispersions( cds, method = method, fitType = fitType,sharingMode=sharingMode )


res <- nbinomTest( cds, control_group, case_group )

# colData            <- data.frame(row.names = colnames(countTable), condition = rep(c(case_group, control_group), times = c(ncol(case_count), ncol(control_count)) ))

# dds <- DESeqDataSetFromMatrix(countData = countTable, colData = colData, design = ~ condition)


res$type <- "Not DEG"
res$type[res$log2FoldChange > log2fc & res$pval < 0.05] <- "Up"
res$type[res$log2FoldChange < -log2fc & res$pval < 0.05] <- "Down"



res  <- cbind(norm,res[,-1])
write.table(res, output, sep="\t", quote=F, col.names=NA)






