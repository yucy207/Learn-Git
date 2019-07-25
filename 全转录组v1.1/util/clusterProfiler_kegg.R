#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=T)


if(length(args) != 3){

        cat(
        "============== do kegg enrich ===========
        Usage: 
        Rscript clusterProfiler_kegg.R de_gene.list db out_dir
                parameters ->
                    de_gene.list: [file -> always de_gene.list];
                    db : [always pathway.db]
                    out_dir: [path -> path for output]; \n")
        options("show.error.messages" = F) 
        stop()
}


gene    <- args[1]
db      <- args[2]
out_dir <- args[3]


dir.create(out_dir, showWarnings = FALSE)

library(clusterProfiler)


gda <- read.table(db,  header = F, sep = "\t" ,stringsAsFactors = F, quote = "")

disease2gene = gda[, c(3,1)]

disease2name = gda[, c(3,4)]

deg = unlist(read.table(gene,  header = F, sep = "\t" ,stringsAsFactors = F))

x = enricher(deg, 
	pvalueCutoff = 1,
	qvalueCutoff = 1,
	minGSSize = 1,
	TERM2GENE = disease2gene, 
	TERM2NAME = disease2name)

res <- data.frame(summary(x))


output <- paste(out_dir, "kegg_enrichment.xls", sep="/")
write.table(res, output, sep="\t", row.names = F, quote = F)






