#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=T)


args = commandArgs(trailingOnly=T)


if(length(args) != 3){

        cat(
        "============== do go enrich ===========
        Usage: 
        Rscript clusterProfiler_go.R de_gene.list db out_dir
                parameters ->
                    de_gene.list: [file -> always de_gene.list];
                    db : [always go.annotation.db]
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

gda_bp <- gda[gda[[4]] == "BP", ]

gda_mf <- gda[gda[[4]] == "MF", ]

gda_cc <- gda[gda[[4]] == "CC", ]


disease2gene_bp <- gda_bp[, c(2,1)]

disease2name_bp <- gda_bp[, c(2,3)]


disease2gene_cc <- gda_cc[, c(2,1)]

disease2name_cc <- gda_cc[, c(2,3)]


disease2gene_mf <- gda_mf[, c(2,1)]

disease2name_mf <- gda_mf[, c(2,3)]


deg <- unlist(read.table(gene,  header = F, sep = "\t" ,stringsAsFactors = F))


deg_bp <- deg[deg %in% gda_bp[[1]]]
deg_mf <- deg[deg %in% gda_mf[[1]]]
deg_cc <- deg[deg %in% gda_cc[[1]]]


result <- list()

if (length(deg_bp) > 0) {
	x_bp <- enricher(deg_bp, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1,TERM2GENE=disease2gene_bp, TERM2NAME=disease2name_bp)
	res_bp <- as.data.frame(x_bp)
	res_bp$type <- "BP"

	result[["BP"]] <- res_bp
}

if (length(deg_cc) > 0) {
	x_cc <- enricher(deg_cc, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, TERM2GENE=disease2gene_cc, TERM2NAME=disease2name_cc)
	res_cc <- as.data.frame(x_cc)
	res_cc$type <- "CC"	

	result[["CC"]] <- res_cc
}	

if (length(deg_mf) > 0) {
	x_mf <- enricher(deg_mf, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1, TERM2GENE=disease2gene_mf, TERM2NAME=disease2name_mf)
	res_mf <- as.data.frame(x_mf)

	res_mf$type <- "MF"

	result[["MF"]] <- res_mf
}


res <- do.call(rbind, result)

output <- paste(out_dir, "go_enrichment.xls", sep="/")
write.table(res, output, sep="\t", row.names = F, quote = F)




