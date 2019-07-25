#!/usr/bin/env Rscript

library(docopt)
"Usage: ballgown_fpkm.R [options] INPUT  OUTPUT

Arguments:
  INPUT          the input file name
  OUTPUT         the output directory name" -> doc


opts          <- docopt(doc)
input         <- opts$INPUT
output        <- opts$OUTPUT


library(ballgown)


samples = unlist(strsplit(input, split = ","))

bg = ballgown(samples = samples, meas='all')

trans <- texpr(bg, 'all')

str(trans)
trans <- trans[, c(6, 9,seq(12, ncol(trans), 2))]

genes <- gexpr(bg)

write.table(genes, paste(output, "genes.fpkm.xls", sep = "/"), quote = F, sep = "\t", col.names = NA)
write.table(trans, paste(output, "trans.fpkm.xls", sep = "/"), quote = F, sep = "\t", row.names = F)

