#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("docopt"))

"Usage: dphyer_test.R  [options]  INPUT RESULT

Options:
   -s --srgan=srgan    the organism, the value in <human, mouse, rno, Pop_tri> [default: human]

Arguments:
   INPUT   the name of input files
   RESULT  the output filename" -> doc

opts   <- docopt(doc)
input  <- opts$INPUT
output <- opts$RESULT


organ <- opts$s
x <- read.table(input, header = F, sep = "\t", stringsAsFactors = F, fill = T)


mir_num <- list(
'human'  = 2588,
'rat'    = 765,
'mouse'  = 1915,
'Pop_tri'= 405
)

pvalue  <- unlist(lapply(1:nrow(x), function(t){ 
	start <- x[t, 7]
	end   <- min(x[t, 8], x[t, 9])

	sum(unlist(lapply(start:end, function(m){

		dhyper(m, x[t, 8],  mir_num[[organ]] -  x[t, 8], x[t, 9], log = FALSE)

	})))

}))


padjust <- p.adjust(pvalue, method = "BH")
result  <- cbind(x[, 1:10], pvalue, padjust , x[,12:ncol(x)])

write.table(result, output, sep = "\t", quote = F, col.names = F, row.names = F)

