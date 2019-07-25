#!/usr/bin/env Rscript

library("docopt")

"Usage: maSigPro.R  [options]  INPUT GROUP RESULT TYPE
Options:
   -w --width=width    the width of viewport  [default: 12]
   -h --height=width   the height of viewport [default: 12]
   -m --minobs=3       the num of min obs [default: 3]
   -sk --skip=0   skip number                 [default: 0]
Arguments:
   INPUT   the name of input files
   GROUP   the design matrix
   RESULT  the out direactory filename
   TYPE    gene or transcript" -> doc


opts   <- docopt(doc)

input  <- opts$INPUT
obs    <- as.integer(opts$m)
design <- opts$GROUP
output <- opts$RESULT
type   <- opts$TYPE

library(maSigPro)

group  <- read.table(design, header = T, sep = "\t", row.names = 1, check.names = F, comment = "")
data   <- read.table(input,  header = T, sep = "\t", row.names = 1, check.names = F, comment = "", stringsAsFactors = F)
if(type == "transcript"){
  data1<- data
  data <- data1[,-1]
}

points <- length(unique(group[[1]]))
design <- make.design.matrix(group, degree = points - 1 )

fit    <- p.vector(data, design, Q = 1, MT.adjust = "BH", min.obs = obs)
tstep  <- T.fit(fit, step.method = "backward", alfa = 1)
sigs   <- get.siggenes(tstep, rsq = 0.6, vars = "each")

#ALL
res <- cbind(sigs$sig.genes$independ$sig.profiles,
	         sigs$sig.genes$independ$coefficients,
             #sigs$sig.genes$independ$group.coeffs,
             sigs$sig.genes$independ$sig.pvalues)

for(i in 1:nrow(res)){

	if(res$"p-value"[i] < 0.01){
		res$type[i] <- "DEG";
	}else{
		res$type[i] <- "Not DEG";
	}
}

if(type == "transcript"){
  geneid <- data1[which(rownames(data1) %in% rownames(res)),]$gene_id
  res    <- cbind(geneid,res)
}

res  = res[order(res$"p-value"), ]
file <- paste(output, "/all.diff.xls", sep = "")
write.table(res, file, quote = F, sep = "\t", col.names = NA)

#every time
for (i in 1:(points - 1)) {
	name <- paste("Time", i, sep = "")
	if (i == 1) name <- "Time"
	res <- cbind(sigs$sig.genes[[name]]$sig.profiles,
	         sigs$sig.genes[[name]]$coefficients,
             #sigs$sig.genes[[name]]$group.coeffs,
             sigs$sig.genes[[name]]$sig.pvalues)

	for(j in 1:nrow(res)){

		if(res$"p-value"[j] < 0.01){
			res$type[j] <- "DEG";
		}else{
			res$type[j] <- "Not DEG";
		}
	}

  if(type == "transcript"){
    geneid <- data1[which(rownames(data1) %in% rownames(res)),]$gene_id
    res    <- cbind(geneid,res)
  }

  res  = res[order(res$"p-value"), ]
	out <- paste(output, "/", name, ".diff.xls", sep = "")
	write.table(res, out, quote = F, sep = "\t", col.names = NA)

}



