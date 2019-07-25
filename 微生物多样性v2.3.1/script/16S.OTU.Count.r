#!/usr/bin/Rscript
rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
	"============== OTU number count ===========
	Usage: 
	Rscript 16S.OTU.Count.r otufile outputpath
		parameters ->
			otufile: [file -> always subsample_otu.tax.0.03.xls];
		     outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

otufile = normalizePath(args[1])
outputpath = normalizePath(args[2])

data = read.table(otufile, header = T, sep = "\t", row.names = 1,check.names = F, comment.char = "", quote = "", fill = T)
loc = grep("OTUsize",colnames(data))

taxon = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

for (t in 1:length(taxon)){

	data1 = data[, c(as.character(colnames(data)[1:loc-1]), taxon[t])]
	data1 = as.matrix(data1[as.character(data1[,ncol(data1)]) != "Unassigned",])
    data1 = as.matrix(data1[as.character(data1[,ncol(data1)]) != "", -ncol(data1) ])
	
	if(t == 1){
		 
		count = apply(data1, 2, function(d) length(which(as.numeric(d) != 0)))
		
	}else{
		
		count = rbind(count, apply(data1, 2, function(d) length(which(as.numeric(d) != 0))))
	
	}
}

rownames(count) = taxon
colnames(count) = colnames(data)[1:loc-1]
count = cbind(rownames(count),count)
colnames(count)[1] = "Taxon"

filename = unlist(strsplit(otufile, "/|\\.xls"))[length(unlist(strsplit(otufile, "/|\\.xls")))]
write.table(count, paste(outputpath, paste(filename, "_count.xls", sep = ""), sep = "/"), sep = "\t", row.names = F, quote = F)

