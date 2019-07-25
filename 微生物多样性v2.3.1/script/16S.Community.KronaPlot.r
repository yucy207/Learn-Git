# https://github.com/marbl/Krona/wiki/KronaTools

rm(list = ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
	"============== KronaPlot using species.taxon.Abundance.xls produced by 16S.CommunityStructure.OTUSplit.r ===========
	Usage: 
	Rscript 16S.Community.KronaPlot.r taxonfile outputpath
		parameters ->
			taxonfile: [file -> always species.taxon.Abundance.xls produced by 16S.CommunityStructure.OTUSplit.r];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

taxonfile = normalizePath(args[1])
outputpath = normalizePath(args[2])

# ====================== prepare data =======================
taxonfile = read.table(taxonfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)

loc = grep("All|No_Rank", rownames(taxonfile))
if(length(loc) > 0) taxon_data = taxonfile[-loc,]

loc2 = which(rownames(taxonfile) == "All")

if(length(loc2) > 0) taxon_data = taxon_data[-loc2,]

index = grep("Abundance|OTUsize", colnames(taxon_data))

sapply(1:(index-1), function(i) {

	A = cbind(taxon_data[,i], taxon_data[,(index+1):ncol(taxon_data)], rownames(taxon_data))
	inputfile = paste(outputpath, paste(colnames(taxon_data)[i], "_data_for_Krona.xls", sep = ""), sep = "/")
	write.table(A, inputfile, sep = "\t", quote = F, row.names = F, col.names = F)
	
	outputfile = paste(outputpath, paste(colnames(taxon_data)[i], ".krona.html", sep = ""), sep = "/")
	
	cmd = paste("/home/panrf/Softwares/KronaTools-2.7/bin/ktImportText ", inputfile, " -o ", outputfile, sep = "")
	system(cmd)
	
})

