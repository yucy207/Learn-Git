rm(list = ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
    "============== do betaNTI ===========
    Usage: 
        Rscript 16S.BetaDiversity.select_stat.RA.r stat.RA groupfile outputpath
			parameters ->
			  otufile: [file -> always otu.tax.0.03.stat.RA.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}


otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])
setwd(outputpath)

# ====================== prepare data =======================
# SampleInfo
SampleInfo = read.table(groupfile, header = F, sep = "\t")
rownames(SampleInfo) = as.character(SampleInfo[,1])
colnames(SampleInfo) = c("Sample", "Group")


OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
otu_data = OTUFile[,rownames(SampleInfo)]
otu_data = otu_data[rowSums(otu_data) != 0,]
otu_data = cbind(rownames(otu_data), otu_data)
colnames(otu_data)[1] = "Taxlevel"
write.table(otu_data,"stat.RA.xls", row.names = F, sep = "\t", quote = F)