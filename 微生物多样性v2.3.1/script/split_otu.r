rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
	"============== OTU_table_split ===========
	Usage: 
	Rscript OTU_table_split otufile groupfile outputpath
		parameters ->
			  otufile: [file -> always subsample_otu_table.txt];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}
otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])


SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]


OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")

otu_table = OTUFile[,rownames(SampleInfo)]
datasum = rowSums(otu_table)
OTUsize = as.matrix(datasum)
otu_table = cbind(otu_table,OTUsize)

otu_table = cbind(rownames(OTUFile),otu_table)
loc = grep("Taxonomy", colnames(OTUFile))
tax = OTUFile[(loc):(loc+7)]
otu = cbind(otu_table,tax)

otu = otu[otu$OTUsize != 0,]
map.ID = cbind(rownames(otu),rownames(otu))

colnames(otu)[1] = "OTU Id"

name = strsplit(otufile,"/")[[1]][length(strsplit(otufile,"/")[[1]])]
write.table(otu, paste(outputpath, paste(name) ,sep = "/"), row.names = F,sep = "\t", quote = F)
write.table(map.ID, paste(outputpath, "map.ID",sep = "/"), row.names = F,col.names = F, sep = "\t", quote = F)

