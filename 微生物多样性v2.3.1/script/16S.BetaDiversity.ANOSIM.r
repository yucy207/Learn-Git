rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat("Usage: 
	Rscript 16S.BetaDiversity.ANOSIM.r otufile groupfile outputpath
		parameters ->
			    otufile: [file -> always otu.tax.0.03.xls];
			  groupfile: [file -> always sample.groups without header and only two column];
		         outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

library(MASS)
library(vegan)

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

setwd(outputpath)

# ====================== prepare data =======================
SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]

OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
otu_data = OTUFile[,rownames(SampleInfo)]

loc = grep("All", rownames(otu_data)) 
if (length(loc) > 0) otu_data = otu_data[-loc, ]

DistBC = vegdist(t(otu_data), method = "bray")

# ANOSIM
an = anosim(DistBC, SampleInfo$Group, permutations = 9999)
permu = paste("Number of permutations", an$permutations, sep = ": ")
dissi = paste("Dissimilarity", an$dissimilarity, sep = ": ")
stati = paste("statistic R", an$statistic, sep = ": ")
signi = paste("Significance", an$signif, sep = ": ")
result = rbind(permu, dissi, stati, signi)
colnames(result) = "ANOSIM"
write.csv(result, "ANOSIM.csv", row.names = FALSE, quote = FALSE)
