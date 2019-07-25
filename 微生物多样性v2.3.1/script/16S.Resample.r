rm(list = ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	"============== do Resample and generate subsample_otu.wang.taxonomy and subsample_otu_table.txt ===========
	Usage: 
	Rscript 16S.Resample.r otufile groupfile trimOTUs outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
			groupfile: [file -> always sample.groups without header and only two column];
			 trimOTUs: [character -> if trim OTUs, option: T or F];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(phyloseq) 
library(gtools)
library(ape)

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
trimOTUs = toupper(as.character(args[3]))
outputpath = normalizePath(args[4])

# ====================== prepare data =======================
SampleInfo = read.table(groupfile, header = F, sep="\t")
rownames(SampleInfo) = as.character(SampleInfo[,1])

OTUData = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
otu_data = OTUData[,rownames(SampleInfo)]


OTU = otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)
SAMPLE = sample_data(SampleInfo)
sample_names(SAMPLE) = as.character(rownames(SampleInfo))

physeq = phyloseq(OTU, SAMPLE)

# prune OTUs that are not present in at least one sample
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

# Add phy_tree data
tree = rtree(ntaxa(physeq), rooted = TRUE, tip.label = taxa_names(physeq))
physeq = merge_phyloseq(physeq, tree)

# subsample
sub = rarefy_even_depth(physeq, replace = F, trimOTUs = trimOTUs, verbose = F)
sub_otu = otu_table(sub)
sub_otu = sub_otu[mixedsort(rownames(sub_otu)),]      # gtools package

loc = grep("Taxonomy", colnames(OTUData))
if (length(loc) > 0) {

	taxonomy = cbind(rownames(sub_otu), as.character(OTUData[rownames(sub_otu),loc]))
	write.table(taxonomy, paste(outputpath, "subsample_otu.wang.taxonomy", sep = "/"), row.names = F, col.names = F, sep = "\t", quote = F)

}

sub_otu = cbind(rownames(sub_otu), sub_otu)
colnames(sub_otu)[1] = "OTUId"

write.table(sub_otu, paste(outputpath, "subsample_otu_table.txt", sep = "/"), row.names = F, sep = "\t", quote = F)

