#!/usr/bin/env Rscript

library(docopt)
"Usage:pca_plot.R [options] INPUT GROUP OUTPUT

Arguments:
  INPUT    the input file name
  GROUP    the group file name
  OUTPUT   the output file name" -> doc

opts   <- docopt(doc)
input  <- opts$INPUT
groupfile  <- opts$GROUP
output <- opts$OUTPUT

library(ggrepel)
library(Rtsne)

SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Group = as.character(SampleInfo[,2])

inputdata <- read.table(input, header=T, sep="\t", row.names=1, check.names=F, stringsAsFactors=F, comment.char="")
inputdata <- inputdata[,rownames(SampleInfo)]
str(inputdata)

inputdata1 <- as.matrix(t(inputdata))
str(inputdata1)

#perplexity should not be bigger than 3 * perplexity < nrow(X) - 1
if (nrow(inputdata1) > 3){
	pd = (nrow(inputdata1)-1)%/%3

}else{
	pd = (nrow(inputdata1)-1)/3
}
data <- Rtsne(inputdata1, dims = 2, perplexity = pd, verbose = TRUE, max_iter = 500)

# ====================== Plot =======================
nm = max(nchar(Group))
height = 8
width = height + 0.1*nm

x = data$Y[,1]
y = data$Y[,2]

sites = as.data.frame(cbind(x, y))
SampleInfo1 = SampleInfo[order(rownames(SampleInfo)),]
sites$sample = SampleInfo1$Sample
rownames(sites) = sites$sample
sites = sites[,-3]

colnames(sites) = c("tsne1", "tsne2")
write.csv(sites, paste(output,"tsne.sites.csv",sep = "/"), quote=F) 

pdf(paste(output,"tsne.pdf",sep = "/"),width = width, height = height)

dataforplot = as.data.frame(sites)
dataforplot$group = SampleInfo1$Group

p = ggplot(dataforplot, aes(tsne1,tsne2)) + labs(x = paste("tsne1 "), y = paste("tsne2 "), title = "TSNE") + 
			geom_point(aes(tsne1,tsne2, color = as.character(dataforplot$group)), size = 5) + theme_bw() +
			theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title=element_blank())

print(p)
dev.off()

