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

#library(ade4)
library(ggrepel)

SampleInfo = read.table(groupfile, header = F, sep="\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]
Group = as.character(SampleInfo[,2])

inputdata <- read.table(input, header=T, sep="\t", row.names=1, check.names=F, stringsAsFactors=F, comment.char="")
inputdata <- inputdata[,rownames(SampleInfo)]
str(inputdata)
inputdata1 <- as.data.frame(t(inputdata))
str(inputdata1)

data <- prcomp(inputdata1)

# ====================== Plot =======================
nm = max(nchar(Group))
height = 8
width = height + 0.1*nm
# set color by group
# cols = c(119,132,147,454,89,404,123,463,552,28,54,84,100,558,43,31,610,477,256,588,99,81,503,104,562,76,96,495,570,616)
# gcol= colors()[rep(cols,20)][1:length(unique(Group))] 

x = data$x[,1]
y = data$x[,2]
sites = cbind(x, y)
colnames(sites) = c("PC1", "PC2")
write.csv(sites, paste(output,"pc.sites.csv",sep = "/"), quote=F)

PC1.proportion = summary(data)$importance[2,1]*100      
PC2.proportion = summary(data)$importance[2,2]*100  

pdf(paste(output,"pca.pdf",sep = "/"),width = width, height = height)

dataforplot = as.data.frame(sites[order(rownames(sites)),])
SampleInfo1 = SampleInfo[order(rownames(SampleInfo)),]
dataforplot$group = SampleInfo1$Group

p = ggplot(dataforplot, aes(PC1,PC2)) + labs(x = paste("PC1 ",sprintf("%.2f%%",PC1.proportion),"of var",sep = ""), y = paste("PC2 ",sprintf("%.2f%%",PC2.proportion),"of var",sep = ""), title = "PCA") + 
			geom_point(aes(PC1, PC2, color = as.character(dataforplot$group)), size = 5) + theme_bw() +
			theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title=element_blank()) 
print(p)
# par(oma=c(5,5,2,4))
# s.class(dfxy = sites,cgrid = 0,fac = as.factor(Group), col = gcol, xax = 1, yax = 2 , cellipse=1 )
# axis(side = 1 ,line =5 )
# axis(side = 2 ,line =4 )
# title(main = "PCA", xlab = paste("PC1 ",sprintf("%.2f%%",PC1.proportion),"of var",sep = ""), ylab = paste("PC2 ",sprintf("%.2f%%",PC2.proportion),"of var",sep = ""), outer = TRUE)

dev.off()

