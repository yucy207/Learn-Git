rm(list=ls())

args = commandArgs(trailingOnly=T) 

if(length(args) != 3){

	cat(
	"============== Samples ClusterTree ===========
	Usage: 
	Rscript 16S.BetaDiversity.SamplesTree.r otufile groupfile outputpath
		parameters ->
			  otufile: [file -> always otu.tax.0.03.xls];
			groupfile: [file -> always sample.groups without header and only two column];
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

library(vegan)

otufile = normalizePath(args[1])
groupfile = normalizePath(args[2])
outputpath = normalizePath(args[3])

# otufile = normalizePath("/home/raoja/project/metagenomic/16B0316B/30sampleV4V5/report/OTU/subsample_otu.tax.0.03.xls")
# groupfile = normalizePath("/home/raoja/project/metagenomic/16B0316B/30sampleV4V5/report/sample.groups")
# outputpath = normalizePath("/home/panrf")

# ====================== prepare data =======================
# sampleinfo
sampleinfo = read.table(groupfile, header = F, sep = "\t")
samples = as.character(sampleinfo[,1])
groups = as.character(sampleinfo[,2])
Group = unique(groups)

OTUFile = read.table(otufile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
otu_data = OTUFile[,samples]

loc = which(rownames(otu_data) == "All")       # 去除 All 行
if(length(loc) > 0) otu_data = otu_data[-loc, ]

# for plot
height = 10
width = height + 0.1*ncol(otu_data)

# ggplot default color set
#gg_color_hue <- function(n) {
	#hues = seq(15, 375, length = n + 1)
	#hcl(h = hues, l = 65, c = 100)[1:n]
#}

n = length(unique(groups))
#cols = gg_color_hue(n)

mycol = c(552,81,53,460,476,551,258,30,426,430,564,561,619,84,576,521,644,10,46,51,53,128,101,107,371,373,386,435,429,528,512,375)
mycol = colors()[rep(mycol,20)]
cols = mycol[1:n]

# set color for groups
ordergroup = sort(levels(factor(groups)))
ordercolor = rep(1, length(groups))

for(i in 1:length(ordergroup)){
	ordercolor[groups%in%ordergroup[i]] = cols[i]
}

fga.dist = vegdist(t(otu_data), method = "bray")		# BC distance
hc = hclust(fga.dist, method = "average")		# UPGMA

y <- rep(hc$height, 2)
x <- as.numeric(hc$merge)
y <- y[which(x < 0)]
x <- x[which(x < 0)]
x <- abs(x)
y <- y[order(x)]
x <- x[order(x)]

# 叶子不齐的树状图
# pdf(paste(outputpath, "SamplesClusterTree.pdf", sep = "/"), width = width, height = height)
# layout(matrix(c(1,2)), heights = c(8, 2), respect = FALSE)
# par(mar = c(0, 4.5, 2, 1))
# plot(hc, labels = FALSE, hang = 0.1, xlab = "", sub = "")
# text(x = x, y = y[hc$order] - (max(hc$height) * 0.1), labels = samples[hc$order], col = ordercolor[hc$order],
        # srt = 90, adj = c(1, 0.5), xpd = NA, cex = 0.7)
# plot.new()
# par(mar = c(0, 0, 0, 0))
# legend("center", legend = ordergroup, fill = cols, ncol = ceiling(width*2/max(nchar(groups))), border = NA, bty = "n")
# dev.off()


# 叶子齐平的树状图
hcd = as.dendrogram(hc)

ncol = ceiling(width*2/max(nchar(groups)))
if( n < ncol ) ncol = n  # 组数小于列数，则根据组数设定列数
nrow = length(Group)/ncol

pdf(paste(outputpath, "SamplesClusterTree.pdf", sep = "/"), width = width, height = height)

kuan = 1.1*nrow
if (kuan > 5) kuan = 5

layout(matrix(c(1,2)), heights = c(8, kuan), respect = FALSE)

#set text size according to the length of the sample names
length = max(nchar(colnames(otu_data)))

par(mar = c(0.4*length, 4.5, 1, 1))
plot(hcd, xlab = "", sub = "", leaflab = "none")

if(length < 3 ){
	text(x = x, y = -length/100, labels = samples[hc$order], col = ordercolor[hc$order], xpd = NA, srt =90, adj = c(1, 0.5), cex = 1.5)
}else{
	text(x = x, y = -length/400, labels = samples[hc$order], col = ordercolor[hc$order], srt =90, adj = c(1, 0.5), xpd = NA, cex = 1.3)
}

par(mar = c(0, 0, 2, 0))
plot.new()

legend("center",legend = ordergroup, fill = cols, ncol = ncol, border = NA, bty = "n", cex = 1.5)

dev.off()
