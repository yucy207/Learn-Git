rm(list = ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 2){

	cat(
	"============== WGCNA ===========
	Usage: 
	Rscript WGCNA.R fpkmfile outputpath
		parameters ->
                fpkmfile: [file -> always genes.fpkm.xls];
                outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()

}

fpkmfile = normalizePath(args[1])
outputpath = normalizePath(args[2])
setwd(outputpath)

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
#读取表达量文件
fpkmdata = read.table(fpkmfile, header = T, sep = "\t", row.names = 1, check.names = F, quote = "", fill = T, comment.char = "")
#转置，行名为样本，列名为基因ID
datExpr0 = as.data.frame(t(fpkmdata))
#数据过滤
gsg = goodSamplesGenes(datExpr0)
if (!gsg$allOK) {datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]}

#Choosing the soft-thresholding power: 乘方值（power）选择
powers = c(1:30)
sft = pickSoftThreshold(datExpr0, powerVector = powers,  RsquaredCut = 0.85, verbose = 5)
softPower = sft$powerEstimate

# 绘制power图
pdf("scale.independence.pdf")
plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  main = paste("Scale independence")
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,
  cex=0.9,
  col="red"
)

abline(h=0.85, col="red")
dev.off()

pdf("mean.connectivity.pdf")
plot(
sft$fitIndices[,1],
sft$fitIndices[,5],
xlab="Soft Threshold (power)",
ylab="Mean Connectivity",
type="n",
main = paste("Mean connectivity")
)

text(
sft$fitIndices[,1],
sft$fitIndices[,5],
labels=powers,
cex=0.9,,
col="red"
)

dev.off()



net = blockwiseModules(datExpr0, power = softPower,
						TOMType = "unsigned", minModuleSize = 30,
						reassignThreshold = 0, mergeCutHeight = 0.25,
						numericLabels = TRUE, pamRespectsDendro = FALSE,
						#saveTOMs = TRUE,
						#saveTOMFileBase = "femaleMouseTOM",
						verbose = 3)

geneTree = net$dendrograms
moduleColors = labels2colors(net$colors)

# 绘制基因聚类的热图
pdf("genes.cluster.pdf")
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

# 绘制modules 聚类图
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MET = orderMEs(MEs0)


pdf("modules.cluster.pdf")
par(mar = c(4, 2, 1, 2), cex = 0.9)
plotEigengeneNetworks(
MET, "",
plotHeatmaps = FALSE,
marDendro = c(0,4,1,2),
marHeatmap = c(3,4,1,2),
cex.lab = 0.8,
xLabelsAngle = 90
)
dev.off()

# 绘制modules 之间相关性的热图

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(
MEs,
MEs,
use = "p"
)

moduleTraitPvalue = corPvalueStudent(
moduleTraitCor,
nSamples
)



textMatrix =  paste(
signif(moduleTraitCor, 2),
"\n(",
signif(moduleTraitPvalue, 1),
")",
sep = ""
)

dim(textMatrix) = dim(moduleTraitCor)

pdf("module-module.relationship.pdf", width = 15, height = 7)
par(mar = c(4, 8, 4, 4))
labeledHeatmap(
Matrix = moduleTraitCor,
xLabels = names(MEs),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-Module relationships")
)

dev.off()



#绘制距离矩阵图热图
dissTOM = 1 - TOMsimilarityFromExpr(datExpr0, power = softPower)
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
png("Network_heatmap.png")
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

#提取输出每个module和对应的gene
data = data.frame(gene_id = names(datExpr0), module_id = (net$colors + 1))
Data = data[order(data$module_id),]
write.table(Data, "module_gene.xls", sep = "\t", quote = F, row.names = F)

#筛选关键基因
datKME = signedKME(datExpr0, net$MEs, outputColumnName="MM.")
colnames(datKME) = gsub("MM.", "", colnames(datKME))
colnames(datKME) = as.numeric(colnames(datKME)) + 1

filtergene = lapply(unique(Data$module_id), function(t) {

	gene_all = Data$gene_id[Data$module_id == t]
	gene_filter = rownames(datKME)[abs(datKME[,which(colnames(datKME) == t)]) > 0.8]
	intersection = intersect(gene_all,gene_filter)

})

filterdata = data.frame(gene_id = unlist(filtergene), 
						module_id = unlist(lapply(1:ncol(datKME), function(t) {rep(t, times = length(filtergene[[t]])) }))
)

write.table(filterdata, "module_filter_gene.xls", sep = "\t", quote = F, row.names = F)

TOM = TOMsimilarityFromExpr(datExpr0, power = softPower)
probes = names(datExpr0)


for(i in unique(moduleColors)) {

  inModule = is.finite(match(moduleColors, i));
  modProbes = probes[inModule];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)

  cyt = exportNetworkToCytoscape(
   modTOM,
   edgeFile = paste(i, ".edges.txt",sep=""),
   nodeFile = paste(i, ".nodes.txt",  sep=""),
   weighted = TRUE,
   threshold = 0.02,
   nodeNames = modProbes,
   nodeAttr = moduleColors[inModule]
 )

  
}