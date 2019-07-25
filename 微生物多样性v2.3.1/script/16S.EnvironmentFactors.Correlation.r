# correlation between EnvironmentFactors and taxon 

rm(list=ls())

args = commandArgs(trailingOnly=T)
if(length(args) != 4){

	cat(
	" ============== correlation between EnvironmentFactors and taxon =========
	Usage: 
	Rscript 16S.EnvironmentFactors.Correlation.r envfile taxonfile N outputpath
		parameters ->
			  envfile: [file -> always env.xls with row env.xls and colum environmental factors];
			taxonfile: [file -> always *.taxon.size.xls with row taxon and column samples];
				N: [number -> N most abundant taxa to plot heatmap]
		       outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

library(vegan)
library(Hmisc)
library(pheatmap)

# An Introduction to corrplot package： http://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

envfile = normalizePath(args[1])
taxonfile = normalizePath(args[2])
N = as.numeric(args[3])
outputpath = normalizePath(args[4])

envdata = read.table(envfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "")
taxondata = read.table(taxonfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
taxondata = taxondata[, rownames(envdata)]

loc = which(rownames(taxondata) == "All")       # 去除 All 行
if(length(loc) > 0) taxondata = taxondata[-loc, ]

loc = which(rownames(taxondata) == "No_Rank")       # 去除 No_Rank 行
if(length(loc) > 0) taxondata = taxondata[-loc, ]
taxondata = taxondata[rowSums(taxondata)!=0,]

# prepare data
rownames(taxondata) = gsub("\"","",rownames(taxondata))     # 去除名字中的双引号,全局替换
colnames(taxondata) = gsub("\"","",colnames(taxondata))
rownames(taxondata) = as.character(strsplit(rownames(taxondata),split="{.*}",perl=T))

envdata = t(envdata)

# correlation between two matrix
# cor_r = t(sapply(1:nrow(taxondata), function(x) {sapply(1:nrow(envdata), function(y) rcorr(as.matrix(taxondata)[x,],envdata[y,])[[1]][1,2])}))
cor_r = t(sapply(1:nrow(taxondata), function(x) {sapply(1:nrow(envdata), function(y) cor.test(as.matrix(taxondata)[x,],envdata[y,], method = "spearman")[[4]])}))

rownames(cor_r) = rownames(taxondata)
colnames(cor_r) = rownames(envdata)
write.csv(cor_r, paste(outputpath, "cor_r.csv", sep="/"))
# cor_p = t(sapply(1:nrow(taxondata), function(x) {sapply(1:nrow(envdata), function(y) rcorr(as.matrix(taxondata)[x,],envdata[y,])[[3]][1,2])}))
cor_p = t(sapply(1:nrow(taxondata), function(x) {sapply(1:nrow(envdata), function(y) cor.test(as.matrix(taxondata)[x,],envdata[y,],method = "spearman")[[3]])}))

rownames(cor_p) = rownames(taxondata)
colnames(cor_p) = rownames(envdata)
write.csv(cor_p, paste(outputpath, "cor_p.csv", sep="/"))

# 提取相关系数和P值
Corresult = data.frame(Taxon = rep(rownames(cor_r), ncol(cor_r)), Factor = rep(colnames(cor_r), each = nrow(cor_r)), Cor = as.vector(cor_r), Pvalue = as.vector(cor_p))
write.csv(Corresult, paste(outputpath, "Correlation_All.csv", sep = "/"), row.names = F)

# 显著
sigres = Corresult[which(abs(Corresult$Cor) > 0.3 & Corresult$Pvalue < 0.05),]
write.csv(sigres, paste(outputpath, "Correlation_cor0.3_p0.05.csv", sep = "/"), row.names = F)


if( nrow(cor_r) <= N ){

	DataForHeatmap = cor_r	   # 小于N,全画
	display = cor_p


}else{

	# get relative.abundance
	taxondata = cbind(taxondata, apply(taxondata, 1, sum))
	taxondata = taxondata[order(taxondata[,ncol(taxondata)], decreasing = T), ]
	DataForHeatmap = cor_r[rownames(taxondata)[1:N],]           #  相对丰度最高的N个
	display = cor_p[rownames(taxondata)[1:N],]
}

# 用于将显著标记在方块上
display[which(display < 0.01)] = "**"
display[which(display > 0.01 & display < 0.05)] = "*"
display[which(display >= 0.05)] = ""

malen = max(nchar(rownames(DataForHeatmap)))
malen2 = max(nchar(colnames(DataForHeatmap)))

width = 0.15*ncol(DataForHeatmap) + 0.1*malen     # 宽度由样本数和分类名长度决定
height = 0.15*nrow(DataForHeatmap) + 0.1*malen2     # 高度由分类数和样本名长度决定

if(width < 4) width = 4
if(height < 4) height = 4
if(height > 20) height = 20

#col <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))(200)
name = unlist(strsplit(taxonfile,"/"))
if(nrow(DataForHeatmap) == 1){
	pdf(paste(outputpath, paste(unlist(strsplit(name[length(name)], "\\."))[1], ".taxon.Correlation.withEnvironmentFactors.pdf", sep = ""), sep="/"), width = 2*width, height = 3.5, onefile = FALSE)  # onefile = FALSE 去除一页空白 
	pheatmap(DataForHeatmap, cluster_rows = F, display_numbers = display)
}else{
	pdf(paste(outputpath, paste(unlist(strsplit(name[length(name)], "\\."))[1], ".taxon.Correlation.withEnvironmentFactors.pdf", sep = ""), sep="/"), width = width, height = height, onefile = FALSE) 
	pheatmap(DataForHeatmap, fontsize_row = 8, fontsize_col = 9, display_numbers = display)
}
dev.off()


#taxon.Correlation_cor0.3_p0.05.withEnvironmentFactors.pdf
DataForHeatmap = cor_r[rownames(cor_r) %in% sigres$Taxon,] 
display = cor_p[rownames(cor_p) %in% sigres$Taxon,]
display[which(display < 0.01)] = "**"
display[which(display > 0.01 & display < 0.05)] = "*"
display[which(display >= 0.05)] = ""

if(length(sigres$Taxon) == 1){
	DataForHeatmap = t(as.matrix(DataForHeatmap))
	rownames(DataForHeatmap) = sigres$Taxon
	colnames(DataForHeatmap) = rownames(envdata)
	display = t(as.matrix(display))
	rownames(display) = sigres$Taxon
	colnames(display) = rownames(envdata)

	pdf(paste(outputpath, paste(unlist(strsplit(name[length(name)], "\\."))[1], ".taxon.Correlation_cor0.3_p0.05.withEnvironmentFactors.pdf", sep = ""), sep="/"), width = 2*width, height = 3.5, onefile = FALSE)  # onefile = FALSE 去除一页空白 
	pheatmap(DataForHeatmap, cluster_rows = F, display_numbers = display)
	dev.off()
}
if(nrow(DataForHeatmap) > 1){
	pdf(paste(outputpath, paste(unlist(strsplit(name[length(name)], "\\."))[1], ".taxon.Correlation_cor0.3_p0.05.withEnvironmentFactors.pdf", sep = ""), sep="/"), width = width, height = height, onefile = FALSE) 
	pheatmap(DataForHeatmap, fontsize_row = 8, fontsize_col = 9, display_numbers = display)
	dev.off()
}

