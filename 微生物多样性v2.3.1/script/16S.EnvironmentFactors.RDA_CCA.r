# ��Ⱥ�뻷������֮�������
# http://blog.sina.com.cn/s/blog_83f77c940102vyqb.html

rm(list=ls())

args = commandArgs(trailingOnly=T)

if(length(args) != 4){

	cat(
	" ============== RDA_CCA ============
	Usage: 
	Rscript 16S.EnvironmentFactors.RDA_CCA.r envfile taxonfile groupfile outputpath
		parameters ->
			envfile: [file -> always env.xls with row env.xls and colum environmental factors];
		      taxonfile: [file -> always *.taxon.Abundance.xls];
		      groupfile: [file -> always sample.groups without header and only two column];
		     outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F) 
	stop()
	
}

library(vegan)
library(maptools)

envfile = normalizePath(args[1])
taxonfile = normalizePath(args[2])
groupfile = normalizePath(args[3])
outputpath = normalizePath(args[4])

# envfile = normalizePath("/home/raoja/project/16B0928B/report/env2.txt")
# taxonfile = normalizePath("/home/raoja/project/16B0928B/report/Community/Community_Structure/genus/genus.taxon.Abundance.xls")
# groupfile = normalizePath("/home/raoja/project/16B0928B/report/sample.groups")
# outputpath = normalizePath("/home/raoja/project/16B0928B/report/EnvironmentFactors/RDA_CCA/genus")

envdata = read.table(envfile, header = T, sep = "\t", row.names = 1)
taxondata = read.table(taxonfile, header = T, sep = "\t", row.names = 1, check.name = F, comment.char = "", quote= "", fill = T)
taxondata = taxondata[, rownames(envdata)]

loc = which(rownames(taxondata) == "All")       # ȥ�� All ��
if(length(loc) > 0) taxondata = taxondata[-loc, ]

# ��ȴ���1%������,���ڻ�ͼ
taxon = cbind(taxondata, apply(taxondata, 1, sum))
taxon1 = setdiff(rownames(taxon)[which(taxon[, ncol(taxon)]/sum(taxon[, ncol(taxon)]) > 0.01)], "No_Rank")

loc = which(rownames(taxondata) == "No_Rank")       # ȥ�� No_Rank ��
if(length(loc) > 0) taxondata = taxondata[-loc, ]

taxondata = taxondata[rowSums(taxondata) != 0 , ]

# prepare data
rownames(taxondata) = gsub("\"","",rownames(taxondata))     # ȥ�������е�˫����,ȫ���滻
colnames(taxondata) = gsub("\"","",colnames(taxondata))
rownames(taxondata) = as.character(strsplit(rownames(taxondata),split="{.*}",perl=T))

SampleInfo = read.table(groupfile, header = F, sep = "\t")
colnames(SampleInfo) = c("Sample", "Group")
rownames(SampleInfo) = SampleInfo[,1]

taxondata = t(taxondata)

source("/home/zhengy/bin/modules/script/decorana_new.r")
otu.dca = decorana_new(taxondata)
capture.output(otu.dca, file = paste(outputpath, "DCA.csv", sep = "/")) 

# �������ݶȳ���4����ѡ�񵥷�ģ��CCA�����ʣ�������ݶ�С��3����ѡ������ģ��RDA�ȽϺ����������3-4֮�䣬�򵥷�ģ�ͺ�����ģ�Ͷ����á�
# http://blog.sina.com.cn/s/blog_670445240101npsa.html
Axis.lengths = capture.output(otu.dca)[11]
maxlen = as.numeric(unlist(strsplit(Axis.lengths," +"))[3])

if( maxlen > 3 ){       # ѡ�� CCA

	result = cca(taxondata, envdata, scale = T)
	#capture.output(result, file = "CCA.csv")	
	select = "CCA"

}else{        # ѡ�� RDA

	result = rda(taxondata, envdata, scale = T)
	#capture.output(result, file = "RDA.csv")
	select = "RDA"

}

# envfit ���黷��������������������	
result.fit = envfit(result, envdata, perm = 999)   	# envfit ���黷��������������������	
capture.output(result.fit, file = paste(outputpath, "envfit.txt", sep = "/"))

Axis1 = summary(result)$concont$importance[2, 1]*100
Axis2 = summary(result)$concont$importance[2, 2]*100

# set color and shape by group
mycol = rep(c("#769C30","#7B0078","#3A89CC","#D99536","#BE0027","#BFBC3B","#4C8B35","#3C68AE","#C10077","#CAAA76","#2E165B","#458B00","#8B4513","#008B8B","#6E8B3D","#8B7D6B","#7FFF00","#FF69B4"),time =2)
mypch = c(21:25,3,4,7,9,8,10,15:18,0:14)

index = sapply(1:nrow(SampleInfo),function(x) which(unique(SampleInfo$Group) %in% SampleInfo[x,2]))
pcol = mycol[index]     
ppch = mypch[index]             

# ========== �������������ӡ����� ========== 

# Samples
sites = summary(result)$sites

if( nrow(sites) < 10 ){

	pdf(paste(outputpath, paste(select, ".pdf", sep = ""), sep = "/"), width = 9, height = 9)
	par(mar = c(3.3, 3.3, 0.5, 0.5), mgp = c(1.9, 0.5, 0))
	plot(sites[,1:2], pch = ppch, col = pcol, bg = pcol, xlab = paste(select, "1  [", round(Axis1, 2), "%]", sep=""), ylab = paste(select, "2  [", round(Axis2, 2), "%]", sep=""), cex.lab = 1.2)
	
	pointLabel(x = sites[,1], y = sites[,2], labels = rownames(sites), cex = 1.1, col = pcol)

}else{

	nm = 0.2*max(nchar(as.character(SampleInfo[,2])))
	nm = ifelse(nm>1, nm, 1)
	pdf(paste(outputpath, paste(select, ".pdf", sep = ""), sep = "/"), width = 9 + nm, height = 9)
	layout(t(matrix(c(1,2))), widths = c(9, nm), respect = FALSE)    # ����
	par(mar = c(3.3, 3.3, 0.5, 0), mgp = c(1.9, 0.5, 0))
	plot(sites[,1:2], pch = ppch, col = pcol, bg = pcol, xlab = paste(select, "1  [", round(Axis1, 2), "%]", sep=""), ylab = paste(select, "2  [", round(Axis2, 2), "%]", sep=""), cex.lab = 1.2)

}

abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

# EnvironmentFactors
ef = summary(result)$biplot
zoom = min(c(abs(max(sites[,1])/max(ef[,1])), abs(min(sites[,1])/min(ef[,1])), abs(max(sites[,2])/max(ef[,2])), abs(min(sites[,2])/min(ef[,2]))))     # �������ӷ�������
if(zoom > 1){
	ef = summary(result)$biplot*floor(zoom)
}else{
	ef = summary(result)$biplot
}
arrows(0, 0, ef[,1], ef[,2], col = "red", angle = 10, length = 0.1)
pointLabel(x = ef[,1], y = ef[,2], labels = rownames(ef), cex = 1, col = "red")

# Taxon
spe = summary(result)$species
if(nrow(spe) < 10){
	spe = spe
}else{
	spe = spe[taxon1,]
}
points(x = spe[,1], y = spe[,2], pch = 8, cex = 0.6, col = "black")
pointLabel(x = spe[,1], y = spe[,2], labels = rownames(spe), cex = 0.7, col = "black")

if( nrow(sites) >= 10 ){

	par(mar = c(0,0,0,0), xpd = TRUE)
	plot.new()
	legend("center", pch = unique(ppch), col = unique(pcol), pt.bg = unique(pcol), legend = unique(SampleInfo$Group), bty = "n")
	
}

dev.off()


# ========== �������������� ========== 

# Samples
sites = summary(result)$sites

if( nrow(sites) < 10 ){

	pdf(paste(outputpath, paste(select, ".WithoutTaxon.pdf", sep = ""), sep = "/"), width = 9, height = 9)
	par(mar = c(3.3, 3.3, 0.5, 0.5), mgp = c(1.9, 0.5, 0))
	plot(sites[,1:2], pch = ppch, col = pcol, bg = pcol, xlab = paste(select, "1  [", round(Axis1, 2), "%]", sep=""), ylab = paste(select, "2  [", round(Axis2, 2), "%]", sep=""), cex.lab = 1.2)
	
	pointLabel(x = sites[,1], y = sites[,2], labels = rownames(sites), cex = 1.1, col = pcol)

}else{

	nm = 0.2*max(nchar(as.character(SampleInfo[,2])))
	nm = ifelse(nm>1, nm, 1)
	pdf(paste(outputpath, paste(select, ".WithoutTaxon.pdf", sep = ""), sep = "/"), width = 9 + nm, height = 9)
	layout(t(matrix(c(1,2))), widths = c(9, nm), respect = FALSE)    # ����
	par(mar = c(3.3, 3.3, 0.5, 0), mgp = c(1.9, 0.5, 0))
	plot(sites[,1:2], pch = ppch, col = pcol, bg = pcol, xlab = paste(select, "1  [", round(Axis1, 2), "%]", sep=""), ylab = paste(select, "2  [", round(Axis2, 2), "%]", sep=""), cex.lab = 1.2)

}

abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

# EnvironmentFactors
ef = summary(result)$biplot
zoom = min(c(abs(max(sites[,1])/max(ef[,1])), abs(min(sites[,1])/min(ef[,1])), abs(max(sites[,2])/max(ef[,2])), abs(min(sites[,2])/min(ef[,2]))))     # �������ӷ�������
if(zoom > 1){
	ef = summary(result)$biplot*floor(zoom)
}else{
	ef = summary(result)$biplot
}
arrows(0, 0, ef[,1], ef[,2], col = "red", angle = 10, length = 0.1)
pointLabel(x = ef[,1], y = ef[,2], labels = rownames(ef), cex = 1, col = "red")

if( nrow(sites) >= 10 ){

	par(mar = c(0,0,0,0), xpd = TRUE)
	plot.new()
	legend("center", pch = unique(ppch), col = unique(pcol), pt.bg = unique(pcol), legend = unique(SampleInfo$Group), bty = "n")
	
}

dev.off()

write.csv(summary(result)$concont$importance, paste(outputpath, "concont.csv", sep = "/"))       # ͳ����Ϣ
write.csv(sites[,1:2], paste(outputpath, "sample.position.csv", sep = "/"))        # ��������
write.csv(spe[,1:2], paste(outputpath, "species.position.csv", sep = "/"))            # ��������
write.csv(ef[,1:2], paste(outputpath, "environmentfactors.position.csv", sep = "/"))         # ������������
