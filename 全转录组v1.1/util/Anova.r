# Differential Analysis Using ANOVA
rm(list=ls())
args = commandArgs(trailingOnly=T)

if(length(args) != 3){

	cat(
    "Usage: 
    Rscript Anova.r countfile groupfile Threshold outputpath
        parameters ->
            countfile: [file -> always *_count_matrix.csv with row gene and column samples];
            groupfile: [file -> always sample.groups without header and only two column];
            outputpath: [path -> path for output]; \n")
	options("show.error.messages" = F)
	stop()
	
}

library(DESeq2)
library(gplots)

countfile  = normalizePath(args[1])
groupfile  = normalizePath(args[2])
outputpath = normalizePath(args[3])

countdata = read.table(countfile, header = T, sep = "\t", row.names = 1, check.names = F, comment.char = "", quote = "", fill = T)

SampleInfo = read.table(groupfile, header = F, sep="\t")
rownames(SampleInfo) = SampleInfo[,1]
colnames(SampleInfo) = c("sample","Group")

# sort by groups
SampleInfo = SampleInfo[order(SampleInfo$Group),]
Groups = levels(SampleInfo$Group)

# prepare data
count = countdata[,rownames(SampleInfo)]
count = count[rowSums(count) != 0 , ]
# re.count = apply(count, 2, function(d) d/sum(as.numeric(d)))
# #筛选相对丰度大于0.001的gene
# index = which(rowSums(re.count) > 0.001)

#标准化
dds <- DESeqDataSetFromMatrix(countData = count, colData = SampleInfo, design = ~ Group)
dds <- DESeq(dds,betaPrior=FALSE)
cnt <- as.data.frame(counts(dds, normalized=TRUE))
sample_num <- ncol(cnt)

DataForDEA = as.matrix(cnt)

res = sapply(1:nrow(DataForDEA), function(g) {
	mod = lm(DataForDEA[g,] ~ SampleInfo$Group)
	am = anova(mod)
	F1 = am[4][1, 1]   # F value
	P1 = am[5][1, 1]   # P value
	FP = cbind(F1, P1)
})

res = t(res)
rownames(res) = rownames(DataForDEA)
colnames(res) = c("Fvalue", "Pvalue")	

fdr = p.adjust(res[,2], method = "BH")
means = sapply(Groups, function(g) apply(DataForDEA[,as.character(SampleInfo[as.character(SampleInfo$Group) %in% g, 1])], 1, mean))
colnames(means) = paste("Mean.In.", Groups, sep = "")
	
result  = data.frame(Gene = rownames(DataForDEA), DataForDEA, means, res, FDR = fdr ,Type = 0,check.names = F)
result  = result[order(result$Pvalue),]
p_index = which(result$Pvalue < 0.05)
result[1:length(p_index),]$Type = "DEG";
result[(length(p_index) + 1):nrow(result),]$Type = "Not DEG";

anova_file = paste(outputpath, "/gene_ANOVA_Test_Result.xls", sep="")
write.table(result, anova_file, row.names = F, quote = F,sep = "\t")

#heatmap
myheatcol    = colorpanel(75, 'green','black','red')
data         = data.matrix(result[result$Type != "Not DEG", 2:sample_num + 1])
heatmap_file = paste(outputpath, "/heatmap.pdf", sep="")
pdf(heatmap_file , width = 12, height = 12)
heatmap.2( data,
           dendrogram   = 'row', 
           Colv         = FALSE,
           col          = myheatcol, 
           scale        = "row", 
           margins      = c(5,10),
           density.info = "none",
           trace        = "none",
           key          = TRUE,
           keysize      = 0.6,
           cexRow       = 0.01,
           cexCol       = 0.8,
           srtCol       = 90
)
dev.off()

#Top50
res = result[result$Type != 'Not DEG', ]
if( nrow(res) < 50 ){

    data = data.matrix(res[ , 2:sample_num + 1])          # 小于50,全画

}else{

    res  = res[order(res$Pvalue), ]
    data = data.matrix(res[1:50, 2:sample_num + 1])       #  相对丰度最高的50个

}

heatmap_file = paste(outputpath, "/Top50.heatmap.pdf", sep="")
pdf(heatmap_file , width = 12, height = 12)
heatmap.2( data,
           dendrogram   = 'both', 
           #Colv         = FALSE,
           col          = myheatcol, 
           scale        = "row", 
           margins      = c(5,10),
           density.info = "none",
           trace        = "none",
           key          = TRUE,
           keysize      = 0.6,
           cexRow       = 0.8,
           cexCol       = 0.8,
           srtCol       = 90
)
dev.off()
